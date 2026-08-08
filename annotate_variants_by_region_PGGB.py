#!/usr/bin/env python3
"""
annotate_variants_by_region.py

For each VCF file (PGGB or per-reference), intersects variants with the
corresponding GFF3 file using bedtools and generates a
"Variant distribution by genomic region" table.

Regions: CDS, exon, inferred intron, UTR5, UTR3, and intergenic
(genomic complement).

Annotation hierarchy (a variant may overlap multiple features):
  CDS > UTR5/UTR3 > exon > intron > intergenic

The most specific region is assigned.

Usage:
    python annotate_variants_by_region.py

Edit the CONFIGURATION section.
"""

import os
import sys
import glob
import subprocess
import pandas as pd
from datetime import datetime

# =========================
# CONFIGURATION FROM ENVIRONMENT VARIABLES
# (defined in the SGE job .sh file)
# =========================

GFF3_DIR     = os.environ.get("GFF3_DIR")
FAI_DIR      = os.environ.get("FAI_DIR", os.environ.get("GFF3_DIR", ""))
PGGB_BASE    = os.environ.get("PGGB_BASE")
GRAPH_HASH   = os.environ.get("GRAPH_HASH", "bf3285f")
PER_REF_BASE = os.environ.get("PER_REF_BASE")
OUTPUT_DIR   = os.environ.get("OUTPUT_DIR")

# COMMUNITY may contain a single community passed by the array job
# or a comma-separated list of communities
_comm_env   = os.environ.get("COMMUNITY", "")
COMMUNITIES = [c.strip() for c in _comm_env.split(",") if c.strip()]

# Validate required variables
for _var, _val in [
    ("GFF3_DIR", GFF3_DIR),
    ("PGGB_BASE", PGGB_BASE),
    ("PER_REF_BASE", PER_REF_BASE),
    ("OUTPUT_DIR", OUTPUT_DIR),
    ("COMMUNITY", _comm_env),
]:
    if not _val:
        print(
            f"[ERROR] Environment variable '{_var}' is not defined.",
            flush=True,
        )
        sys.exit(1)

# Region hierarchy
# Higher values indicate more specific regions and therefore take priority
REGION_PRIORITY = {
    "intergenic": 0,
    "intron":     1,
    "exon":       2,
    "UTR5":       3,
    "UTR3":       4,
    "CDS":        5,
}

# =========================
# SETUP
# =========================

os.makedirs(OUTPUT_DIR, exist_ok=True)

tmp_dir = os.path.join(OUTPUT_DIR, "tmp")
os.makedirs(tmp_dir, exist_ok=True)

start_time = datetime.now()

print("=" * 65, flush=True)
print(f"Start      : {start_time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
print(f"OUTPUT_DIR : {OUTPUT_DIR}", flush=True)
print("=" * 65, flush=True)

# =========================
# GFF3 FUNCTIONS
# =========================

def find_file_for_sample(sample, search_dir, extensions):
    """
    Find a file for the given sample by testing several extensions.

    Uses regular-expression boundaries to prevent sample names such as CO8
    from incorrectly matching names such as CO84.
    """
    import re

    pattern = re.compile(
        r"(?<![A-Za-z0-9])"
        + re.escape(sample)
        + r"(?![A-Za-z0-9])"
    )

    candidates = []

    for ext in extensions:
        candidates += glob.glob(os.path.join(search_dir, ext))

    candidates_sorted = sorted(
        candidates,
        key=lambda p: len(os.path.basename(p)),
        reverse=True,
    )

    for candidate in candidates_sorted:
        if pattern.search(os.path.basename(candidate)):
            return candidate

    return None


def find_gff3_for_sample(sample, gff3_dir):
    """Find the GFF3 file associated with a sample."""
    return find_file_for_sample(sample, gff3_dir, ["*.gff3"])


def find_fai_for_sample(sample, fai_dir):
    """
    Find the sample's .fai file in fai_dir.

    Returns:
        tuple: (fai_path, fasta_base_path)
    """
    fai = find_file_for_sample(sample, fai_dir, ["*.fai"])

    if fai:
        # The FASTA path is the .fai path without the final ".fai"
        fasta_base = fai[:-4]
        return fai, fasta_base

    return None, None


def normalize_chrom(chrom, sample):
    """
    Convert chromosome names from the GFF3 or FAI file to the PanSN format
    used in the VCF.

    Conversion examples:
      GFF3: MrorC26_Group1    → C26#0#Group1
      GFF3: MrorCO8_Group1    → CO8#0#Group1
      FAI:  Mror_C26_Group1   → C26#0#Group1
      FAI:  Mror_CO84_Group1  → CO84#0#Group1

    Strategy:
      1. Remove the "Mror_" or "Mror" prefix, with or without an underscore.
      2. Extract the sample name, which appears immediately after the prefix.
      3. Treat the sequence after the sample name and "_" as the contig name.
      4. Format the result as SAMPLE#0#CONTIG.
    """
    import re

    # Remove the Mror_ or Mror prefix, with or without an underscore
    cleaned = re.sub(r"^Mror_?", "", chrom)

    # At this point, cleaned should look like:
    # "C26_Group1" or "CO84_Group1"
    if "_" in cleaned:
        parts  = cleaned.split("_", 1)
        smp    = parts[0]  # C26, CO84, B3, etc.
        contig = parts[1]  # Group1, Group7, Ungrouped1, etc.

        return f"{smp}#0#{contig}"

    # Return the original name if the pattern is not recognized
    return chrom


def get_chrom_sizes_normalized(fai_path, sample, tmp_dir):
    """
    Read a .fai file and generate a chrom.sizes file with chromosome names
    normalized to the PanSN format.

    Returns the path to the normalized chrom.sizes file.
    """
    if not fai_path or not os.path.exists(fai_path):
        print(f"    [WARNING] .fai file not found for {sample}", flush=True)
        return None

    sizes_path = os.path.join(tmp_dir, f"{sample}_chrom.sizes")

    written = 0
    norm_name = None

    with open(fai_path) as fin, open(sizes_path, "w") as fout:
        for line in fin:
            cols = line.strip().split("\t")

            if len(cols) < 2:
                continue

            norm_name = normalize_chrom(cols[0], sample)
            fout.write(f"{norm_name}\t{cols[1]}\n")
            written += 1

    if written == 0:
        print(
            f"    [ERROR] chrom.sizes is empty for {sample} — check the .fai file",
            flush=True,
        )
        return None

    print(
        f"    chrom.sizes: {written} normalized chromosomes "
        f"(example: {norm_name})",
        flush=True,
    )

    return sizes_path


def gff3_to_region_beds(
    gff3_path,
    sample,
    tmp_dir,
    fasta_sizes_path=None,
):
    """
    Parse the GFF3 file and generate BED files for each genomic region:
      CDS, exon, UTR5, UTR3, and gene.

    Gene coordinates are used to infer intronic and intergenic regions.

    Introns:
      gene coordinates - exon coordinates

    Intergenic regions:
      complement of gene coordinates

    Returns:
        dict: {region_name: bed_path}
    """
    records = {
        "CDS":  [],
        "exon": [],
        "UTR5": [],
        "UTR3": [],
        "gene": [],
    }

    print(
        f"    Parsing GFF3: {os.path.basename(gff3_path)}",
        flush=True,
    )

    with open(gff3_path) as gff3_file:
        for line in gff3_file:
            if line.startswith("#") or not line.strip():
                continue

            cols = line.strip().split("\t")

            if len(cols) < 9:
                continue

            raw_chrom = cols[0]

            # Normalize chromosome names to the PanSN format:
            # SAMPLE#0#Contig
            chrom = normalize_chrom(raw_chrom, sample)

            feat  = cols[2]
            start = cols[3]
            end   = cols[4]

            # Convert GFF3 coordinates to 0-based BED coordinates
            try:
                s = int(start) - 1
                e = int(end)
            except ValueError:
                continue

            feat_lower = feat.lower()

            if feat_lower == "cds":
                records["CDS"].append((chrom, s, e))

            elif feat_lower == "exon":
                records["exon"].append((chrom, s, e))

            elif feat_lower in (
                "five_prime_utr",
                "5utr",
                "utr5",
                "five_prime_utr_region",
            ):
                records["UTR5"].append((chrom, s, e))

            elif feat_lower in (
                "three_prime_utr",
                "3utr",
                "utr3",
                "three_prime_utr_region",
            ):
                records["UTR3"].append((chrom, s, e))

            elif feat_lower == "gene":
                records["gene"].append((chrom, s, e))

    bed_paths = {}

    # Write the basic BED files and sort them
    for region, rows in records.items():
        if not rows:
            continue

        bed_path = os.path.join(
            tmp_dir,
            f"{sample}_{region}.bed",
        )

        with open(bed_path, "w") as bed_file:
            for chrom, start, end in sorted(rows):
                bed_file.write(f"{chrom}\t{start}\t{end}\n")

        # Merge overlapping regions
        merged_path = os.path.join(
            tmp_dir,
            f"{sample}_{region}_merged.bed",
        )

        subprocess.run(
            f"bedtools sort -i {bed_path} "
            f"| bedtools merge -i - > {merged_path}",
            shell=True,
            check=True,
        )

        bed_paths[region] = merged_path

    # Infer introns as gene coordinates minus exon coordinates
    if "gene" in bed_paths and "exon" in bed_paths:
        intron_path = os.path.join(
            tmp_dir,
            f"{sample}_intron_merged.bed",
        )

        subprocess.run(
            f"bedtools subtract "
            f"-a {bed_paths['gene']} "
            f"-b {bed_paths['exon']} "
            f"> {intron_path}",
            shell=True,
            check=True,
        )

        bed_paths["intron"] = intron_path

    elif "gene" in bed_paths:
        # If exon annotations are unavailable, treat the complete gene
        # coordinates as intronic
        bed_paths["intron"] = bed_paths["gene"]

    # Infer intergenic regions as the complement of genes in the genome.
    # This requires a chromosome-sizes file.
    if (
        "gene" in bed_paths
        and fasta_sizes_path
        and os.path.exists(fasta_sizes_path)
    ):
        intergenic_path = os.path.join(
            tmp_dir,
            f"{sample}_intergenic_merged.bed",
        )

        subprocess.run(
            f"bedtools complement "
            f"-i {bed_paths['gene']} "
            f"-g {fasta_sizes_path} "
            f"> {intergenic_path}",
            shell=True,
            check=True,
        )

        bed_paths["intergenic"] = intergenic_path

    else:
        print(
            "    [WARNING] No chromosome-sizes file was found — "
            "intergenic variants will be assigned based on the absence "
            "of gene overlap",
            flush=True,
        )

    # Remove "gene" from the final dictionary because it was only used
    # to infer intronic and intergenic regions
    bed_paths.pop("gene", None)

    print(
        f"    Generated regions: {list(bed_paths.keys())}",
        flush=True,
    )

    return bed_paths


def get_chrom_sizes(fasta_path, tmp_dir, sample):
    """
    Generate a chromosome-sizes file using samtools faidx or an existing
    FASTA index.
    """
    sizes_path = os.path.join(
        tmp_dir,
        f"{sample}_chrom.sizes",
    )

    if os.path.exists(sizes_path):
        return sizes_path

    fai_path = fasta_path + ".fai"

    if os.path.exists(fai_path):
        # Extract chromosome name and size from columns 1 and 2 of the .fai file
        subprocess.run(
            f"cut -f1,2 {fai_path} > {sizes_path}",
            shell=True,
            check=True,
        )
        return sizes_path

    # Attempt to generate the .fai file with samtools
    result = subprocess.run(
        f"samtools faidx {fasta_path}",
        shell=True,
        capture_output=True,
    )

    if result.returncode == 0 and os.path.exists(fai_path):
        subprocess.run(
            f"cut -f1,2 {fai_path} > {sizes_path}",
            shell=True,
            check=True,
        )
        return sizes_path

    print(
        f"    [WARNING] Could not generate chrom.sizes for {sample}",
        flush=True,
    )

    return None

# =========================
# INTERSECTION FUNCTIONS
# =========================

def intersect_vcf_with_regions(
    vcf_path,
    bed_paths,
    sample,
    tmp_dir,
):
    """
    Run bedtools intersect between the VCF and each genomic-region BED file.

    Returns:
        dict: {region: set_of_(chrom, pos)}
    """
    hits = {
        region: set()
        for region in bed_paths
    }

    for region, bed_path in bed_paths.items():
        out_path = os.path.join(
            tmp_dir,
            f"{sample}_{region}_hits.bed",
        )

        result = subprocess.run(
            f"bedtools intersect "
            f"-a {vcf_path} "
            f"-b {bed_path} "
            f"-u > {out_path}",
            shell=True,
            capture_output=True,
            text=True,
        )

        if result.returncode != 0:
            print(
                f"    [ERROR] bedtools intersect failed for {region}: "
                f"{result.stderr.strip()}",
                flush=True,
            )
            continue

        with open(out_path) as hit_file:
            for line in hit_file:
                if line.startswith("#"):
                    continue

                cols = line.strip().split("\t")

                if len(cols) >= 2:
                    hits[region].add(
                        (cols[0], cols[1])
                    )

    return hits


def assign_region_by_priority(
    hits,
    region_priority,
    all_positions,
):
    """
    Assign the highest-priority genomic region to each variant position.

    Variants that do not overlap any annotated region are assigned to
    "intergenic".
    """
    assigned = {}

    for position in all_positions:
        best_region   = "intergenic"
        best_priority = -1

        for region, position_set in hits.items():
            if position in position_set:
                priority = region_priority.get(region, 0)

                if priority > best_priority:
                    best_priority = priority
                    best_region   = region

        assigned[position] = best_region

    return assigned


def count_variants_from_vcf(vcf_path):
    """
    Read a VCF file and return a list of:
      (chrom, pos, ref, alt, type)
    """
    import gzip as gz

    def _open(path):
        if path.endswith(".gz"):
            return gz.open(path, "rt")

        return open(path)

    def _classify(ref, alt):
        if len(ref) == 1 and len(alt) == 1:
            return "SNP"
        elif len(ref) == len(alt):
            return "MNP"
        elif len(ref) < len(alt):
            return "INS"
        elif len(ref) > len(alt):
            return "DEL"
        else:
            return "COMPLEX"

    variants = []

    with _open(vcf_path) as vcf_file:
        for line in vcf_file:
            if line.startswith("#"):
                continue

            cols = line.strip().split("\t")

            if len(cols) < 5:
                continue

            chrom = cols[0]
            pos   = cols[1]
            ref   = cols[3]
            alt   = cols[4]

            for alt_allele in alt.split(","):
                variants.append(
                    (
                        chrom,
                        pos,
                        ref,
                        alt_allele,
                        _classify(ref, alt_allele),
                    )
                )

    return variants


def build_merged_beds_for_all_samples(
    samples,
    gff3_dir,
    fai_dir,
    label,
    tmp_dir,
):
    """
    Generate genomic-region BED files for every sample in a community
    and combine them into one BED file per region for the PGGB VCF.

    This is necessary because the PGGB VCF contains chromosomes from all
    genomes in the same file, for example:
      B3#0#Group1
      C26#0#Group1
      CO8#0#Group1

    Returns:
        dict: {region: merged_bed_path}
    """
    # Store BED lines from every sample by genomic region
    region_lines = {}

    for sample in samples:
        gff3 = find_gff3_for_sample(
            sample,
            gff3_dir,
        )

        fai_path, _ = find_fai_for_sample(
            sample,
            fai_dir,
        )

        if not gff3:
            print(
                f"    [WARNING] No GFF3 file found for sample "
                f"'{sample}' — skipping",
                flush=True,
            )
            continue

        sample_label = f"{label}_{sample}"

        fasta_sizes = get_chrom_sizes_normalized(
            fai_path,
            sample_label,
            tmp_dir,
        )

        sample_beds = gff3_to_region_beds(
            gff3,
            sample_label,
            tmp_dir,
            fasta_sizes,
        )

        for region, bed_path in sample_beds.items():
            with open(bed_path) as bed_file:
                lines = [
                    line
                    for line in bed_file
                    if line.strip()
                ]

            region_lines.setdefault(
                region,
                [],
            ).extend(lines)

    if not region_lines:
        return {}

    # Write one combined BED per region and merge it again
    combined_beds = {}

    for region, lines in region_lines.items():
        combined_raw = os.path.join(
            tmp_dir,
            f"{label}_ALL_{region}_raw.bed",
        )

        combined_path = os.path.join(
            tmp_dir,
            f"{label}_ALL_{region}_merged.bed",
        )

        with open(combined_raw, "w") as combined_file:
            combined_file.writelines(lines)

        subprocess.run(
            f"bedtools sort -i {combined_raw} "
            f"| bedtools merge -i - > {combined_path}",
            shell=True,
            check=True,
        )

        combined_beds[region] = combined_path

    available_regions = list(combined_beds.keys())

    print(
        f"    Combined BED files ({len(samples)} samples): "
        f"{available_regions}",
        flush=True,
    )

    return combined_beds

# =========================
# MAIN FUNCTION FOR EACH VCF
# =========================

def process_one_vcf(
    vcf_path,
    gff3_path,
    fasta_path,
    label,
    tmp_dir,
):
    """
    Process one VCF file by generating genomic-region BED files,
    intersecting the variants with those regions, and returning a DataFrame.

    Parameters:
        label:
            Identifier string, for example:
            "community.9_PGGB" or "community.9_B3_ref"
    """
    sample = (
        label
        .replace("/", "_")
        .replace(" ", "_")
    )

    print(f"\n  [{label}]", flush=True)

    # 1. Chromosome sizes:
    # Use the .fai file normalized to the PanSN format
    fai_path = fasta_path

    fasta_sizes = get_chrom_sizes_normalized(
        fai_path,
        sample,
        tmp_dir,
    )

    # 2. Generate genomic-region BED files from the GFF3
    bed_paths = gff3_to_region_beds(
        gff3_path,
        sample,
        tmp_dir,
        fasta_sizes,
    )

    if not bed_paths:
        print(
            f"    [ERROR] No BED files were generated for {label}",
            flush=True,
        )
        return None

    # 3. Read variants from the VCF
    variants = count_variants_from_vcf(vcf_path)

    if not variants:
        print(
            f"    [WARNING] Empty VCF: {vcf_path}",
            flush=True,
        )
        return None

    all_positions = {
        (chrom, pos)
        for chrom, pos, _, _, _ in variants
    }

    print(
        f"    Variants in VCF: {len(variants)} "
        f"({len(all_positions)} unique positions)",
        flush=True,
    )

    # 4. Run bedtools intersections
    hits = intersect_vcf_with_regions(
        vcf_path,
        bed_paths,
        sample,
        tmp_dir,
    )

    # 5. Assign genomic regions according to priority
    assigned = assign_region_by_priority(
        hits,
        REGION_PRIORITY,
        all_positions,
    )

    # 6. Build a DataFrame containing variant type and genomic region
    rows = []

    for chrom, pos, ref, alt, variant_type in variants:
        region = assigned.get(
            (chrom, pos),
            "intergenic",
        )

        rows.append(
            {
                "label":  label,
                "chrom":  chrom,
                "pos":    pos,
                "ref":    ref,
                "alt":    alt,
                "type":   variant_type,
                "region": region,
            }
        )

    return pd.DataFrame(rows)

# =========================
# MAIN LOOP
# =========================

# List of DataFrames used to build the final comparison table
all_results = []

for community in COMMUNITIES:
    print(f"\n{'=' * 65}", flush=True)
    print(f"Community: {community}", flush=True)
    print(f"{'=' * 65}", flush=True)

    community_output = os.path.join(
        OUTPUT_DIR,
        community,
    )

    os.makedirs(
        community_output,
        exist_ok=True,
    )

    # --------------------------------------------------
    # A) PGGB VCF: variants.vcf from the community
    # --------------------------------------------------

    community_dir = os.path.join(
        PGGB_BASE,
        f"all_pacbio_pansn.fasta.{GRAPH_HASH}.{community}",
    )

    pggb_vcf_matches = glob.glob(
        os.path.join(
            community_dir,
            "variants.vcf",
        )
    )

    pggb_vcf = (
        pggb_vcf_matches[0]
        if pggb_vcf_matches
        else None
    )

    if pggb_vcf:
        # Extract every sample from the PGGB GFA file and generate combined
        # BED files. The PGGB VCF contains chromosomes from all genomes:
        # B3#0#Group1, C26#0#Group1, CO8#0#Group1, etc.
        gfa_files = glob.glob(
            os.path.join(
                community_dir,
                "*.gfa",
            )
        )

        all_samples = []

        if gfa_files:
            with open(gfa_files[0]) as gfa_file:
                seen = set()

                for line in gfa_file:
                    if line.startswith("P\t"):
                        sample = (
                            line
                            .split("\t")[1]
                            .split("#")[0]
                        )

                        if sample not in seen:
                            all_samples.append(sample)
                            seen.add(sample)

        if all_samples:
            print(
                f"  Samples in GFA: {all_samples}",
                flush=True,
            )

            label = f"{community}_PGGB"

            pggb_beds = build_merged_beds_for_all_samples(
                all_samples,
                GFF3_DIR,
                FAI_DIR,
                label,
                tmp_dir,
            )

            if pggb_beds:
                variants = count_variants_from_vcf(
                    pggb_vcf
                )

                if variants:
                    all_positions = {
                        (chrom, pos)
                        for chrom, pos, _, _, _ in variants
                    }

                    print(
                        f"    Variants in VCF: {len(variants)} "
                        f"({len(all_positions)} unique positions)",
                        flush=True,
                    )

                    hits = intersect_vcf_with_regions(
                        pggb_vcf,
                        pggb_beds,
                        label,
                        tmp_dir,
                    )

                    assigned = assign_region_by_priority(
                        hits,
                        REGION_PRIORITY,
                        all_positions,
                    )

                    rows = []

                    for (
                        chrom,
                        pos,
                        ref,
                        alt,
                        variant_type,
                    ) in variants:
                        region = assigned.get(
                            (chrom, pos),
                            "intergenic",
                        )

                        rows.append(
                            {
                                "label":  label,
                                "chrom":  chrom,
                                "pos":    pos,
                                "ref":    ref,
                                "alt":    alt,
                                "type":   variant_type,
                                "region": region,
                            }
                        )

                    df = pd.DataFrame(rows)

                    annotated_output = os.path.join(
                        community_output,
                        f"{label}_annotated.tsv",
                    )

                    df.to_csv(
                        annotated_output,
                        sep="\t",
                        index=False,
                    )

                    all_results.append(df)

                    summary = (
                        df
                        .groupby(["type", "region"])
                        .size()
                        .unstack(fill_value=0)
                    )

                    print(
                        f"  {label}:\n{summary.to_string()}\n",
                        flush=True,
                    )

            else:
                print(
                    f"  [SKIP] No BED files were generated for {label}",
                    flush=True,
                )

        else:
            print(
                f"  [SKIP] Could not extract samples from the GFA "
                f"in {community_dir}",
                flush=True,
            )

    else:
        print(
            f"  [SKIP] variants.vcf was not found in {community_dir}",
            flush=True,
        )

    # --------------------------------------------------
    # B) Per-reference VCF files: *_as_ref.vcf
    # --------------------------------------------------

    per_ref_dir = os.path.join(
        PER_REF_BASE,
        community,
    )

    per_ref_vcfs = sorted(
        glob.glob(
            os.path.join(
                per_ref_dir,
                "*_as_ref.vcf",
            )
        )
    )

    if not per_ref_vcfs:
        print(
            f"  [SKIP] No *_as_ref.vcf files were found in "
            f"{per_ref_dir}",
            flush=True,
        )

    else:
        for vcf_path in per_ref_vcfs:
            ref_sample = (
                os.path.basename(vcf_path)
                .replace("_as_ref.vcf", "")
            )

            gff3 = find_gff3_for_sample(
                ref_sample,
                GFF3_DIR,
            )

            fai_path, _ = find_fai_for_sample(
                ref_sample,
                FAI_DIR,
            )

            if not gff3:
                print(
                    f"  [SKIP] No GFF3 file found for '{ref_sample}'",
                    flush=True,
                )
                continue

            label = f"{community}_{ref_sample}_ref"

            df = process_one_vcf(
                vcf_path,
                gff3,
                fai_path,
                label,
                tmp_dir,
            )

            if df is not None:
                annotated_output = os.path.join(
                    community_output,
                    f"{label}_annotated.tsv",
                )

                df.to_csv(
                    annotated_output,
                    sep="\t",
                    index=False,
                )

                all_results.append(df)

# =========================
# FINAL SUMMARY TABLES
# =========================

print(f"\n{'=' * 65}", flush=True)
print("Generating summary tables...", flush=True)

if not all_results:
    print(
        "[ERROR] No VCF files were processed.",
        flush=True,
    )
    sys.exit(1)

full_df = pd.concat(
    all_results,
    ignore_index=True,
)

# Table 1:
# Variant distribution by genomic region using absolute counts
region_order = [
    "CDS",
    "UTR5",
    "UTR3",
    "exon",
    "intron",
    "intergenic",
]

for label_id, group in full_df.groupby("label"):
    pivot = (
        group
        .groupby(["type", "region"])
        .size()
        .unstack(fill_value=0)
        .reindex(
            columns=[
                region
                for region in region_order
                if region in group["region"].unique()
            ],
            fill_value=0,
        )
    )

    pivot.loc["TOTAL"] = pivot.sum()

    output_table = os.path.join(
        OUTPUT_DIR,
        f"{label_id}_region_distribution.tsv",
    )

    pivot.to_csv(
        output_table,
        sep="\t",
    )

    print(
        f"\n  {label_id}:",
        flush=True,
    )

    print(
        pivot.to_string(),
        flush=True,
    )

# Table 2:
# Comparison between references using the percentage of variants
# assigned to each genomic region
summary_rows = []

for label_id, group in full_df.groupby("label"):
    total = len(group)

    row = {
        "label": label_id,
    }

    for region in region_order:
        count = (
            group["region"] == region
        ).sum()

        row[region] = (
            round(count / total * 100, 2)
            if total > 0
            else 0
        )

    summary_rows.append(row)

summary_df = (
    pd.DataFrame(summary_rows)
    .set_index("label")
)

summary_output = os.path.join(
    OUTPUT_DIR,
    "all_labels_region_summary_pct.tsv",
)

summary_df.to_csv(
    summary_output,
    sep="\t",
)

print(f"\n{'=' * 65}", flush=True)
print("Percentage summary across labels:", flush=True)
print(summary_df.to_string(), flush=True)
print(f"\nSaved to: {summary_output}", flush=True)

end_time = datetime.now()

print(
    f"\nTotal runtime: {end_time - start_time}",
    flush=True,
)

print(
    f"End: {end_time.strftime('%Y-%m-%d %H:%M:%S')}",
    flush=True,
)
