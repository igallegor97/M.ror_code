#!/usr/bin/env python3
"""
build_snp_matrix.py

Builds a unified binary SNP matrix from partitioned PGGB or minigraph-cactus
VCF files.

Verified naming conventions:
  PGGB:   C26, CO8, CO84, E7
  Cactus: Mror_C26_GroupN, Mror_CO8_GroupN, Mror_CO84_GroupN, Mror_E7_GroupN

Cactus sample names are normalized to C26, CO8, CO84, and E7. B3 is absent
from the FORMAT columns because it defines the linear reference. The script
verifies that every VCF uses B3 contigs and then reconstructs B3 as 0 for each
retained SNP ALT allele.

Encoding:
  0  = ALT allele absent
  1  = ALT allele present
  NA = missing genotype

Multiallelic SNP records are decomposed into one binary feature per ALT allele.
"""

import argparse
import gzip
import os
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True)
    parser.add_argument("--vcf-list", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--sample-name-mode", choices=["pggb", "cactus"], required=True)
    parser.add_argument("--reference-sample", default="B3")
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--max-missing", type=float, default=0.20)
    return parser.parse_args()


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r", encoding="utf-8", errors="replace")


def normalize_sample(name, mode):
    if mode == "pggb":
        return name
    match = re.fullmatch(r"Mror_(.+)_Group\d+", name)
    if not match:
        raise ValueError(f"Unexpected Cactus sample name: {name}")
    return match.group(1)


def read_metadata(path):
    metadata = pd.read_csv(path, sep="\t", dtype=str)
    required = {"sample_id", "region"}
    missing = required - set(metadata.columns)
    if missing:
        raise ValueError(f"Metadata is missing columns: {sorted(missing)}")
    metadata = metadata.dropna(subset=["sample_id"]).copy()
    if metadata["sample_id"].duplicated().any():
        raise ValueError("Metadata contains duplicate sample_id values")
    return metadata


def read_vcf_list(path):
    vcfs = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            value = line.strip()
            if not value or value.startswith("#"):
                continue
            if not os.path.isfile(value):
                raise FileNotFoundError(value)
            vcfs.append(value)
    if not vcfs:
        raise ValueError(f"No VCF paths found in {path}")
    return vcfs


def parse_header(vcf_path, mode):
    contigs = []
    raw_samples = []
    with open_text(vcf_path) as handle:
        for line in handle:
            if line.startswith("##contig=<ID="):
                contigs.append(line.split("##contig=<ID=", 1)[1].split(",", 1)[0].split(">", 1)[0])
            elif line.startswith("#CHROM"):
                raw_samples = line.rstrip("\n").split("\t")[9:]
                break
    normalized = [normalize_sample(sample, mode) for sample in raw_samples]
    if len(normalized) != len(set(normalized)):
        raise ValueError(f"Duplicate normalized sample names in {vcf_path}: {normalized}")
    return raw_samples, normalized, contigs


def validate_reference(vcf_path, contigs, mode, reference):
    expected = f"{reference}#" if mode == "pggb" else f"Mror_{reference}_Group"
    bad = [contig for contig in contigs if not contig.startswith(expected)]
    if bad:
        raise ValueError(f"VCF does not consistently use {reference} as reference: {vcf_path}; unexpected contigs: {bad[:10]}")


def parse_gt(format_field, sample_field):
    keys = format_field.split(":")
    if "GT" not in keys:
        return None
    index = keys.index("GT")
    values = sample_field.split(":")
    if index >= len(values):
        return None
    gt = values[index]
    if gt in {".", "./.", ".|."}:
        return None
    alleles = []
    for token in gt.replace("|", "/").split("/"):
        if token == ".":
            return None
        try:
            alleles.append(int(token))
        except ValueError:
            return None
    return alleles


def build_matrix(vcfs, metadata, mode, reference, max_missing):
    samples = metadata["sample_id"].tolist()
    if reference not in samples:
        raise ValueError(f"Reference sample {reference} is absent from metadata")
    non_reference = [sample for sample in samples if sample != reference]

    rows, feature_ids, feature_info, file_qc = [], [], [], []
    total_records = 0

    for index, vcf_path in enumerate(vcfs, start=1):
        print(f"  [{index}/{len(vcfs)}] {vcf_path}", flush=True)
        raw_samples, normalized_samples, contigs = parse_header(vcf_path, mode)
        validate_reference(vcf_path, contigs, mode, reference)

        missing = [sample for sample in non_reference if sample not in normalized_samples]
        unexpected = [sample for sample in normalized_samples if sample not in non_reference]
        if missing or unexpected:
            raise ValueError(f"Sample mismatch in {vcf_path}; missing={missing}; unexpected={unexpected}")

        column_by_sample = {
            normalized: 9 + raw_samples.index(raw)
            for raw, normalized in zip(raw_samples, normalized_samples)
        }

        prefix = Path(vcf_path).name.removesuffix(".gz").removesuffix(".vcf")
        file_records = 0
        file_features = 0

        with open_text(vcf_path) as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 10:
                    continue
                total_records += 1
                file_records += 1
                chrom, pos, record_id, ref = fields[0], fields[1], fields[2], fields[3]
                if len(ref) != 1:
                    continue
                format_field = fields[8]

                for alt_index, alt in enumerate(fields[4].split(","), start=1):
                    if len(alt) != 1 or alt in {".", "*"}:
                        continue

                    values_by_sample = {reference: 0.0}
                    for sample in non_reference:
                        alleles = parse_gt(format_field, fields[column_by_sample[sample]])
                        values_by_sample[sample] = np.nan if alleles is None else float(alt_index in alleles)

                    feature_id = f"{prefix}|{chrom}:{pos}:{ref}>{alt}|ALT{alt_index}"
                    feature_ids.append(feature_id)
                    rows.append([values_by_sample[sample] for sample in samples])
                    feature_info.append({
                        "feature_id": feature_id,
                        "vcf_file": vcf_path,
                        "chrom": chrom,
                        "pos": pos,
                        "record_id": record_id,
                        "ref": ref,
                        "alt": alt,
                        "alt_index": alt_index,
                        "reference_sample_reconstructed": reference,
                    })
                    file_features += 1

        file_qc.append({
            "vcf_file": vcf_path,
            "raw_vcf_samples": ",".join(raw_samples),
            "normalized_vcf_samples": ",".join(normalized_samples),
            "reference_sample_added": reference,
            "vcf_records": file_records,
            "snp_allele_features": file_features,
        })

    if not rows:
        raise ValueError("No SNP features were extracted")

    raw = pd.DataFrame(rows, index=feature_ids, columns=samples, dtype=float)
    missing_fraction = raw.isna().mean(axis=1)
    polymorphic = raw.min(axis=1, skipna=True).ne(raw.max(axis=1, skipna=True))
    passes_missing = missing_fraction.le(max_missing)
    retained = polymorphic & passes_missing
    filtered = raw.loc[retained].copy()
    if filtered.empty:
        raise ValueError("No SNP features remained after filtering")
    imputed = filtered.apply(lambda row: row.fillna(row.mean()), axis=1)

    features = pd.DataFrame(feature_info).set_index("feature_id")
    features["missing_fraction"] = missing_fraction
    features["polymorphic"] = polymorphic
    features["passes_missing_filter"] = passes_missing
    features["retained"] = retained

    summary = {
        "vcf_files": len(vcfs),
        "vcf_records_total": total_records,
        "snp_allele_features_total": len(raw),
        "polymorphic_features": int(polymorphic.sum()),
        "features_passing_missing_filter": int(passes_missing.sum()),
        "features_retained": int(retained.sum()),
        "samples_total": raw.shape[1],
        "samples_from_vcf": len(non_reference),
        "reference_sample": reference,
        "reference_samples_reconstructed": 1,
        "max_missing": max_missing,
    }
    return raw, filtered, imputed, features, pd.DataFrame(file_qc), summary


def main():
    args = parse_args()
    start = datetime.now()
    if not 0 <= args.max_missing <= 1:
        raise ValueError("--max-missing must be between 0 and 1")

    output_prefix = Path(args.output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    print("=" * 65)
    print("SNP MATRIX CONSTRUCTION")
    print(f"Source             : {args.source}")
    print(f"Sample-name mode   : {args.sample_name_mode}")
    print(f"Reference sample   : {args.reference_sample}")
    print(f"Maximum missingness: {args.max_missing:.2f}")
    print("=" * 65)

    try:
        metadata = read_metadata(args.metadata)
        vcfs = read_vcf_list(args.vcf_list)
        raw, filtered, imputed, features, qc, summary = build_matrix(
            vcfs, metadata, args.sample_name_mode, args.reference_sample, args.max_missing
        )
    except Exception as error:
        print(f"[ERROR] {error}", flush=True)
        sys.exit(1)

    raw.to_csv(f"{output_prefix}_snp_matrix_raw.tsv", sep="\t", na_rep="NA")
    filtered.to_csv(f"{output_prefix}_snp_matrix_filtered.tsv", sep="\t", na_rep="NA")
    imputed.to_csv(f"{output_prefix}_snp_matrix_imputed.tsv", sep="\t")
    features.to_csv(f"{output_prefix}_snp_feature_metadata.tsv", sep="\t")
    qc.to_csv(f"{output_prefix}_vcf_qc.tsv", sep="\t", index=False)

    with open(f"{output_prefix}_matrix_summary.tsv", "w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"source\t{args.source}\n")
        for key, value in summary.items():
            handle.write(f"{key}\t{value}\n")

    end = datetime.now()
    print("\n" + "=" * 65)
    print("SNP MATRIX SUMMARY")
    for key, value in summary.items():
        print(f"{key:<36}: {value}")
    print(f"Total runtime                       : {end - start}")
    print("=" * 65)


if __name__ == "__main__":
    main()
