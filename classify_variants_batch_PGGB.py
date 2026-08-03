#!/usr/bin/env python3
"""
classify_variants_batch.py

Variant classification by community (PGGB output).
Reads configuration from environment variables defined in the SLURM job.

Required environment variables:
  BASE_DIR     - Root directory containing community folders
  VCF_NAME     - Name of the VCF file inside each folder (default: variants.vcf.gz)
  OUTPUT_DIR   - Output directory
  COMM_PATTERN - Glob pattern for community folders (default: *.community.*)
"""

import pandas as pd
import gzip
import os
import glob
import sys
from datetime import datetime

# =========================
# CONFIGURATION FROM ENV
# =========================

BASE_DIR     = os.environ.get("BASE_DIR")
VCF_NAME     = os.environ.get("VCF_NAME",     "variants.vcf.gz")
OUTPUT_DIR   = os.environ.get("OUTPUT_DIR",   "variant_classification_results")
COMM_PATTERN = os.environ.get("COMM_PATTERN", "*.community.*")

if not BASE_DIR:
    print("[ERROR] Environment variable BASE_DIR is not defined.", flush=True)
    sys.exit(1)

if not os.path.isdir(BASE_DIR):
    print(f"[ERROR] BASE_DIR does not exist or is not a directory: {BASE_DIR}", flush=True)
    sys.exit(1)

os.makedirs(OUTPUT_DIR, exist_ok=True)

start_time = datetime.now()
print("=" * 60, flush=True)
print(f"Start        : {start_time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
print(f"BASE_DIR     : {BASE_DIR}",     flush=True)
print(f"VCF_NAME     : {VCF_NAME}",     flush=True)
print(f"OUTPUT_DIR   : {OUTPUT_DIR}",   flush=True)
print(f"COMM_PATTERN : {COMM_PATTERN}", flush=True)
print("=" * 60, flush=True)

# =========================
# FIND COMMUNITIES
# =========================

community_dirs = sorted(glob.glob(os.path.join(BASE_DIR, COMM_PATTERN)))

if not community_dirs:
    print(f"[ERROR] No community folders matching '{COMM_PATTERN}' were found in: {BASE_DIR}", flush=True)
    sys.exit(1)

print(f"\nCommunities found: {len(community_dirs)}\n", flush=True)

# =========================
# FUNCTIONS
# =========================

def classify_variant(ref, alt_allele):
    """Classify a variant by comparing the lengths of REF and ALT."""
    if len(ref) == 1 and len(alt_allele) == 1:
        return "SNP"
    elif len(ref) == len(alt_allele):
        return "MNP"
    elif len(ref) < len(alt_allele):
        return "INS"
    elif len(ref) > len(alt_allele):
        return "DEL"
    else:
        return "COMPLEX"


def open_vcf(vcf_path):
    """Open a compressed (.gz) or plain-text VCF file."""
    if vcf_path.endswith(".gz"):
        return gzip.open(vcf_path, "rt")
    return open(vcf_path, "r")


def process_community(comm_dir, vcf_name):
    """Process the VCF of a community. Returns (comm_id, DataFrame | None)."""
    comm_name = os.path.basename(comm_dir)
    comm_id   = comm_name.split("community.")[-1] if "community." in comm_name else comm_name
    vcf_path  = os.path.join(comm_dir, vcf_name)

    if not os.path.exists(vcf_path):
        print(f"  [SKIP]  comm_{comm_id}: '{vcf_name}' not found", flush=True)
        return comm_id, None

    print(f"  [OK]    comm_{comm_id} → {vcf_path}", flush=True)

    variants = []
    try:
        with open_vcf(vcf_path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                cols = line.strip().split("\t")
                if len(cols) < 5:
                    continue
                chrom = cols[0]
                pos   = int(cols[1])
                ref   = cols[3]
                alt   = cols[4]
                for alt_allele in alt.split(","):
                    variants.append({
                        "community": comm_id,
                        "chrom":     chrom,
                        "pos":       pos,
                        "ref":       ref,
                        "alt":       alt_allele,
                        "type":      classify_variant(ref, alt_allele)
                    })
    except Exception as e:
        print(f"  [ERROR] comm_{comm_id}: {e}", flush=True)
        return comm_id, None

    if not variants:
        print(f"  [WARNING] comm_{comm_id}: no variants found.", flush=True)
        return comm_id, None

    return comm_id, pd.DataFrame(variants)

# =========================
# COMMUNITY LOOP
# =========================

all_summaries = []
failed        = []

for comm_dir in community_dirs:
    comm_id, df = process_community(comm_dir, VCF_NAME)

    if df is None:
        failed.append(comm_id)
        continue

    # Complete variant table for the community
    out_tsv = os.path.join(OUTPUT_DIR, f"comm_{comm_id}_variants_classified.tsv")
    df.to_csv(out_tsv, sep="\t", index=False)

    # Summary by variant type
    summary = df["type"].value_counts().sort_index()
    out_sum = os.path.join(OUTPUT_DIR, f"comm_{comm_id}_summary.tsv")
    summary.to_csv(out_sum, sep="\t", header=["count"])

    print(f"           Variants: {len(df)} | {summary.to_dict()}", flush=True)

    row = {"community": comm_id}
    row.update(summary.to_dict())
    all_summaries.append(row)

# =========================
# FINAL COMPARISON TABLE
# =========================

print("\n" + "=" * 60, flush=True)
print("Comparison table across communities:", flush=True)

if all_summaries:
    comp_df = (
        pd.DataFrame(all_summaries)
        .fillna(0)
        .set_index("community")
    )
    ordered_cols = [c for c in ["SNP", "MNP", "INS", "DEL", "COMPLEX"] if c in comp_df.columns]
    comp_df = comp_df[ordered_cols].astype(int)
    comp_df.loc["TOTAL"] = comp_df.sum()

    out_comp = os.path.join(OUTPUT_DIR, "all_communities_comparison.tsv")
    comp_df.to_csv(out_comp, sep="\t")

    print(comp_df.to_string(), flush=True)
    print(f"\nSaved to: {out_comp}", flush=True)
else:
    print("[WARNING] No communities were successfully processed.", flush=True)

# =========================
# FINAL REPORT
# =========================

end_time = datetime.now()
elapsed  = end_time - start_time

print("\n" + "=" * 60, flush=True)
print(f"Communities processed : {len(all_summaries)}", flush=True)
print(f"Communities failed    : {len(failed)} {failed}", flush=True)
print(f"Total runtime         : {elapsed}", flush=True)
print(f"End : {end_time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
print("=" * 60, flush=True)

if failed and not all_summaries:
    sys.exit(1)
