#!/usr/bin/env python3
"""
sample_balanced_snps.py

Creates a partition-balanced SNP matrix from the complete imputed SNP matrix.

The script samples an exact total number of SNP features while distributing
them as evenly as possible across genomic partitions.

Supported partitions:
  PGGB: community.0 ... community.9
  Cactus: group1 ... group11

Example with 1,000 total SNPs:
  PGGB: 100 SNPs per community
  Cactus: 90 or 91 SNPs per chromosome group

Usage:
    python sample_balanced_snps.py \
        --source PGGB \
        --matrix PGGB_snp_matrix_imputed.tsv \
        --features PGGB_snp_feature_metadata.tsv \
        --total-snps 1000 \
        --seed 42 \
        --output-prefix PGGB_balanced
"""

import argparse
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Sample an exact number of SNPs evenly across partitions."
    )
    parser.add_argument("--source", required=True, choices=["PGGB", "Cactus"])
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--features", required=True)
    parser.add_argument("--total-snps", required=True, type=int)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output-prefix", required=True)
    return parser.parse_args()


def extract_partition(vcf_file, source):
    path = str(vcf_file)

    if source == "PGGB":
        match = re.search(r"community\.(\d+)", path)
        if match:
            return f"community.{int(match.group(1))}"
    else:
        match = re.search(r"(?:^|/)group(\d+)(?:/|$)", path)
        if not match:
            match = re.search(r"cactus_group(\d+)_", path)
        if match:
            return f"group{int(match.group(1))}"

    raise ValueError(f"Could not extract a {source} partition from: {path}")


def natural_key(label):
    match = re.search(r"(\d+)$", label)
    return int(match.group(1)) if match else label


def allocate_exact_total(partition_sizes, total_snps, seed):
    partitions = sorted(partition_sizes, key=natural_key)
    available_total = sum(partition_sizes.values())

    if total_snps > available_total:
        raise ValueError(
            f"Requested {total_snps} SNPs, but only {available_total} are available."
        )

    allocation = {partition: 0 for partition in partitions}
    remaining = total_snps
    rng = np.random.default_rng(seed)

    while remaining > 0:
        eligible = [
            partition
            for partition in partitions
            if allocation[partition] < partition_sizes[partition]
        ]

        if not eligible:
            raise RuntimeError("Unable to complete SNP allocation.")

        minimum = min(allocation[partition] for partition in eligible)
        lowest = [
            partition
            for partition in eligible
            if allocation[partition] == minimum
        ]

        chosen = str(rng.choice(lowest))
        allocation[chosen] += 1
        remaining -= 1

    return allocation


def main():
    args = parse_arguments()
    start_time = datetime.now()

    if args.total_snps < 1:
        print("[ERROR] --total-snps must be positive.", flush=True)
        sys.exit(1)

    matrix_path = Path(args.matrix)
    feature_path = Path(args.features)
    output_prefix = Path(args.output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    if not matrix_path.is_file():
        print(f"[ERROR] Matrix not found: {matrix_path}", flush=True)
        sys.exit(1)

    if not feature_path.is_file():
        print(f"[ERROR] Feature metadata not found: {feature_path}", flush=True)
        sys.exit(1)

    print("=" * 65, flush=True)
    print("PARTITION-BALANCED SNP SAMPLING", flush=True)
    print(f"Start         : {start_time:%Y-%m-%d %H:%M:%S}", flush=True)
    print(f"Source        : {args.source}", flush=True)
    print(f"Total SNPs    : {args.total_snps}", flush=True)
    print(f"Random seed   : {args.seed}", flush=True)
    print(f"Output prefix : {output_prefix}", flush=True)
    print("=" * 65, flush=True)

    matrix = pd.read_csv(matrix_path, sep="\t", index_col=0)
    features = pd.read_csv(feature_path, sep="\t", index_col=0)

    if "vcf_file" not in features.columns:
        print("[ERROR] Feature metadata lacks the 'vcf_file' column.", flush=True)
        sys.exit(1)

    missing = matrix.index.difference(features.index)
    if len(missing) > 0:
        print(
            f"[ERROR] {len(missing)} matrix features are missing from metadata.",
            flush=True,
        )
        sys.exit(1)

    features = features.loc[matrix.index].copy()

    try:
        features["partition"] = features["vcf_file"].map(
            lambda value: extract_partition(value, args.source)
        )
    except ValueError as error:
        print(f"[ERROR] {error}", flush=True)
        sys.exit(1)

    partition_sizes = features["partition"].value_counts().to_dict()

    try:
        allocation = allocate_exact_total(
            partition_sizes,
            args.total_snps,
            args.seed,
        )
    except (ValueError, RuntimeError) as error:
        print(f"[ERROR] {error}", flush=True)
        sys.exit(1)

    rng = np.random.default_rng(args.seed)
    selected_features = []

    for partition in sorted(allocation, key=natural_key):
        candidates = features.index[
            features["partition"] == partition
        ].to_numpy()

        selected = rng.choice(
            candidates,
            size=allocation[partition],
            replace=False,
        )

        selected_features.extend(selected.tolist())

    selected_features = rng.permutation(selected_features).tolist()

    balanced_matrix = matrix.loc[selected_features].copy()
    selected_metadata = features.loc[selected_features].copy()

    if balanced_matrix.shape[0] != args.total_snps:
        print("[ERROR] Balanced matrix has an unexpected number of rows.", flush=True)
        sys.exit(1)

    sampling_rows = []
    for partition in sorted(partition_sizes, key=natural_key):
        available = partition_sizes[partition]
        selected = allocation[partition]
        sampling_rows.append(
            {
                "partition": partition,
                "available_snps": available,
                "selected_snps": selected,
                "selected_pct": round(selected / available * 100, 4),
            }
        )

    sampling_summary = pd.DataFrame(sampling_rows)

    matrix_output = f"{output_prefix}_snp_matrix_imputed.tsv"
    metadata_output = f"{output_prefix}_selected_feature_metadata.tsv"
    list_output = f"{output_prefix}_selected_features.txt"
    partition_output = f"{output_prefix}_partition_sampling.tsv"
    summary_output = f"{output_prefix}_sampling_summary.tsv"

    balanced_matrix.to_csv(matrix_output, sep="\t")
    selected_metadata.to_csv(metadata_output, sep="\t")
    sampling_summary.to_csv(partition_output, sep="\t", index=False)

    with open(list_output, "w", encoding="utf-8") as handle:
        for feature_id in selected_features:
            handle.write(f"{feature_id}\n")

    with open(summary_output, "w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"source\t{args.source}\n")
        handle.write(f"random_seed\t{args.seed}\n")
        handle.write(f"partitions\t{len(partition_sizes)}\n")
        handle.write(f"available_snps\t{matrix.shape[0]}\n")
        handle.write(f"selected_snps\t{balanced_matrix.shape[0]}\n")
        handle.write(f"samples\t{balanced_matrix.shape[1]}\n")

    end_time = datetime.now()

    print("\nPer-partition sampling:", flush=True)
    print(sampling_summary.to_string(index=False), flush=True)
    print("\n" + "=" * 65, flush=True)
    print("BALANCED SNP SAMPLING SUMMARY", flush=True)
    print(f"Partitions     : {len(partition_sizes)}", flush=True)
    print(f"Available SNPs : {matrix.shape[0]}", flush=True)
    print(f"Selected SNPs  : {balanced_matrix.shape[0]}", flush=True)
    print(f"Samples        : {balanced_matrix.shape[1]}", flush=True)
    print(f"Total runtime  : {end_time - start_time}", flush=True)
    print(f"End            : {end_time:%Y-%m-%d %H:%M:%S}", flush=True)
    print("=" * 65, flush=True)


if __name__ == "__main__":
    main()
