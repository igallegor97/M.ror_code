#!/usr/bin/env python3
"""
build_master_region_comparison.py

Combines all variant-by-region distribution tables into a single master table.

For each community and reference, the script:
  1. Reads the corresponding *_region_distribution.tsv file.
  2. Removes the TOTAL row.
  3. Transposes the table so that variant types become rows.
  4. Calculates absolute counts and percentages for each genomic region.
  5. Saves the combined results as TSV and Excel files.

Expected input filename format:
  community.<number>_<reference>_region_distribution.tsv

Example:
  community.9_B3_ref_region_distribution.tsv

Usage:
    python build_master_region_comparison.py

Edit the CONFIGURATION section to change input and output paths.
"""

import glob
import os
import re
import sys
from datetime import datetime

import pandas as pd

# =========================
# CONFIGURATION
# =========================

INPUT_DIR = (
    "/home/isabella_gallego/OneDrive/Documentos/Maestria/"
    "PGGB/var_analysis_comms/region_annotation"
)

OUTPUT_TSV = "master_region_comparison.tsv"
OUTPUT_XLSX = "master_region_comparison.xlsx"

FILE_PATTERN = "*_region_distribution.tsv"

# =========================
# SETUP
# =========================

start_time = datetime.now()

print("=" * 65, flush=True)
print(f"Start      : {start_time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
print(f"INPUT_DIR  : {INPUT_DIR}", flush=True)
print(f"OUTPUT_TSV : {OUTPUT_TSV}", flush=True)
print(f"OUTPUT_XLSX: {OUTPUT_XLSX}", flush=True)
print("=" * 65, flush=True)

if not os.path.isdir(INPUT_DIR):
    print(
        f"[ERROR] INPUT_DIR does not exist or is not a directory: {INPUT_DIR}",
        flush=True,
    )
    sys.exit(1)

# =========================
# FIND INPUT FILES
# =========================

files = sorted(
    glob.glob(
        os.path.join(INPUT_DIR, FILE_PATTERN)
    )
)

if not files:
    print(
        f"[ERROR] No files matching '{FILE_PATTERN}' were found in: "
        f"{INPUT_DIR}",
        flush=True,
    )
    sys.exit(1)

print(f"\nFiles found: {len(files)}\n", flush=True)

# =========================
# FUNCTIONS
# =========================

def extract_metadata(file_path):
    """
    Extract the community and reference labels from an input filename.

    Expected format:
      community.<number>_<reference>_region_distribution.tsv

    Example:
      community.9_B3_ref_region_distribution.tsv

    Returns:
        tuple: (community, reference_label)

        Returns (None, None) when the filename does not match the
        expected format.
    """
    basename = os.path.basename(file_path)

    match = re.fullmatch(
        r"(community\.\d+)_(.+?)_region_distribution\.tsv",
        basename,
    )

    if not match:
        return None, None

    community = match.group(1)
    reference_label = match.group(2)

    return community, reference_label


def process_region_distribution(file_path):
    """
    Read and transform one variant-by-region distribution table.

    The original table is expected to contain:
      - variant types as rows
      - genomic regions as columns
      - an optional TOTAL row

    Returns the transposed DataFrame, with genomic regions as rows and
    variant types as columns.
    """
    df = pd.read_csv(
        file_path,
        sep="\t",
        index_col=0,
    )

    # Remove the summary row before transposing
    df = df.drop(
        index="TOTAL",
        errors="ignore",
    )

    # Transpose so that genomic regions become rows and variant types
    # become columns
    df = df.T

    # Remove a TOTAL row if one is present after transposition
    df = df.drop(
        index="TOTAL",
        errors="ignore",
    )

    # Convert all table values to numeric values
    df = df.apply(
        pd.to_numeric,
        errors="coerce",
    ).fillna(0)

    return df


def build_rows(
    df,
    community,
    reference_label,
):
    """
    Convert one transformed distribution table into master-table rows.

    For each genomic region, absolute counts and percentages are calculated
    for every variant type.
    """
    rows = []

    for variant_type in df.index:
        row = {
            "community": community,
            "reference": reference_label,
            "variant_type": variant_type,
        }

        values = df.loc[variant_type]
        total = values.sum()

        for region in df.columns:
            value = values[region]

            row[region] = value

            if total > 0:
                row[f"pct_{region}"] = round(
                    value / total * 100,
                    2,
                )
            else:
                row[f"pct_{region}"] = 0

        row["total_variants"] = total
        rows.append(row)

    return rows

# =========================
# PROCESS INPUT FILES
# =========================

all_rows = []
skipped_files = []
failed_files = []

for file_path in files:
    community, reference_label = extract_metadata(file_path)

    if community is None:
        print(
            f"  [SKIP] Unexpected filename format: "
            f"{os.path.basename(file_path)}",
            flush=True,
        )
        skipped_files.append(file_path)
        continue

    print(
        f"  [OK] {community} | {reference_label} "
        f"→ {os.path.basename(file_path)}",
        flush=True,
    )

    try:
        region_df = process_region_distribution(file_path)

        if region_df.empty:
            print(
                "       [WARNING] The transformed table is empty.",
                flush=True,
            )
            failed_files.append(file_path)
            continue

        rows = build_rows(
            region_df,
            community,
            reference_label,
        )

        all_rows.extend(rows)

        print(
            f"       Rows added: {len(rows)}",
            flush=True,
        )

    except Exception as error:
        print(
            f"       [ERROR] Could not process file: {error}",
            flush=True,
        )
        failed_files.append(file_path)

# =========================
# BUILD MASTER TABLE
# =========================

if not all_rows:
    print(
        "\n[ERROR] No valid rows were generated.",
        flush=True,
    )
    sys.exit(1)

master_df = pd.DataFrame(all_rows)

master_df = master_df.sort_values(
    [
        "community",
        "reference",
        "variant_type",
    ]
).reset_index(drop=True)

# =========================
# SAVE RESULTS
# =========================

master_df.to_csv(
    OUTPUT_TSV,
    sep="\t",
    index=False,
)

try:
    master_df.to_excel(
        OUTPUT_XLSX,
        index=False,
    )
except ImportError:
    print(
        "[ERROR] An Excel writer is not installed.",
        flush=True,
    )
    print(
        "        Install openpyxl with: pip install openpyxl --user",
        flush=True,
    )
    sys.exit(1)

# =========================
# FINAL REPORT
# =========================

end_time = datetime.now()
elapsed = end_time - start_time

print("\n" + "=" * 65, flush=True)
print("Master table preview:", flush=True)
print(master_df.head().to_string(index=False), flush=True)

print("\n" + "=" * 65, flush=True)
print(f"Files found     : {len(files)}", flush=True)
print(f"Files processed : {len(files) - len(skipped_files) - len(failed_files)}", flush=True)
print(f"Files skipped   : {len(skipped_files)}", flush=True)
print(f"Files failed    : {len(failed_files)}", flush=True)
print(f"Rows generated  : {len(master_df)}", flush=True)
print(f"TSV output      : {os.path.abspath(OUTPUT_TSV)}", flush=True)
print(f"Excel output    : {os.path.abspath(OUTPUT_XLSX)}", flush=True)
print(f"Total runtime   : {elapsed}", flush=True)
print(f"End             : {end_time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
print("=" * 65, flush=True)
