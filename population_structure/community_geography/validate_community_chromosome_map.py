#!/usr/bin/env python3
"""
validate_community_chromosome_map.py

Validates the community-to-chromosome mapping before chromosome-level
interpretation and plotting.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


REQUIRED_COLUMNS = {
    "community",
    "chromosome",
    "syntenic_group",
    "mapping_type",
    "chromosome_length_bp",
    "start_bp",
    "end_bp",
    "note",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--map", required=True)
    parser.add_argument("--ranking", required=True)
    parser.add_argument("--output-dir", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    mapping = pd.read_csv(args.map, sep="\t", dtype=str).fillna("")
    ranking = pd.read_csv(args.ranking, sep="\t", dtype=str).fillna("")

    missing_columns = REQUIRED_COLUMNS - set(mapping.columns)
    if missing_columns:
        raise ValueError(
            f"Mapping table missing columns: {sorted(missing_columns)}"
        )

    if mapping["community"].duplicated().any():
        duplicates = mapping.loc[
            mapping["community"].duplicated(), "community"
        ].tolist()
        raise ValueError(f"Duplicated communities in mapping: {duplicates}")

    analyzed = set(ranking["community"])
    mapped = set(mapping["community"])

    issues = []

    for community in sorted(analyzed - mapped):
        issues.append(
            {
                "community": community,
                "issue": "analyzed community absent from mapping table",
            }
        )

    analyzed_mapping = mapping.loc[
        mapping["community"].isin(analyzed)
    ].copy()

    for _, row in analyzed_mapping.iterrows():
        if not row["chromosome"].strip():
            issues.append(
                {
                    "community": row["community"],
                    "issue": "chromosome field is empty",
                }
            )

        if row["mapping_type"] in {
            "partial_or_ambiguous",
            "mixed_assignment",
        }:
            issues.append(
                {
                    "community": row["community"],
                    "issue": (
                        "mapping is ambiguous and must not be interpreted "
                        "as a unique chromosome"
                    ),
                }
            )

        numeric_fields = [
            row["chromosome_length_bp"],
            row["start_bp"],
            row["end_bp"],
        ]

        populated = [value.strip() != "" for value in numeric_fields]

        if any(populated) and not all(populated):
            issues.append(
                {
                    "community": row["community"],
                    "issue": (
                        "chromosome_length_bp/start_bp/end_bp must either "
                        "all be filled or all be empty"
                    ),
                }
            )

        if all(populated):
            try:
                chromosome_length = int(float(row["chromosome_length_bp"]))
                start = int(float(row["start_bp"]))
                end = int(float(row["end_bp"]))
            except ValueError:
                issues.append(
                    {
                        "community": row["community"],
                        "issue": "non-numeric chromosome coordinates",
                    }
                )
                continue

            if not (0 <= start < end <= chromosome_length):
                issues.append(
                    {
                        "community": row["community"],
                        "issue": (
                            "invalid coordinates; require "
                            "0 <= start < end <= chromosome_length"
                        ),
                    }
                )

    issue_df = pd.DataFrame(
        issues,
        columns=["community", "issue"],
    )

    issue_df.to_csv(
        output_dir / "community_chromosome_map_issues.tsv",
        sep="\t",
        index=False,
    )

    validated = mapping.loc[
        mapping["community"].isin(analyzed)
    ].copy()

    validated.to_csv(
        output_dir / "community_chromosome_map_validated.tsv",
        sep="\t",
        index=False,
    )

    summary = pd.DataFrame(
        [
            ["analyzed_communities", len(analyzed)],
            [
                "analyzed_communities_with_mapping_rows",
                len(analyzed & mapped),
            ],
            [
                "analyzed_communities_with_nonempty_chromosome",
                int(
                    validated["chromosome"]
                    .astype(str)
                    .str.strip()
                    .ne("")
                    .sum()
                ),
            ],
            ["mapping_issues", len(issue_df)],
        ],
        columns=["metric", "value"],
    )

    summary.to_csv(
        output_dir / "community_chromosome_map_validation_summary.tsv",
        sep="\t",
        index=False,
    )

    print(summary.to_string(index=False))

    if (
        validated["chromosome"]
        .astype(str)
        .str.strip()
        .eq("")
        .any()
    ):
        print(
            "\n[WARNING] Some analyzed communities lack a validated "
            "chromosome assignment. Chromosome figures will label them "
            "as unmapped."
        )


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"[ERROR] {error}", file=sys.stderr, flush=True)
        raise
