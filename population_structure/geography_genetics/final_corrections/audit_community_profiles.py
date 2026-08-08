#!/usr/bin/env python3
"""
audit_community_profiles.py

Audits the community labels used by the all11 and conservative9 profiles.

For every expected community, the script reports:
  - whether it appears in the profile manifest;
  - whether it appears in the VCF QC table;
  - raw VCF records;
  - polymorphic biallelic SNPs before final filtering;
  - retained SNPs in the final feature metadata;
  - final status and explanation.

This resolves apparent discrepancies between profile names (11 and 9)
and the number of communities visible in final SNP-contribution tables.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument("--manifest", required=True)
    parser.add_argument("--all11-vcf-qc", required=True)
    parser.add_argument("--all11-features", required=True)
    parser.add_argument("--conservative9-vcf-qc", required=True)
    parser.add_argument("--conservative9-features", required=True)
    parser.add_argument("--output-dir", required=True)

    return parser.parse_args()


def read_tsv(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str)


def numeric_column(
    dataframe: pd.DataFrame,
    column_name: str,
) -> pd.Series:
    if column_name not in dataframe.columns:
        return pd.Series(
            [0] * len(dataframe),
            index=dataframe.index,
            dtype="int64",
        )

    return (
        pd.to_numeric(
            dataframe[column_name],
            errors="coerce",
        )
        .fillna(0)
        .astype(int)
    )


def audit_profile(
    manifest: pd.DataFrame,
    manifest_column: str,
    profile_name: str,
    vcf_qc_path: str,
    feature_path: str,
) -> pd.DataFrame:
    expected = manifest.loc[
        manifest[manifest_column].str.upper() == "YES",
        ["community", "note"],
    ].copy()

    vcf_qc = read_tsv(vcf_qc_path)
    features = read_tsv(feature_path)

    if "community" not in vcf_qc.columns:
        raise ValueError(
            f"{profile_name} VCF QC lacks a community column."
        )

    if "community" not in features.columns:
        raise ValueError(
            f"{profile_name} feature metadata lacks a community column."
        )

    vcf_qc["vcf_records_numeric"] = numeric_column(
        vcf_qc,
        "vcf_records",
    )

    snp_column_candidates = [
        "polymorphic_biallelic_snps",
        "polymorphic_features",
        "snp_allele_features",
    ]

    detected_snp_column = next(
        (
            column
            for column in snp_column_candidates
            if column in vcf_qc.columns
        ),
        None,
    )

    if detected_snp_column is None:
        vcf_qc["pre_filter_snps_numeric"] = 0
    else:
        vcf_qc["pre_filter_snps_numeric"] = numeric_column(
            vcf_qc,
            detected_snp_column,
        )

    retained_counts = (
        features.groupby("community")
        .size()
        .rename("retained_final_snps")
        .reset_index()
    )

    qc_summary = (
        vcf_qc.groupby("community", as_index=False)
        .agg(
            vcf_files=("community", "size"),
            raw_vcf_records=(
                "vcf_records_numeric",
                "sum",
            ),
            polymorphic_biallelic_snps_before_final_filter=(
                "pre_filter_snps_numeric",
                "sum",
            ),
        )
    )

    audit = expected.merge(
        qc_summary,
        on="community",
        how="left",
    )

    audit = audit.merge(
        retained_counts,
        on="community",
        how="left",
    )

    fill_columns = [
        "vcf_files",
        "raw_vcf_records",
        "polymorphic_biallelic_snps_before_final_filter",
        "retained_final_snps",
    ]

    for column in fill_columns:
        audit[column] = (
            pd.to_numeric(
                audit[column],
                errors="coerce",
            )
            .fillna(0)
            .astype(int)
        )

    audit.insert(0, "profile", profile_name)
    audit["expected_in_profile"] = True
    audit["present_in_vcf_qc"] = audit["vcf_files"] > 0
    audit["contributed_final_snps"] = (
        audit["retained_final_snps"] > 0
    )

    def classify(row: pd.Series) -> str:
        if not row["present_in_vcf_qc"]:
            return "expected but absent from VCF QC"

        if (
            row[
                "polymorphic_biallelic_snps_before_final_filter"
            ]
            == 0
        ):
            return "processed but no polymorphic biallelic SNPs"

        if row["retained_final_snps"] == 0:
            return (
                "SNPs detected before filtering but none retained "
                "in final matrix"
            )

        return "contributed retained SNPs"

    audit["final_status"] = audit.apply(
        classify,
        axis=1,
    )

    total_retained = audit["retained_final_snps"].sum()

    audit["retained_snp_percentage"] = (
        audit["retained_final_snps"]
        .div(total_retained if total_retained else 1)
        .mul(100)
        .round(4)
    )

    return audit


def main() -> None:
    args = parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    manifest = read_tsv(args.manifest)

    required_manifest_columns = {
        "community",
        "analysis_all",
        "analysis_conservative",
        "note",
    }

    missing = required_manifest_columns - set(manifest.columns)

    if missing:
        raise ValueError(
            "Manifest missing columns: "
            f"{sorted(missing)}"
        )

    all11_audit = audit_profile(
        manifest=manifest,
        manifest_column="analysis_all",
        profile_name="all11",
        vcf_qc_path=args.all11_vcf_qc,
        feature_path=args.all11_features,
    )

    conservative9_audit = audit_profile(
        manifest=manifest,
        manifest_column="analysis_conservative",
        profile_name="conservative9",
        vcf_qc_path=args.conservative9_vcf_qc,
        feature_path=args.conservative9_features,
    )

    audit = pd.concat(
        [all11_audit, conservative9_audit],
        ignore_index=True,
    )

    audit.to_csv(
        output_dir / "community_profile_audit.tsv",
        sep="\t",
        index=False,
    )

    summary_rows = []

    for profile_name, group in audit.groupby("profile"):
        summary_rows.append(
            {
                "profile": profile_name,
                "communities_expected": len(group),
                "communities_in_vcf_qc": int(
                    group["present_in_vcf_qc"].sum()
                ),
                "communities_contributing_final_snps": int(
                    group["contributed_final_snps"].sum()
                ),
                "communities_without_final_snps": int(
                    (~group["contributed_final_snps"]).sum()
                ),
                "raw_vcf_records": int(
                    group["raw_vcf_records"].sum()
                ),
                "pre_filter_polymorphic_snps": int(
                    group[
                        "polymorphic_biallelic_snps_before_final_filter"
                    ].sum()
                ),
                "retained_final_snps": int(
                    group["retained_final_snps"].sum()
                ),
            }
        )

    summary = pd.DataFrame(summary_rows)

    summary.to_csv(
        output_dir / "community_profile_audit_summary.tsv",
        sep="\t",
        index=False,
    )

    print(summary.to_string(index=False))


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(
            f"[ERROR] {error}",
            file=sys.stderr,
            flush=True,
        )
        raise
