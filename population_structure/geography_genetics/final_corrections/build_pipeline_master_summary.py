#!/usr/bin/env python3
"""
build_pipeline_master_summary.py

Creates a single master summary for the corrected PGGB24
geography-versus-genetics analysis.

Outputs:
  pipeline_master_summary.tsv
  key_results_summary.tsv
  extreme_pair_summary.tsv
  pipeline_master_summary.md
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument("--coordinate-summary", required=True)
    parser.add_argument("--pair-table", required=True)
    parser.add_argument("--mantel", required=True)
    parser.add_argument("--mrm-coefficients", required=True)
    parser.add_argument("--mrm-models", required=True)
    parser.add_argument("--regression", required=True)
    parser.add_argument("--all11-robustness", required=True)
    parser.add_argument("--conservative9-robustness", required=True)
    parser.add_argument("--all11-matrix-summary", required=True)
    parser.add_argument(
        "--conservative9-matrix-summary",
        required=True,
    )
    parser.add_argument("--community-audit", required=True)
    parser.add_argument("--colocation-summary", required=True)
    parser.add_argument("--output-dir", required=True)

    return parser.parse_args()


def read_tsv(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def scalar(value: Any) -> Any:
    if pd.isna(value):
        return ""

    if isinstance(value, float):
        return round(value, 8)

    return value


def append_dataframe(
    rows: list[dict[str, Any]],
    section: str,
    dataframe: pd.DataFrame,
    identifier_columns: list[str],
    value_columns: list[str],
    note: str = "",
) -> None:
    for _, record in dataframe.iterrows():
        identifiers = "; ".join(
            f"{column}={record[column]}"
            for column in identifier_columns
            if column in record.index
        )

        for value_column in value_columns:
            if value_column not in record.index:
                continue

            rows.append(
                {
                    "section": section,
                    "record": identifiers,
                    "metric": value_column,
                    "value": scalar(record[value_column]),
                    "note": note,
                }
            )


def metric_value(
    summary: pd.DataFrame,
    metric_name: str,
) -> Any:
    row = summary.loc[summary["metric"] == metric_name]

    if row.empty:
        return ""

    return row.iloc[0]["value"]


def extract_robust_metric(
    dataframe: pd.DataFrame,
    metric: str,
    statistic: str,
) -> Any:
    row = dataframe.loc[
        (dataframe["metric"] == metric)
        & (dataframe["statistic"] == statistic)
    ]

    if row.empty:
        return ""

    return row.iloc[0]["value"]


def main() -> None:
    args = parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    coordinate_summary = read_tsv(args.coordinate_summary)
    pairs = read_tsv(args.pair_table)
    mantel = read_tsv(args.mantel)
    mrm_coefficients = read_tsv(args.mrm_coefficients)
    mrm_models = read_tsv(args.mrm_models)
    regression = read_tsv(args.regression)
    all11_robustness = read_tsv(args.all11_robustness)
    conservative9_robustness = read_tsv(
        args.conservative9_robustness
    )
    all11_matrix = read_tsv(args.all11_matrix_summary)
    conservative9_matrix = read_tsv(
        args.conservative9_matrix_summary
    )
    community_audit = read_tsv(args.community_audit)
    colocation_summary = read_tsv(
        args.colocation_summary
    )

    rows: list[dict[str, Any]] = []

    append_dataframe(
        rows,
        "coordinate_validation",
        coordinate_summary,
        identifier_columns=[],
        value_columns=["value"],
    )

    append_dataframe(
        rows,
        "mantel",
        mantel,
        identifier_columns=[
            "profile",
            "subset",
            "method",
            "n_samples",
        ],
        value_columns=[
            "statistic",
            "p_value",
            "permutations",
        ],
        note=(
            "Permutation-based matrix correlation; "
            "the first distance matrix is permuted."
        ),
    )

    append_dataframe(
        rows,
        "mrm_coefficients",
        mrm_coefficients,
        identifier_columns=[
            "profile",
            "subset",
            "n_samples",
            "term",
        ],
        value_columns=[
            "estimate",
            "p_value",
            "permutations",
        ],
        note=(
            "MRM coefficient and permutation p-value; "
            "model includes geographic and technology distances."
        ),
    )

    append_dataframe(
        rows,
        "mrm_model_summary",
        mrm_models,
        identifier_columns=[
            "profile",
            "subset",
            "n_samples",
        ],
        value_columns=[
            "r_squared",
            "r_squared_p_value",
            "f_statistic",
            "f_test_p_value",
            "permutations",
        ],
    )

    append_dataframe(
        rows,
        "descriptive_regression",
        regression,
        identifier_columns=[
            "profile",
            "subset",
            "n_pairs",
        ],
        value_columns=[
            "intercept",
            "slope_log10_km",
            "r_squared",
            "adjusted_r_squared",
        ],
        note=(
            "Effect-size summary only; ordinary LM p-values "
            "are not used because pair rows are non-independent."
        ),
    )

    append_dataframe(
        rows,
        "colocation_sensitivity",
        colocation_summary,
        identifier_columns=[
            "profile",
            "method",
            "combinations",
        ],
        value_columns=[
            "statistic_median",
            "statistic_q025",
            "statistic_q975",
            "p_value_median",
            "significant_fraction_0.05",
        ],
    )

    append_dataframe(
        rows,
        "community_audit",
        community_audit,
        identifier_columns=[
            "profile",
            "community",
            "final_status",
        ],
        value_columns=[
            "raw_vcf_records",
            (
                "polymorphic_biallelic_snps_"
                "before_final_filter"
            ),
            "retained_final_snps",
            "retained_snp_percentage",
        ],
    )

    master = pd.DataFrame(rows)

    master.to_csv(
        output_dir / "pipeline_master_summary.tsv",
        sep="\t",
        index=False,
    )

    key_rows = []

    for profile_name, matrix_summary in [
        ("all11", all11_matrix),
        ("conservative9", conservative9_matrix),
    ]:
        robustness = (
            all11_robustness
            if profile_name == "all11"
            else conservative9_robustness
        )

        primary_mantel = mantel.loc[
            (mantel["profile"] == profile_name)
            & (
                mantel["subset"]
                == "all_valid_coordinates"
            )
            & (mantel["method"] == "pearson")
        ]

        high_precision_mantel = mantel.loc[
            (mantel["profile"] == profile_name)
            & (mantel["subset"] == "high_precision")
            & (mantel["method"] == "pearson")
        ]

        geographic_mrm = mrm_coefficients.loc[
            (mrm_coefficients["profile"] == profile_name)
            & (
                mrm_coefficients["subset"]
                == "all_valid_coordinates"
            )
            & mrm_coefficients["term"].str.contains(
                "geographic|log10",
                case=False,
                regex=True,
            )
        ]

        mrm_model = mrm_models.loc[
            (mrm_models["profile"] == profile_name)
            & (
                mrm_models["subset"]
                == "all_valid_coordinates"
            )
        ]

        audit_profile = community_audit.loc[
            community_audit["profile"] == profile_name
        ]

        key_rows.append(
            {
                "profile": profile_name,
                "samples_in_genetic_matrix": metric_value(
                    matrix_summary,
                    "samples",
                ),
                "retained_snps": metric_value(
                    matrix_summary,
                    "features_retained_final",
                )
                or metric_value(
                    matrix_summary,
                    "raw_polymorphic_features",
                ),
                "spatial_samples": int(
                    len(
                        set(pairs["sample_1"])
                        | set(pairs["sample_2"])
                    )
                ),
                "pairwise_comparisons": len(pairs),
                "zero_km_pairs": int(
                    pairs["zero_km_pair"].sum()
                ),
                "expected_communities": len(
                    audit_profile
                ),
                "communities_contributing_final_snps": int(
                    audit_profile[
                        "contributed_final_snps"
                    ].sum()
                ),
                "mantel_pearson_r": (
                    primary_mantel.iloc[0]["statistic"]
                    if not primary_mantel.empty
                    else ""
                ),
                "mantel_pearson_p": (
                    primary_mantel.iloc[0]["p_value"]
                    if not primary_mantel.empty
                    else ""
                ),
                "high_precision_mantel_r": (
                    high_precision_mantel.iloc[0][
                        "statistic"
                    ]
                    if not high_precision_mantel.empty
                    else ""
                ),
                "high_precision_mantel_p": (
                    high_precision_mantel.iloc[0][
                        "p_value"
                    ]
                    if not high_precision_mantel.empty
                    else ""
                ),
                "mrm_geographic_estimate": (
                    geographic_mrm.iloc[0]["estimate"]
                    if not geographic_mrm.empty
                    else ""
                ),
                "mrm_geographic_p": (
                    geographic_mrm.iloc[0]["p_value"]
                    if not geographic_mrm.empty
                    else ""
                ),
                "mrm_r_squared": (
                    mrm_model.iloc[0]["r_squared"]
                    if not mrm_model.empty
                    else ""
                ),
                "mrm_r_squared_p": (
                    mrm_model.iloc[0][
                        "r_squared_p_value"
                    ]
                    if not mrm_model.empty
                    else ""
                ),
                "balanced_pearson_median": (
                    extract_robust_metric(
                        robustness,
                        "pearson",
                        "median",
                    )
                ),
                "balanced_pearson_q025": (
                    extract_robust_metric(
                        robustness,
                        "pearson",
                        "q025.2.5%",
                    )
                ),
                "balanced_pearson_q975": (
                    extract_robust_metric(
                        robustness,
                        "pearson",
                        "q975.97.5%",
                    )
                ),
                "balanced_slope_median": (
                    extract_robust_metric(
                        robustness,
                        "slope_log10_km",
                        "median",
                    )
                ),
                "balanced_r_squared_median": (
                    extract_robust_metric(
                        robustness,
                        "r_squared",
                        "median",
                    )
                ),
            }
        )

    key_summary = pd.DataFrame(key_rows)

    key_summary.to_csv(
        output_dir / "key_results_summary.tsv",
        sep="\t",
        index=False,
    )

    extreme_rows = []

    for profile_name, residual_column in [
        ("all11", "standardized_residual_all11"),
        (
            "conservative9",
            "standardized_residual_conservative9",
        ),
    ]:
        ordered = pairs.sort_values(residual_column)

        selected = pd.concat(
            [
                ordered.head(5).assign(
                    extreme_type=(
                        "genetically closer than fitted model"
                    )
                ),
                ordered.tail(5).assign(
                    extreme_type=(
                        "genetically farther than fitted model"
                    )
                ),
            ],
            ignore_index=True,
        )

        selected.insert(0, "profile", profile_name)

        extreme_rows.append(
            selected[
                [
                    "profile",
                    "extreme_type",
                    "sample_1",
                    "sample_2",
                    "geographic_distance_km",
                    (
                        "genetic_distance_all11"
                        if profile_name == "all11"
                        else (
                            "genetic_distance_"
                            "conservative9"
                        )
                    ),
                    residual_column,
                    "zero_km_pair",
                    "zero_km_interpretation",
                ]
            ]
        )

    extreme_pairs = pd.concat(
        extreme_rows,
        ignore_index=True,
    )

    extreme_pairs.to_csv(
        output_dir / "extreme_pair_summary.tsv",
        sep="\t",
        index=False,
    )

    markdown_lines = [
        "# PGGB24 geography-versus-genetics master summary",
        "",
        "## Core profile metrics",
        "",
        "```text",
        key_summary.to_csv(sep="\t", index=False).strip(),
        "```",
        "",
        "## Interpretation rules",
        "",
        (
            "- Mantel and MRM p-values are based on "
            "permutation tests."
        ),
        (
            "- Ordinary regression is retained only as an "
            "effect-size description."
        ),
        (
            "- Zero kilometres means co-location at the "
            "available municipal/locality resolution, "
            "not necessarily the same farm."
        ),
        (
            "- Residual-ranked pairs are relative to the "
            "fitted distance-decay model."
        ),
        "",
        "## Community audit",
        "",
        "```text",
        community_audit.to_csv(sep="\t", index=False).strip(),
        "```",
    ]

    (
        output_dir / "pipeline_master_summary.md"
    ).write_text(
        "\n".join(markdown_lines),
        encoding="utf-8",
    )

    print(key_summary.to_string(index=False))


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
