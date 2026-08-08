#!/usr/bin/env python3
"""
rank_communities.py

Integrates per-community Mantel, MRM, regression, SNP contribution, and
chromosome mapping into a transparent community ranking.

Primary ranking criteria:
  1. MRM geographic coefficient (effect size)
  2. MRM model R2
  3. Mantel Pearson r
  4. Statistical support after Benjamini-Hochberg FDR correction

An exploratory composite score is also reported. It is the mean percentile
rank of positive geographic-effect metrics and should not replace the
individual effect sizes in biological interpretation.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser()

    parser.add_argument("--mantel", required=True)
    parser.add_argument("--mrm-coefficients", required=True)
    parser.add_argument("--mrm-models", required=True)
    parser.add_argument("--regression", required=True)
    parser.add_argument("--community-summary", required=True)
    parser.add_argument("--community-map", required=True)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--subset", default="high_precision")
    parser.add_argument("--output-dir", required=True)

    return parser.parse_args()


def bh_fdr(values: pd.Series) -> pd.Series:
    """
    Apply the Benjamini-Hochberg false-discovery-rate correction.

    Missing values are preserved as NaN.
    """

    numeric = pd.to_numeric(
        values,
        errors="coerce",
    )

    result = pd.Series(
        np.nan,
        index=values.index,
        dtype=float,
    )

    valid = numeric.notna()

    if not valid.any():
        return result

    p_values = numeric.loc[valid].to_numpy()

    order = np.argsort(p_values)
    ranked_p_values = p_values[order]
    number_tests = len(ranked_p_values)

    adjusted = (
        ranked_p_values
        * number_tests
        / np.arange(1, number_tests + 1)
    )

    adjusted = np.minimum.accumulate(
        adjusted[::-1]
    )[::-1]

    adjusted = np.minimum(
        adjusted,
        1.0,
    )

    restored = np.empty(number_tests)
    restored[order] = adjusted

    result.loc[valid] = restored

    return result


def percentile_rank(series: pd.Series) -> pd.Series:
    """Return percentile ranks after coercing values to numeric."""

    numeric = pd.to_numeric(
        series,
        errors="coerce",
    )

    return numeric.rank(
        method="average",
        pct=True,
    )


def require_columns(
    dataframe: pd.DataFrame,
    required_columns: set[str],
    dataframe_name: str,
) -> None:
    """Raise an informative error if required columns are missing."""

    missing_columns = (
        required_columns
        - set(dataframe.columns)
    )

    if missing_columns:
        raise ValueError(
            f"{dataframe_name} is missing required columns: "
            f"{sorted(missing_columns)}. "
            f"Available columns: {sorted(dataframe.columns)}"
        )


def main() -> None:
    """Run the community-ranking workflow."""

    args = parse_args()

    output_dir = Path(args.output_dir)

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    # ========================================================================
    # Load input tables
    # ========================================================================

    mantel = pd.read_csv(
        args.mantel,
        sep="\t",
    )

    mrm_coefficients = pd.read_csv(
        args.mrm_coefficients,
        sep="\t",
    )

    mrm_models = pd.read_csv(
        args.mrm_models,
        sep="\t",
    )

    regression = pd.read_csv(
        args.regression,
        sep="\t",
    )

    community_summary = pd.read_csv(
        args.community_summary,
        sep="\t",
    )

    community_map = pd.read_csv(
        args.community_map,
        sep="\t",
        dtype=str,
    )

    # ========================================================================
    # Validate input tables
    # ========================================================================

    require_columns(
        mantel,
        {
            "profile",
            "community",
            "subset",
            "method",
            "statistic",
            "p_value",
        },
        "Mantel results",
    )

    require_columns(
        mrm_coefficients,
        {
            "profile",
            "community",
            "subset",
            "term",
            "estimate",
            "p_value",
        },
        "MRM coefficient results",
    )

    require_columns(
        mrm_models,
        {
            "profile",
            "community",
            "subset",
            "r_squared",
            "r_squared_p_value",
            "f_statistic",
            "f_test_p_value",
        },
        "MRM model results",
    )

    require_columns(
        regression,
        {
            "profile",
            "community",
            "subset",
            "slope_log10_km",
            "r_squared",
        },
        "Regression results",
    )

    require_columns(
        community_summary,
        {
            "profile",
            "community",
            "retained_snps",
            "variable_snps",
            "samples",
            "analysis_status",
        },
        "Community analysis summary",
    )

    require_columns(
        community_map,
        {
            "community",
            "chromosome",
            "syntenic_group",
            "mapping_type",
        },
        "Community-to-chromosome map",
    )

    # ========================================================================
    # Select primary Mantel result
    # ========================================================================

    mantel_primary = mantel.loc[
        (mantel["subset"] == args.subset)
        & (mantel["method"] == "pearson")
    ].copy()

    # retained_snps and variable_snps already exist in community_summary.
    # They must NOT be imported again from the Mantel table.
    #
    # Otherwise pandas generates:
    #   retained_snps_x
    #   retained_snps_y
    #
    # and the canonical retained_snps column disappears.
    mantel_primary = mantel_primary[
        [
            "profile",
            "community",
            "statistic",
            "p_value",
        ]
    ].rename(
        columns={
            "statistic": "mantel_pearson_r",
            "p_value": "mantel_pearson_p",
        }
    )

    # ========================================================================
    # Select MRM geographic coefficient
    # ========================================================================

    geo_terms = mrm_coefficients.loc[
        (mrm_coefficients["subset"] == args.subset)
        & mrm_coefficients["term"].str.contains(
            "geographic|log10",
            case=False,
            regex=True,
            na=False,
        )
    ].copy()

    geo_terms = geo_terms[
        [
            "profile",
            "community",
            "estimate",
            "p_value",
        ]
    ].rename(
        columns={
            "estimate": "mrm_geographic_beta",
            "p_value": "mrm_geographic_p",
        }
    )

    # ========================================================================
    # Select MRM technology coefficient
    # ========================================================================

    technology_terms = mrm_coefficients.loc[
        (mrm_coefficients["subset"] == args.subset)
        & mrm_coefficients["term"].str.contains(
            "technology",
            case=False,
            regex=True,
            na=False,
        )
    ].copy()

    technology_terms = technology_terms[
        [
            "profile",
            "community",
            "estimate",
            "p_value",
        ]
    ].rename(
        columns={
            "estimate": "mrm_technology_beta",
            "p_value": "mrm_technology_p",
        }
    )

    # ========================================================================
    # Select MRM model-level statistics
    # ========================================================================

    model_primary = mrm_models.loc[
        mrm_models["subset"] == args.subset
    ].copy()

    model_primary = model_primary[
        [
            "profile",
            "community",
            "r_squared",
            "r_squared_p_value",
            "f_statistic",
            "f_test_p_value",
        ]
    ].rename(
        columns={
            "r_squared": "mrm_r_squared",
            "r_squared_p_value": "mrm_r_squared_p",
            "f_statistic": "mrm_f_statistic",
            "f_test_p_value": "mrm_f_test_p",
        }
    )

    # ========================================================================
    # Select descriptive regression results
    # ========================================================================

    regression_primary = regression.loc[
        regression["subset"] == args.subset
    ].copy()

    regression_primary = regression_primary[
        [
            "profile",
            "community",
            "slope_log10_km",
            "r_squared",
        ]
    ].rename(
        columns={
            "slope_log10_km": "descriptive_slope_log10_km",
            "r_squared": "descriptive_r_squared",
        }
    )

    # ========================================================================
    # Start ranking from canonical community summary
    # ========================================================================

    ranking = community_summary.loc[
        community_summary["analysis_status"]
        == "completed"
    ].copy()

    if ranking.empty:
        raise ValueError(
            "No communities with analysis_status='completed' "
            "were found in the community analysis summary."
        )

    ranking["retained_snps"] = pd.to_numeric(
        ranking["retained_snps"],
        errors="raise",
    )

    ranking["variable_snps"] = pd.to_numeric(
        ranking["variable_snps"],
        errors="raise",
    )

    ranking["samples"] = pd.to_numeric(
        ranking["samples"],
        errors="raise",
    )

    # ========================================================================
    # Merge statistical tables
    # ========================================================================

    statistical_tables = [
        (
            "Mantel results",
            mantel_primary,
        ),
        (
            "MRM geographic coefficients",
            geo_terms,
        ),
        (
            "MRM technology coefficients",
            technology_terms,
        ),
        (
            "MRM model results",
            model_primary,
        ),
        (
            "Descriptive regression results",
            regression_primary,
        ),
    ]

    for table_name, dataframe in statistical_tables:

        duplicated = dataframe.duplicated(
            subset=[
                "profile",
                "community",
            ]
        )

        if duplicated.any():

            duplicated_rows = dataframe.loc[
                dataframe.duplicated(
                    subset=[
                        "profile",
                        "community",
                    ],
                    keep=False,
                ),
                [
                    "profile",
                    "community",
                ],
            ]

            raise ValueError(
                f"{table_name} contains duplicated "
                f"profile/community rows: "
                f"{duplicated_rows.to_dict(orient='records')}"
            )

        ranking = ranking.merge(
            dataframe,
            on=[
                "profile",
                "community",
            ],
            how="left",
            validate="one_to_one",
        )

    # ========================================================================
    # Merge chromosome mapping
    # ========================================================================

    if community_map[
        "community"
    ].duplicated().any():

        duplicated_communities = community_map.loc[
            community_map[
                "community"
            ].duplicated(
                keep=False,
            ),
            "community",
        ].tolist()

        raise ValueError(
            "Community-to-chromosome map contains "
            "duplicated communities: "
            f"{duplicated_communities}"
        )

    ranking = ranking.merge(
        community_map,
        on="community",
        how="left",
        validate="many_to_one",
    )

    # ========================================================================
    # Validate merged ranking table
    # ========================================================================

    required_ranking_columns = {
        "retained_snps",
        "variable_snps",
        "mantel_pearson_r",
        "mantel_pearson_p",
        "mrm_geographic_beta",
        "mrm_geographic_p",
        "mrm_r_squared",
        "mrm_r_squared_p",
    }

    missing_ranking_columns = (
        required_ranking_columns
        - set(ranking.columns)
    )

    if missing_ranking_columns:

        raise ValueError(
            "Ranking table is missing required columns "
            "after merging: "
            f"{sorted(missing_ranking_columns)}. "
            f"Available columns: "
            f"{sorted(ranking.columns)}"
        )

    # Detect accidental duplicated columns generated by merges.
    suffix_columns = [
        column
        for column in ranking.columns
        if column.endswith("_x")
        or column.endswith("_y")
    ]

    if suffix_columns:

        raise ValueError(
            "Unexpected merge suffix columns were created: "
            f"{suffix_columns}. "
            "This indicates duplicated non-key columns "
            "between input tables."
        )

    numeric_result_columns = [
        "mantel_pearson_r",
        "mantel_pearson_p",
        "mrm_geographic_beta",
        "mrm_geographic_p",
        "mrm_technology_beta",
        "mrm_technology_p",
        "mrm_r_squared",
        "mrm_r_squared_p",
        "mrm_f_statistic",
        "mrm_f_test_p",
        "descriptive_slope_log10_km",
        "descriptive_r_squared",
    ]

    for column in numeric_result_columns:

        if column in ranking.columns:

            ranking[column] = pd.to_numeric(
                ranking[column],
                errors="coerce",
            )

    # ========================================================================
    # Multiple-testing correction
    # ========================================================================

    ranking["mantel_pearson_q"] = bh_fdr(
        ranking["mantel_pearson_p"]
    )

    ranking["mrm_geographic_q"] = bh_fdr(
        ranking["mrm_geographic_p"]
    )

    ranking["mrm_model_q"] = bh_fdr(
        ranking["mrm_r_squared_p"]
    )

    ranking["mrm_technology_q"] = bh_fdr(
        ranking["mrm_technology_p"]
    )

    # ========================================================================
    # Exploratory effect-size score
    # ========================================================================

    ranking[
        "mantel_effect_percentile"
    ] = percentile_rank(
        ranking["mantel_pearson_r"]
    )

    ranking[
        "mrm_beta_percentile"
    ] = percentile_rank(
        ranking["mrm_geographic_beta"]
    )

    ranking[
        "mrm_r2_percentile"
    ] = percentile_rank(
        ranking["mrm_r_squared"]
    )

    ranking[
        "exploratory_composite_score"
    ] = ranking[
        [
            "mantel_effect_percentile",
            "mrm_beta_percentile",
            "mrm_r2_percentile",
        ]
    ].mean(
        axis=1,
        skipna=True,
    )

    # ========================================================================
    # Evidence classification
    # ========================================================================

    ranking[
        "supported_by_mantel_fdr"
    ] = (
        ranking["mantel_pearson_q"]
        <= 0.05
    )

    ranking[
        "supported_by_mrm_coefficient_fdr"
    ] = (
        ranking["mrm_geographic_q"]
        <= 0.05
    )

    ranking[
        "supported_by_mrm_model_fdr"
    ] = (
        ranking["mrm_model_q"]
        <= 0.05
    )

    ranking[
        "positive_geographic_effect"
    ] = (
        (
            ranking["mantel_pearson_r"]
            > 0
        )
        & (
            ranking["mrm_geographic_beta"]
            > 0
        )
    )

    support_columns = [
        "supported_by_mantel_fdr",
        "supported_by_mrm_coefficient_fdr",
        "supported_by_mrm_model_fdr",
    ]

    ranking["evidence_tier"] = np.select(
        [
            (
                ranking[
                    support_columns
                ].all(axis=1)
            )
            & ranking[
                "positive_geographic_effect"
            ],
            (
                ranking[
                    support_columns
                ].sum(axis=1)
                >= 2
            )
            & ranking[
                "positive_geographic_effect"
            ],
            ranking[
                "positive_geographic_effect"
            ],
        ],
        [
            (
                "Tier 1: concordant FDR-supported "
                "geographic signal"
            ),
            (
                "Tier 2: partially FDR-supported "
                "geographic signal"
            ),
            (
                "Tier 3: positive exploratory signal"
            ),
        ],
        default=(
            "No positive geographic signal"
        ),
    )

    # ========================================================================
    # Sort and rank communities
    # ========================================================================

    ranking = ranking.sort_values(
        [
            "exploratory_composite_score",
            "mrm_geographic_beta",
            "mantel_pearson_r",
        ],
        ascending=[
            False,
            False,
            False,
        ],
        na_position="last",
    ).reset_index(
        drop=True,
    )

    ranking.insert(
        0,
        "rank",
        np.arange(
            1,
            len(ranking) + 1,
        ),
    )

    total_snps = ranking[
        "retained_snps"
    ].sum()

    ranking[
        "retained_snp_pct_among_analyzed"
    ] = (
        ranking["retained_snps"]
        / total_snps
        * 100
        if total_snps > 0
        else 0
    )

    # ========================================================================
    # Save complete ranking
    # ========================================================================

    ranking_output = (
        output_dir
        / (
            f"{args.profile}_"
            "community_geographic_ranking.tsv"
        )
    )

    ranking.to_csv(
        ranking_output,
        sep="\t",
        index=False,
    )

    # ========================================================================
    # Select top communities
    # ========================================================================

    top = ranking.loc[
        ranking[
            "evidence_tier"
        ].str.startswith(
            (
                "Tier 1",
                "Tier 2",
            ),
            na=False,
        )
    ].copy()

    if top.empty:

        top = ranking.head(
            min(
                3,
                len(ranking),
            )
        ).copy()

        top[
            "selection_reason"
        ] = (
            "Top exploratory communities; "
            "no Tier 1/2 community found"
        )

    else:

        top[
            "selection_reason"
        ] = top[
            "evidence_tier"
        ]

    top_output = (
        output_dir
        / (
            f"{args.profile}_"
            "top_geographic_communities.tsv"
        )
    )

    top.to_csv(
        top_output,
        sep="\t",
        index=False,
    )

    # ========================================================================
    # Print summary
    # ========================================================================

    display_columns = [
        "rank",
        "community",
        "chromosome",
        "retained_snps",
        "mantel_pearson_r",
        "mantel_pearson_q",
        "mrm_geographic_beta",
        "mrm_geographic_q",
        "mrm_r_squared",
        "mrm_model_q",
        "evidence_tier",
    ]

    available_display_columns = [
        column
        for column in display_columns
        if column in ranking.columns
    ]

    print(
        ranking[
            available_display_columns
        ].to_string(
            index=False,
        )
    )

    print(
        (
            "\n[OK] Complete ranking saved to: "
            f"{ranking_output}"
        ),
        flush=True,
    )

    print(
        (
            "[OK] Top-community table saved to: "
            f"{top_output}"
        ),
        flush=True,
    )


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
