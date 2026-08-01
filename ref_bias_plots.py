#!/usr/bin/env python3
"""
plot_reference_bias.py

Generates publication-ready figures to evaluate reference-dependent variation
in variant counts.

The script reads the comparison table produced by
classify_per_reference_vcfs.py and creates four figures:

  1. Fig1_grouped_barplot.pdf
     Absolute variant counts by reference genome and variant type.

  2. Fig2_heatmap_normalized.pdf
     Percentage contribution of each variant type within each reference.

  3. Fig3_stacked_barplot.pdf
     Relative variant composition for each reference genome.

  4. Fig4_reference_bias.pdf
     Two-panel summary showing:
       A. Mean variant counts and the observed minimum-to-maximum range.
       B. Relative variation across references, calculated as:
          (maximum - minimum) / mean × 100.

The script also saves the reference-bias statistics as a TSV file.

Expected input table:
  - Rows represent reference genomes.
  - Columns represent variant types.
  - An optional TOTAL row may be present.
  - Expected variant-type columns include:
      SNP, MNP, INS, DEL, and COMPLEX.

Usage:
    python plot_reference_bias.py

Configuration may be provided through environment variables or by editing
the default values in the CONFIGURATION section.
"""

import os
import sys
from datetime import datetime

import matplotlib

# Use a non-interactive backend suitable for HPC and headless environments.
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import seaborn as sns

# =========================
# CONFIGURATION
# =========================

INPUT_TSV = os.environ.get(
    "INPUT_TSV",
    (
        "/home/isabella_gallego/OneDrive/Documentos/Maestria/PGGB/"
        "comunidades/per_reference_vcfs/community.9/"
        "classification_results/community.9_reference_comparison.tsv"
    ),
)

OUTPUT_DIR = os.environ.get(
    "OUTPUT_DIR",
    (
        "/home/isabella_gallego/OneDrive/Documentos/Maestria/PGGB/"
        "comunidades/per_reference_vcfs/community.9/"
        "classification_results/figures"
    ),
)

COMMUNITY = os.environ.get(
    "COMMUNITY",
    "community.9",
)

# Preferred variant-type order
VARIANT_TYPE_ORDER = [
    "SNP",
    "MNP",
    "INS",
    "DEL",
    "COMPLEX",
]

# Accessible color palette based on Wong (2011)
VARIANT_COLORS = {
    "SNP":     "#0072B2",
    "MNP":     "#56B4E9",
    "INS":     "#009E73",
    "DEL":     "#E69F00",
    "COMPLEX": "#CC79A7",
}

# Optional visual reference line for the relative-variation panel.
# This is not a formal statistical threshold.
REFERENCE_LINE_PCT = 10.0

# =========================
# PLOT SETTINGS
# =========================

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": [
            "Arial",
            "Helvetica",
            "DejaVu Sans",
        ],
        "font.size": 11,
        "axes.titlesize": 12,
        "axes.titleweight": "bold",
        "axes.labelsize": 11,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "legend.frameon": False,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.format": "pdf",
    }
)

# =========================
# FUNCTIONS
# =========================


def load_reference_comparison(input_path):
    """
    Read and validate the reference-comparison table.

    The optional TOTAL row is removed before plotting. Variant-count columns
    are converted to numeric values.

    Parameters
    ----------
    input_path : str
        Path to the tab-separated comparison table.

    Returns
    -------
    tuple
        A tuple containing:
          - validated DataFrame
          - ordered list of variant types present in the table
          - list of reference labels
    """
    if not os.path.isfile(input_path):
        print(
            f"[ERROR] Input table not found: {input_path}",
            flush=True,
        )
        sys.exit(1)

    try:
        raw_df = pd.read_csv(
            input_path,
            sep="\t",
            index_col=0,
        )
    except Exception as error:
        print(
            f"[ERROR] Could not read the input table: {error}",
            flush=True,
        )
        sys.exit(1)

    if raw_df.empty:
        print(
            f"[ERROR] Input table is empty: {input_path}",
            flush=True,
        )
        sys.exit(1)

    # Remove the summary row before plotting
    data_df = raw_df.drop(
        index="TOTAL",
        errors="ignore",
    ).copy()

    if data_df.empty:
        print(
            "[ERROR] No reference rows remain after removing the TOTAL row.",
            flush=True,
        )
        sys.exit(1)

    variant_types = [
        variant_type
        for variant_type in VARIANT_TYPE_ORDER
        if variant_type in data_df.columns
    ]

    if not variant_types:
        print(
            "[ERROR] No recognized variant-type columns were found.",
            flush=True,
        )
        print(
            "        Expected one or more of: "
            + ", ".join(VARIANT_TYPE_ORDER),
            flush=True,
        )
        sys.exit(1)

    # Convert recognized variant columns to numeric values
    data_df[variant_types] = data_df[variant_types].apply(
        pd.to_numeric,
        errors="coerce",
    )

    invalid_mask = data_df[variant_types].isna()

    if invalid_mask.any().any():
        invalid_columns = invalid_mask.any(axis=0)
        invalid_columns = invalid_columns[
            invalid_columns
        ].index.tolist()

        print(
            "[ERROR] Non-numeric or missing values were found in columns: "
            + ", ".join(invalid_columns),
            flush=True,
        )
        sys.exit(1)

    if (data_df[variant_types] < 0).any().any():
        print(
            "[ERROR] Negative variant counts were found in the input table.",
            flush=True,
        )
        sys.exit(1)

    references = data_df.index.astype(str).tolist()
    data_df.index = references

    return data_df, variant_types, references


def calculate_normalized_composition(data_df, variant_types):
    """
    Calculate the percentage contribution of each variant type per reference.

    Rows with a total count of zero are assigned zero percentages rather than
    producing missing or infinite values.

    Parameters
    ----------
    data_df : pandas.DataFrame
        Absolute variant-count table.

    variant_types : list of str
        Variant-type columns to include.

    Returns
    -------
    pandas.DataFrame
        Row-normalized percentage table.
    """
    counts = data_df[variant_types].copy()
    row_totals = counts.sum(axis=1)

    normalized_df = counts.div(
        row_totals.replace(0, np.nan),
        axis=0,
    ) * 100

    return normalized_df.fillna(0)


def calculate_reference_bias(data_df, variant_types):
    """
    Calculate reference-dependent variation for each variant type.

    The relative-variation index is defined as:

        (maximum - minimum) / mean × 100

    This is a descriptive range-based metric. It should not be interpreted as
    a formal statistical test of reference bias.

    Parameters
    ----------
    data_df : pandas.DataFrame
        Absolute variant-count table.

    variant_types : list of str
        Variant types to evaluate.

    Returns
    -------
    pandas.DataFrame
        Summary containing minimum, maximum, mean, range, and relative
        variation for each variant type.
    """
    bias_records = []

    for variant_type in variant_types:
        values = data_df[variant_type].astype(float)

        minimum = values.min()
        maximum = values.max()
        mean = values.mean()
        count_range = maximum - minimum

        relative_variation = (
            count_range / mean * 100
            if mean > 0
            else 0.0
        )

        bias_records.append(
            {
                "variant_type": variant_type,
                "minimum": minimum,
                "maximum": maximum,
                "mean": mean,
                "range": count_range,
                "relative_variation_pct": relative_variation,
            }
        )

    return (
        pd.DataFrame(bias_records)
        .set_index("variant_type")
    )


def save_figure(figure, output_path):
    """
    Save and close a Matplotlib figure.

    Parameters
    ----------
    figure : matplotlib.figure.Figure
        Figure to save.

    output_path : str
        Destination PDF path.
    """
    figure.savefig(
        output_path,
        format="pdf",
    )

    plt.close(figure)

    print(
        f"  [OK] Saved: {output_path}",
        flush=True,
    )


def plot_grouped_barplot(
    data_df,
    variant_types,
    references,
    output_path,
):
    """
    Generate a grouped bar plot of absolute variant counts.
    """
    number_of_references = len(references)
    number_of_types = len(variant_types)

    x_positions = np.arange(
        number_of_references
    )

    bar_width = 0.8 / number_of_types

    figure, axis = plt.subplots(
        figsize=(7, 4.5)
    )

    for index, variant_type in enumerate(variant_types):
        offset = (
            index
            - number_of_types / 2
            + 0.5
        ) * bar_width

        axis.bar(
            x_positions + offset,
            data_df[variant_type].to_numpy(),
            width=bar_width * 0.9,
            color=VARIANT_COLORS[variant_type],
            label=variant_type,
            zorder=3,
        )

    axis.set_xticks(
        x_positions
    )

    axis.set_xticklabels(
        references,
        rotation=0,
    )

    axis.set_xlabel(
        "Reference genome"
    )

    axis.set_ylabel(
        "Number of variants"
    )

    axis.set_title(
        f"Variant counts per reference — {COMMUNITY}"
    )

    axis.yaxis.set_major_formatter(
        mticker.FuncFormatter(
            lambda value, _: f"{int(value):,}"
        )
    )

    axis.grid(
        axis="y",
        linestyle="--",
        linewidth=0.5,
        alpha=0.6,
        zorder=0,
    )

    axis.legend(
        title="Variant type",
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
    )

    figure.tight_layout()

    save_figure(
        figure,
        output_path,
    )


def plot_normalized_heatmap(
    normalized_df,
    output_path,
):
    """
    Generate a heatmap of variant-type percentages by reference genome.
    """
    figure_width = max(
        5,
        0.9 * normalized_df.shape[1] + 2,
    )

    figure_height = max(
        3.5,
        0.55 * normalized_df.shape[0] + 1.5,
    )

    figure, axis = plt.subplots(
        figsize=(
            figure_width,
            figure_height,
        )
    )

    sns.heatmap(
        normalized_df,
        ax=axis,
        cmap="YlOrRd",
        annot=True,
        fmt=".1f",
        linewidths=0.5,
        linecolor="white",
        vmin=0,
        vmax=100,
        cbar_kws={
            "label": "% of total variants",
            "shrink": 0.8,
        },
        annot_kws={
            "size": 9,
        },
    )

    axis.set_title(
        f"Variant composition per reference — {COMMUNITY}"
    )

    axis.set_xlabel(
        "Variant type"
    )

    axis.set_ylabel(
        "Reference genome"
    )

    axis.tick_params(
        axis="x",
        rotation=0,
    )

    axis.tick_params(
        axis="y",
        rotation=0,
    )

    figure.tight_layout()

    save_figure(
        figure,
        output_path,
    )


def plot_stacked_barplot(
    normalized_df,
    variant_types,
    references,
    output_path,
):
    """
    Generate a stacked bar plot of relative variant composition.
    """
    number_of_references = len(
        references
    )

    x_positions = np.arange(
        number_of_references
    )

    figure, axis = plt.subplots(
        figsize=(6, 4.5)
    )

    bottom_values = np.zeros(
        number_of_references
    )

    for variant_type in variant_types:
        percentages = normalized_df[
            variant_type
        ].to_numpy()

        axis.bar(
            x_positions,
            percentages,
            bottom=bottom_values,
            color=VARIANT_COLORS[variant_type],
            label=variant_type,
            width=0.55,
            zorder=3,
        )

        # Add labels only when the segment is large enough to remain legible
        for reference_index, (
            percentage,
            bottom,
        ) in enumerate(
            zip(
                percentages,
                bottom_values,
            )
        ):
            if percentage > 4:
                axis.text(
                    reference_index,
                    bottom + percentage / 2,
                    f"{percentage:.1f}%",
                    ha="center",
                    va="center",
                    fontsize=8,
                    color="white",
                    fontweight="bold",
                )

        bottom_values += percentages

    axis.set_xticks(
        x_positions
    )

    axis.set_xticklabels(
        references,
        rotation=0,
    )

    axis.set_ylim(
        0,
        100,
    )

    axis.set_xlabel(
        "Reference genome"
    )

    axis.set_ylabel(
        "Proportion of variants (%)"
    )

    axis.set_title(
        f"Variant composition per reference — {COMMUNITY}"
    )

    axis.grid(
        axis="y",
        linestyle="--",
        linewidth=0.5,
        alpha=0.6,
        zorder=0,
    )

    axis.legend(
        title="Variant type",
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
    )

    figure.tight_layout()

    save_figure(
        figure,
        output_path,
    )


def plot_reference_bias(
    bias_df,
    output_path,
):
    """
    Generate a two-panel reference-bias summary.

    Panel A shows the mean count and observed minimum-to-maximum range.

    Panel B shows the descriptive relative-variation index:
        (maximum - minimum) / mean × 100
    """
    variant_types = bias_df.index.tolist()

    x_positions = np.arange(
        len(variant_types)
    )

    bar_colors = [
        VARIANT_COLORS[variant_type]
        for variant_type in variant_types
    ]

    figure, (
        count_axis,
        variation_axis,
    ) = plt.subplots(
        1,
        2,
        figsize=(10, 4.5),
    )

    figure.suptitle(
        f"Reference-dependent variation — {COMMUNITY}",
        fontsize=13,
        fontweight="bold",
        y=1.01,
    )

    # ---------------------------------------------------------
    # Panel A: Mean counts with minimum-to-maximum ranges
    # ---------------------------------------------------------

    lower_errors = (
        bias_df["mean"]
        - bias_df["minimum"]
    ).to_numpy()

    upper_errors = (
        bias_df["maximum"]
        - bias_df["mean"]
    ).to_numpy()

    count_axis.bar(
        x_positions,
        bias_df["mean"].to_numpy(),
        color=bar_colors,
        width=0.55,
        zorder=3,
    )

    count_axis.errorbar(
        x=x_positions,
        y=bias_df["mean"].to_numpy(),
        yerr=np.vstack(
            [
                lower_errors,
                upper_errors,
            ]
        ),
        fmt="none",
        color="black",
        capsize=5,
        linewidth=1.5,
        zorder=4,
    )

    count_axis.set_xticks(
        x_positions
    )

    count_axis.set_xticklabels(
        variant_types
    )

    count_axis.set_xlabel(
        "Variant type"
    )

    count_axis.set_ylabel(
        "Number of variants (mean and observed range)"
    )

    count_axis.set_title(
        "A   Absolute counts",
        loc="left",
        fontweight="bold",
    )

    count_axis.yaxis.set_major_formatter(
        mticker.FuncFormatter(
            lambda value, _: f"{int(value):,}"
        )
    )

    count_axis.grid(
        axis="y",
        linestyle="--",
        linewidth=0.5,
        alpha=0.6,
        zorder=0,
    )

    maximum_count = bias_df[
        "maximum"
    ].max()

    annotation_offset = (
        maximum_count * 0.025
        if maximum_count > 0
        else 1
    )

    for index, (_, row) in enumerate(
        bias_df.iterrows()
    ):
        count_axis.text(
            index,
            row["maximum"] + annotation_offset,
            f"{row['mean']:,.0f}",
            ha="center",
            va="bottom",
            fontsize=8.5,
            color="black",
        )

    # Add headroom for labels above the error bars
    count_axis.set_ylim(
        bottom=0,
        top=(
            maximum_count + annotation_offset * 4
            if maximum_count > 0
            else 1
        ),
    )

    # ---------------------------------------------------------
    # Panel B: Relative variation
    # ---------------------------------------------------------

    variation_axis.bar(
        x_positions,
        bias_df["relative_variation_pct"].to_numpy(),
        color=bar_colors,
        width=0.55,
        zorder=3,
    )

    variation_axis.axhline(
        REFERENCE_LINE_PCT,
        color="gray",
        linestyle="--",
        linewidth=1,
        zorder=2,
        label=f"{REFERENCE_LINE_PCT:g}% visual reference",
    )

    variation_axis.set_xticks(
        x_positions
    )

    variation_axis.set_xticklabels(
        variant_types
    )

    maximum_variation = bias_df[
        "relative_variation_pct"
    ].max()

    variation_offset = max(
        maximum_variation * 0.025,
        0.5,
    )

    for index, (_, row) in enumerate(
        bias_df.iterrows()
    ):
        variation_axis.text(
            index,
            row["relative_variation_pct"] + variation_offset,
            f"{row['relative_variation_pct']:.1f}%",
            ha="center",
            va="bottom",
            fontsize=9,
            fontweight="bold",
            color="black",
        )

    variation_axis.set_ylim(
        bottom=0,
        top=max(
            maximum_variation + variation_offset * 4,
            REFERENCE_LINE_PCT * 1.25,
            1,
        ),
    )

    variation_axis.set_xlabel(
        "Variant type"
    )

    variation_axis.set_ylabel(
        "Relative variation (%)"
    )

    variation_axis.set_title(
        "B   Range-based variation index",
        loc="left",
        fontweight="bold",
    )

    variation_axis.grid(
        axis="y",
        linestyle="--",
        linewidth=0.5,
        alpha=0.6,
        zorder=0,
    )

    variation_axis.legend(
        fontsize=9,
    )

    figure.text(
        0.75,
        -0.02,
        "Relative variation = (maximum − minimum) / mean × 100",
        ha="center",
        va="top",
        fontsize=9,
    )

    figure.tight_layout()

    save_figure(
        figure,
        output_path,
    )


# =========================
# MAIN WORKFLOW
# =========================


def main():
    """Run the complete reference-bias plotting workflow."""
    start_time = datetime.now()

    os.makedirs(
        OUTPUT_DIR,
        exist_ok=True,
    )

    print("=" * 65, flush=True)
    print("REFERENCE-BIAS FIGURE GENERATION", flush=True)
    print(
        f"Start       : {start_time.strftime('%Y-%m-%d %H:%M:%S')}",
        flush=True,
    )
    print(f"INPUT_TSV   : {INPUT_TSV}", flush=True)
    print(f"OUTPUT_DIR  : {OUTPUT_DIR}", flush=True)
    print(f"COMMUNITY   : {COMMUNITY}", flush=True)
    print("=" * 65, flush=True)

    data_df, variant_types, references = load_reference_comparison(
        INPUT_TSV
    )

    normalized_df = calculate_normalized_composition(
        data_df,
        variant_types,
    )

    bias_df = calculate_reference_bias(
        data_df,
        variant_types,
    )

    print(
        f"\nReferences found : {len(references)}",
        flush=True,
    )
    print(
        f"Reference labels : {', '.join(references)}",
        flush=True,
    )
    print(
        f"Variant types    : {', '.join(variant_types)}",
        flush=True,
    )

    print(
        "\nInput counts:",
        flush=True,
    )
    print(
        data_df[variant_types].to_string(),
        flush=True,
    )

    print(
        "\nGenerating figures...",
        flush=True,
    )

    grouped_output = os.path.join(
        OUTPUT_DIR,
        "Fig1_grouped_barplot.pdf",
    )

    heatmap_output = os.path.join(
        OUTPUT_DIR,
        "Fig2_heatmap_normalized.pdf",
    )

    stacked_output = os.path.join(
        OUTPUT_DIR,
        "Fig3_stacked_barplot.pdf",
    )

    bias_figure_output = os.path.join(
        OUTPUT_DIR,
        "Fig4_reference_bias.pdf",
    )

    bias_table_output = os.path.join(
        OUTPUT_DIR,
        "reference_bias_summary.tsv",
    )

    plot_grouped_barplot(
        data_df,
        variant_types,
        references,
        grouped_output,
    )

    plot_normalized_heatmap(
        normalized_df,
        heatmap_output,
    )

    plot_stacked_barplot(
        normalized_df,
        variant_types,
        references,
        stacked_output,
    )

    plot_reference_bias(
        bias_df,
        bias_figure_output,
    )

    bias_df.to_csv(
        bias_table_output,
        sep="\t",
        index=True,
        index_label="variant_type",
        float_format="%.4f",
    )

    print(
        f"  [OK] Saved: {bias_table_output}",
        flush=True,
    )

    end_time = datetime.now()
    elapsed_time = end_time - start_time

    print("\n" + "=" * 65, flush=True)
    print("REFERENCE-BIAS SUMMARY", flush=True)
    print(
        bias_df[
            [
                "minimum",
                "maximum",
                "mean",
                "range",
                "relative_variation_pct",
            ]
        ]
        .round(2)
        .to_string(),
        flush=True,
    )

    print("\n" + "=" * 65, flush=True)
    print(f"References analyzed : {len(references)}", flush=True)
    print(f"Variant types       : {len(variant_types)}", flush=True)
    print(f"Figures generated   : 4", flush=True)
    print(f"Results saved to    : {OUTPUT_DIR}", flush=True)
    print(f"Total runtime       : {elapsed_time}", flush=True)
    print(
        f"End                 : {end_time.strftime('%Y-%m-%d %H:%M:%S')}",
        flush=True,
    )
    print("=" * 65, flush=True)


if __name__ == "__main__":
    main()
