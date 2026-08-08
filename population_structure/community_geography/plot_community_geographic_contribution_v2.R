#!/usr/bin/env Rscript

# =============================================================================
# plot_community_geographic_contribution_v2.R
#
# Publication-ready visualization of chromosome/community contributions to
# geographic population structure in Moniliophthora roreri.
#
# Inputs:
#   1. Community geographic ranking TSV
#   2. Complete community-to-chromosome map TSV
#   3. Output directory
#   4. Profile name: all11 or conservative9
#
# Outputs:
#   - Integrated four-panel figure (PDF and PNG)
#   - Individual publication-ready panels
#   - Plotting-data TSV used for reproducibility
#
# Important:
#   This script does not recalculate Mantel, MRM, PCA, FDR, or SNP counts.
#   It only visualizes the final corrected outputs.
# =============================================================================

required_packages <- c(
  "ggplot2",
  "scales",
  "gridExtra",
  "ggrepel"
)

for (package_name in required_packages) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    stop(sprintf("Required R package not installed: %s", package_name))
  }
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(gridExtra)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop(
    paste(
      "Usage: Rscript plot_community_geographic_contribution_v2.R",
      "RANKING COMPLETE_MAP OUTPUT_DIR PROFILE"
    )
  )
}

RANKING_FILE <- args[1]
MAP_FILE <- args[2]
OUTPUT_DIR <- args[3]
PROFILE <- args[4]

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

message(
  sprintf(
    "Graphics capabilities: PDF=%s; PNG=%s; Cairo=%s",
    capabilities("pdf"),
    capabilities("png"),
    capabilities("cairo")
  )
)

# -----------------------------------------------------------------------------
# Read and validate inputs
# -----------------------------------------------------------------------------

ranking <- read.delim(
  RANKING_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

mapping <- read.delim(
  MAP_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_ranking_columns <- c(
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
  "retained_snp_pct_among_analyzed"
)

required_mapping_columns <- c(
  "community",
  "chromosome",
  "syntenic_group",
  "mapping_type",
  "chromosome_length_bp",
  "start_bp",
  "end_bp",
  "note"
)

missing_ranking <- setdiff(
  required_ranking_columns,
  colnames(ranking)
)

missing_mapping <- setdiff(
  required_mapping_columns,
  colnames(mapping)
)

if (length(missing_ranking) > 0) {
  stop(
    paste(
      "Ranking table is missing columns:",
      paste(missing_ranking, collapse = ", ")
    )
  )
}

if (length(missing_mapping) > 0) {
  stop(
    paste(
      "Community map is missing columns:",
      paste(missing_mapping, collapse = ", ")
    )
  )
}

# -----------------------------------------------------------------------------
# Normalize data types
# -----------------------------------------------------------------------------

numeric_ranking_columns <- c(
  "retained_snps",
  "mantel_pearson_r",
  "mantel_pearson_q",
  "mrm_geographic_beta",
  "mrm_geographic_q",
  "mrm_r_squared",
  "mrm_model_q",
  "retained_snp_pct_among_analyzed"
)

for (column_name in numeric_ranking_columns) {
  ranking[[column_name]] <- suppressWarnings(
    as.numeric(ranking[[column_name]])
  )
}

numeric_mapping_columns <- c(
  "chromosome_length_bp",
  "start_bp",
  "end_bp"
)

for (column_name in numeric_mapping_columns) {
  mapping[[column_name]] <- suppressWarnings(
    as.numeric(mapping[[column_name]])
  )
}

# -----------------------------------------------------------------------------
# Build complete chromosome plotting table
# -----------------------------------------------------------------------------

ranking_statistics <- ranking[
  ,
  setdiff(
    colnames(ranking),
    c(
      "chromosome",
      "syntenic_group",
      "mapping_type",
      "chromosome_length_bp",
      "start_bp",
      "end_bp",
      "note"
    )
  ),
  drop = FALSE
]

plot_data <- merge(
  mapping,
  ranking_statistics,
  by = "community",
  all.x = TRUE,
  sort = FALSE
)

plot_data$analyzed <- !is.na(
  plot_data$mrm_geographic_beta
)

plot_data$analysis_label <- ifelse(
  plot_data$analyzed,
  "Analyzed",
  "No retained SNPs / not analyzed"
)

plot_data$retained_snps[
  is.na(plot_data$retained_snps)
] <- 0

plot_data$retained_snp_pct_among_analyzed[
  is.na(plot_data$retained_snp_pct_among_analyzed)
] <- 0

plot_data$chromosome_length_mb <- (
  plot_data$chromosome_length_bp / 1e6
)

plot_data$start_mb <- plot_data$start_bp / 1e6
plot_data$end_mb <- plot_data$end_bp / 1e6

# Natural chromosome order: Group1 ... Group11.
extract_group_number <- function(value) {
  number <- suppressWarnings(
    as.integer(
      sub(
        "^.*?([0-9]+).*$",
        "\\1",
        value
      )
    )
  )

  ifelse(
    is.na(number),
    999L,
    number
  )
}

plot_data$chromosome_number <- extract_group_number(
  plot_data$chromosome
)

plot_data <- plot_data[
  order(
    plot_data$chromosome_number,
    plot_data$community
  ),
  ,
  drop = FALSE
]

chromosome_levels <- unique(
  plot_data$chromosome
)

plot_data$chromosome <- factor(
  plot_data$chromosome,
  levels = rev(chromosome_levels)
)

# Ranking-only table.
analyzed_data <- plot_data[
  plot_data$analyzed,
  ,
  drop = FALSE
]

# Short labels for publication figures.
analyzed_data$community_chromosome_label <- paste0(
  as.character(analyzed_data$chromosome),
  " / ",
  analyzed_data$community
)

analyzed_data$community_chromosome_label <- factor(
  analyzed_data$community_chromosome_label,
  levels = rev(
    analyzed_data$community_chromosome_label[
      order(analyzed_data$mrm_geographic_beta)
    ]
  )
)

# -----------------------------------------------------------------------------
# Shared style
# -----------------------------------------------------------------------------

evidence_palette <- c(
  "Tier 1: concordant FDR-supported geographic signal" = "#0072B2",
  "Tier 2: partially FDR-supported geographic signal" = "#E69F00",
  "Tier 3: positive exploratory signal" = "#009E73",
  "No positive geographic signal" = "grey65"
)

mapping_linetypes <- c(
  "chromosome_scale" = "solid",
  "mixed_assignment" = "longdash",
  "partial_or_ambiguous" = "dotted",
  "unmapped" = "dotdash"
)

profile_title <- ifelse(
  PROFILE == "all11",
  "All 11 selected communities",
  ifelse(
    PROFILE == "conservative9",
    "Conservative 9 selected communities",
    PROFILE
  )
)

theme_publication <- theme_classic(
  base_size = 11
) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 12.5,
      margin = margin(b = 4)
    ),
    plot.subtitle = element_text(
      size = 9.5,
      color = "grey30",
      margin = margin(b = 7)
    ),
    axis.title = element_text(size = 10.5),
    axis.text = element_text(size = 9.3),
    legend.title = element_text(
      face = "bold",
      size = 9.5
    ),
    legend.text = element_text(size = 8.6),
    plot.margin = margin(8, 10, 8, 8)
  )

# -----------------------------------------------------------------------------
# Panel A: community/chromosome MRM effect ranking
# -----------------------------------------------------------------------------

panel_a <- ggplot(
  analyzed_data,
  aes(
    x = community_chromosome_label,
    y = mrm_geographic_beta,
    fill = evidence_tier
  )
) +
  geom_col(
    width = 0.72,
    color = "grey20",
    linewidth = 0.3
  ) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.45
  ) +
  geom_text(
    aes(
      label = sprintf(
        "beta = %.3f",
        mrm_geographic_beta
      )
    ),
    hjust = -0.08,
    size = 3
  ) +
  coord_flip(
    clip = "off"
  ) +
  scale_fill_manual(
    values = evidence_palette,
    drop = FALSE,
    name = "Evidence tier"
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(0.02, 0.20)
    )
  ) +
  labs(
    title = "A. Chromosome-specific geographic effect",
    subtitle = paste(
      "MRM geographic coefficient after controlling",
      "for sequencing technology."
    ),
    x = NULL,
    y = "MRM geographic coefficient"
  ) +
  theme_publication +
  theme(
    legend.position = "bottom",
    legend.direction = "vertical",
    plot.margin = margin(8, 22, 8, 8)
  )

# -----------------------------------------------------------------------------
# Panel B: clean effect-size heatmap with FDR markers
# -----------------------------------------------------------------------------

heat_data <- rbind(
  data.frame(
    chromosome = analyzed_data$chromosome,
    community = analyzed_data$community,
    metric = "Mantel r",
    value = analyzed_data$mantel_pearson_r,
    q_value = analyzed_data$mantel_pearson_q
  ),
  data.frame(
    chromosome = analyzed_data$chromosome,
    community = analyzed_data$community,
    metric = "MRM beta",
    value = analyzed_data$mrm_geographic_beta,
    q_value = analyzed_data$mrm_geographic_q
  ),
  data.frame(
    chromosome = analyzed_data$chromosome,
    community = analyzed_data$community,
    metric = "MRM R2",
    value = analyzed_data$mrm_r_squared,
    q_value = analyzed_data$mrm_model_q
  )
)

heat_data$row_label <- paste0(
  as.character(heat_data$chromosome),
  " / ",
  heat_data$community
)

rank_order <- analyzed_data$community[
  order(analyzed_data$rank)
]

heat_levels <- paste0(
  as.character(
    analyzed_data$chromosome[
      match(
        rank_order,
        analyzed_data$community
      )
    ]
  ),
  " / ",
  rank_order
)

heat_data$row_label <- factor(
  heat_data$row_label,
  levels = rev(heat_levels)
)

heat_data$significance_marker <- ifelse(
  is.na(heat_data$q_value),
  "",
  ifelse(
    heat_data$q_value <= 0.01,
    "**",
    ifelse(
      heat_data$q_value <= 0.05,
      "*",
      ""
    )
  )
)

maximum_absolute_effect <- max(
  abs(heat_data$value),
  na.rm = TRUE
)

panel_b <- ggplot(
  heat_data,
  aes(
    x = metric,
    y = row_label,
    fill = value
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.6
  ) +
  geom_text(
    aes(label = significance_marker),
    size = 5,
    fontface = "bold",
    color = "black"
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(
      -maximum_absolute_effect,
      maximum_absolute_effect
    ),
    name = "Effect size"
  ) +
  labs(
    title = "B. Concordance of geographic-effect metrics",
    subtitle = "* FDR q <= 0.05; ** FDR q <= 0.01.",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(
    base_size = 11
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 12.5
    ),
    plot.subtitle = element_text(
      size = 9.5,
      color = "grey30"
    ),
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 20,
      hjust = 1
    ),
    axis.text.y = element_text(size = 9),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(8, 10, 8, 8)
  )

# -----------------------------------------------------------------------------
# Panel C: chromosome ideogram using real chromosome lengths
# -----------------------------------------------------------------------------

# Background represents the full chromosome.
# Colored overlay represents the mapped community segment.
# In the current map, chromosome-scale communities span the full chromosome,
# while mixed/ambiguous mappings retain a distinct border style.

plot_data$chromosome_numeric <- as.numeric(
  plot_data$chromosome
)

panel_c <- ggplot() +
  geom_rect(
    data = plot_data,
    aes(
      xmin = 0,
      xmax = chromosome_length_mb,
      ymin = chromosome_numeric - 0.30,
      ymax = chromosome_numeric + 0.30
    ),
    fill = "grey93",
    color = "grey45",
    linewidth = 0.35
  ) +
  geom_rect(
    data = plot_data[
      plot_data$analyzed,
      ,
      drop = FALSE
    ],
    aes(
      xmin = start_mb,
      xmax = end_mb,
      ymin = chromosome_numeric - 0.30,
      ymax = chromosome_numeric + 0.30,
      fill = mrm_geographic_beta,
      linetype = mapping_type
    ),
    color = "black",
    linewidth = 0.75
  ) +
  geom_text(
    data = plot_data,
    aes(
      x = chromosome_length_mb + 0.15,
      y = chromosome_numeric,
      label = ifelse(
        analyzed,
        paste0(
          community,
          " | ",
          comma(retained_snps),
          " SNPs"
        ),
        paste0(
          community,
          " | no retained SNPs"
        )
      )
    ),
    hjust = 0,
    size = 2.75
  ) +
  scale_y_continuous(
    breaks = seq_along(
      levels(plot_data$chromosome)
    ),
    labels = levels(plot_data$chromosome),
    expand = expansion(add = 0.6)
  ) +
  scale_x_continuous(
    labels = label_number(
      accuracy = 1,
      suffix = " Mb"
    ),
    expand = expansion(
      mult = c(0.01, 0.32)
    )
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = "MRM geographic beta",
    na.value = "grey85"
  ) +
  scale_linetype_manual(
    values = mapping_linetypes,
    name = "Mapping type",
    drop = FALSE
  ) +
  labs(
    title = "C. Chromosomal distribution of geographic structure",
    subtitle = paste(
      "Bars are proportional to C26 chromosome length.",
      "Grey chromosomes had no retained SNPs in the final profile."
    ),
    x = "Reference chromosome length",
    y = NULL
  ) +
  theme_publication +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.box = "vertical",
    plot.margin = margin(8, 90, 8, 8)
  ) +
  coord_cartesian(
    clip = "off"
  )

# -----------------------------------------------------------------------------
# Panel D: SNP abundance versus geographic effect
# -----------------------------------------------------------------------------

# Bubble size represents Mantel r, as requested.
# Negative values, if present in future analyses, are truncated to zero only
# for size mapping; their signed value remains visible on the y axis.

analyzed_data$mantel_size <- pmax(
  analyzed_data$mantel_pearson_r,
  0
)

panel_d <- ggplot(
  analyzed_data,
  aes(
    x = retained_snp_pct_among_analyzed,
    y = mrm_geographic_beta,
    size = mantel_size,
    fill = evidence_tier,
    label = paste0(
      chromosome,
      "\n",
      community
    )
  )
) +
  geom_point(
    shape = 21,
    color = "black",
    stroke = 0.75,
    alpha = 0.88
  ) +
  geom_text_repel(
    size = 2.9,
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.30,
    min.segment.length = 0,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = evidence_palette,
    drop = FALSE,
    name = "Evidence tier"
  ) +
  scale_size_continuous(
    name = "Mantel Pearson r",
    range = c(3.5, 10),
    breaks = pretty_breaks(n = 4)
  ) +
  scale_x_continuous(
    labels = label_percent(
      scale = 1,
      accuracy = 1
    ),
    expand = expansion(
      mult = c(0.05, 0.12)
    )
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(0.08, 0.20)
    )
  ) +
  labs(
    title = "D. SNP abundance does not fully determine geographic effect",
    subtitle = paste(
      "Bubble size represents Mantel correlation;",
      "color represents statistical evidence."
    ),
    x = "Retained SNP contribution among analyzed communities",
    y = "MRM geographic coefficient"
  ) +
  theme_publication +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    plot.margin = margin(8, 14, 8, 8)
  )

# -----------------------------------------------------------------------------
# Save reproducibility table
# -----------------------------------------------------------------------------

write.table(
  plot_data,
  file.path(
    OUTPUT_DIR,
    paste0(
      PROFILE,
      "_community_chromosome_plotting_data.tsv"
    )
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Save individual panels
# -----------------------------------------------------------------------------

save_panel <- function(
  plot_object,
  file_stub,
  width,
  height
) {
  pdf_path <- file.path(
    OUTPUT_DIR,
    paste0(
      PROFILE,
      "_",
      file_stub,
      ".pdf"
    )
  )

  ggsave(
    pdf_path,
    plot_object,
    width = width,
    height = height,
    dpi = 300,
    device = "pdf"
  )

  if (capabilities("png")) {
    png_path <- file.path(
      OUTPUT_DIR,
      paste0(
        PROFILE,
        "_",
        file_stub,
        ".png"
      )
    )

    ggsave(
      png_path,
      plot_object,
      width = width,
      height = height,
      dpi = 300,
      device = "png"
    )
  } else {
    message(
      sprintf(
        "[WARNING] PNG device unavailable in this R build; skipped PNG for %s.",
        file_stub
      )
    )
  }
}

save_panel(
  panel_a,
  "panel_A_chromosome_effect_ranking",
  8.5,
  5.8
)

save_panel(
  panel_b,
  "panel_B_effect_heatmap",
  7.4,
  5.8
)

save_panel(
  panel_c,
  "panel_C_chromosome_ideogram",
  10.5,
  7.0
)

save_panel(
  panel_d,
  "panel_D_snp_effect_bubble",
  8.3,
  6.3
)

# -----------------------------------------------------------------------------
# Save integrated figure
# -----------------------------------------------------------------------------

integrated_title <- paste(
  "Chromosomal contributions to geographic population structure -",
  profile_title
)

pdf(
  file.path(
    OUTPUT_DIR,
    paste0(
      PROFILE,
      "_integrated_genomic_geography_figure_v2.pdf"
    )
  ),
  width = 15.5,
  height = 11.5,
  onefile = FALSE
)

grid.arrange(
  panel_a,
  panel_b,
  panel_c,
  panel_d,
  ncol = 2,
  widths = c(1.05, 1),
  heights = c(1, 1.12),
  top = grid::textGrob(
    integrated_title,
    gp = grid::gpar(
      fontsize = 17,
      fontface = "bold"
    )
  )
)

dev.off()

if (capabilities("png")) {
  png(
    file.path(
      OUTPUT_DIR,
      paste0(
        PROFILE,
        "_integrated_genomic_geography_figure_v2.png"
      )
    ),
    width = 15.5,
    height = 11.5,
    units = "in",
    res = 300
  )

  grid.arrange(
    panel_a,
    panel_b,
    panel_c,
    panel_d,
    ncol = 2,
    widths = c(1.05, 1),
    heights = c(1, 1.12),
    top = grid::textGrob(
      integrated_title,
      gp = grid::gpar(
        fontsize = 17,
        fontface = "bold"
      )
    )
  )

  dev.off()
} else {
  message(
    paste(
      "[WARNING] This R installation has no PNG device support.",
      "The integrated PNG was skipped; the PDF was generated successfully."
    )
  )
}

cat(
  sprintf(
    "Publication figures completed for profile %s.\n",
    PROFILE
  )
)
