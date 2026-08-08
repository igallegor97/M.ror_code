#!/usr/bin/env Rscript

# =============================================================================
# plot_community_geographic_contribution.R
#
# Produces the integrated figure answering:
# "Which genomic regions explain the geographic structure of M. roreri?"
#
# Panels:
#   A. Community ranking by MRM geographic effect
#   B. Mantel r and MRM R2 heatmap
#   C. Chromosome/community ideogram colored by geographic contribution
#   D. SNP contribution versus geographic-effect size
#
# Usage:
# Rscript plot_community_geographic_contribution.R \
#   RANKING VALIDATED_MAP OUTPUT_DIR PROFILE
# =============================================================================

required_packages <- c("ggplot2", "scales", "gridExtra")
for (package_name in required_packages) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    stop(sprintf("Required R package not installed: %s", package_name))
  }
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(gridExtra)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript plot_community_geographic_contribution.R RANKING MAP OUTPUT PROFILE")
}

RANKING_FILE <- args[1]
MAP_FILE <- args[2]
OUTPUT_DIR <- args[3]
PROFILE <- args[4]
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

ranking <- read.delim(RANKING_FILE, stringsAsFactors = FALSE, check.names = FALSE)
mapping <- read.delim(MAP_FILE, stringsAsFactors = FALSE, check.names = FALSE)

data <- merge(
  ranking,
  mapping,
  by = "community",
  all.x = TRUE,
  suffixes = c("", "_map"),
  sort = FALSE
)

if ("chromosome_map" %in% colnames(data)) {
  data$chromosome <- ifelse(
    is.na(data$chromosome) | data$chromosome == "",
    data$chromosome_map,
    data$chromosome
  )
}

data$chromosome_display <- ifelse(
  is.na(data$chromosome) | data$chromosome == "",
  paste0("Unmapped: ", data$community),
  data$chromosome
)

data$community_label <- paste0(
  data$community,
  " (",
  format(data$retained_snps, big.mark = ","),
  " SNPs)"
)

data$community_label <- factor(
  data$community_label,
  levels = rev(data$community_label[order(data$mrm_geographic_beta)])
)

evidence_palette <- c(
  "Tier 1: concordant FDR-supported geographic signal" = "#0072B2",
  "Tier 2: partially FDR-supported geographic signal" = "#E69F00",
  "Tier 3: positive exploratory signal" = "#009E73",
  "No positive geographic signal" = "grey65"
)

p_rank <- ggplot(
  data,
  aes(
    x = community_label,
    y = mrm_geographic_beta,
    fill = evidence_tier
  )
) +
  geom_col(color = "grey20", linewidth = 0.25) +
  coord_flip() +
  geom_hline(yintercept = 0, linewidth = 0.45) +
  scale_fill_manual(values = evidence_palette, drop = FALSE) +
  labs(
    title = "A. Community-specific geographic effect",
    x = NULL,
    y = "MRM geographic coefficient",
    fill = "Evidence tier"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

heat_long <- rbind(
  data.frame(
    community = data$community,
    metric = "Mantel Pearson r",
    value = data$mantel_pearson_r
  ),
  data.frame(
    community = data$community,
    metric = "MRM R-squared",
    value = data$mrm_r_squared
  ),
  data.frame(
    community = data$community,
    metric = "MRM geographic beta",
    value = data$mrm_geographic_beta
  )
)

heat_long$community <- factor(
  heat_long$community,
  levels = data$community[order(data$rank)]
)

p_heat <- ggplot(
  heat_long,
  aes(x = metric, y = community, fill = value)
) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f", value)), size = 3) +
  scale_fill_gradient2(
    low = "#0072B2",
    mid = "white",
    high = "#D55E00",
    midpoint = 0,
    name = "Effect"
  ) +
  labs(
    title = "B. Geographic signal by community",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

# Ideogram mode:
# - if start/end/length are available, draw true segments;
# - otherwise draw one normalized bar per mapped chromosome/community.
numeric_fields <- c("chromosome_length_bp", "start_bp", "end_bp")
has_coordinates <- all(numeric_fields %in% colnames(data)) &&
  all(
    !is.na(suppressWarnings(as.numeric(data$chromosome_length_bp))) &
    !is.na(suppressWarnings(as.numeric(data$start_bp))) &
    !is.na(suppressWarnings(as.numeric(data$end_bp)))
  )

if (has_coordinates) {
  data$chromosome_length_bp <- as.numeric(data$chromosome_length_bp)
  data$start_bp <- as.numeric(data$start_bp)
  data$end_bp <- as.numeric(data$end_bp)
  data$start_fraction <- data$start_bp / data$chromosome_length_bp
  data$end_fraction <- data$end_bp / data$chromosome_length_bp
} else {
  data$start_fraction <- 0
  data$end_fraction <- 1
}

chrom_order <- unique(data$chromosome_display)
data$chromosome_display <- factor(
  data$chromosome_display,
  levels = rev(chrom_order)
)

p_ideo <- ggplot(
  data,
  aes(
    y = chromosome_display,
    xmin = start_fraction,
    xmax = end_fraction,
    ymin = as.numeric(chromosome_display) - 0.32,
    ymax = as.numeric(chromosome_display) + 0.32,
    fill = mrm_geographic_beta
  )
) +
  geom_rect(color = "grey20", linewidth = 0.4) +
  geom_text(
    aes(
      x = (start_fraction + end_fraction) / 2,
      label = community
    ),
    size = 3,
    color = "black"
  ) +
  scale_fill_gradient2(
    low = "#0072B2",
    mid = "white",
    high = "#D55E00",
    midpoint = 0,
    name = "MRM geographic beta"
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = percent_format(accuracy = 1)
  ) +
  labs(
    title = "C. Chromosomal distribution of geographic signal",
    x = ifelse(
      has_coordinates,
      "Relative chromosome position",
      "Normalized chromosome span"
    ),
    y = "Chromosome / mapped region"
  ) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

p_snp <- ggplot(
  data,
  aes(
    x = retained_snp_pct_among_analyzed,
    y = mrm_geographic_beta,
    size = mrm_r_squared,
    color = evidence_tier,
    label = community
  )
) +
  geom_point(alpha = 0.85) +
  geom_text(
    nudge_y = 0.003,
    size = 3,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  scale_color_manual(values = evidence_palette, drop = FALSE) +
  scale_size_continuous(
    name = "MRM R-squared",
    range = c(3, 9)
  ) +
  labs(
    title = "D. SNP abundance versus geographic effect",
    x = "Retained SNP contribution among analyzed communities (%)",
    y = "MRM geographic coefficient",
    color = "Evidence tier"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

pdf(
  file.path(
    OUTPUT_DIR,
    paste0(PROFILE, "_integrated_genomic_geography_figure.pdf")
  ),
  width = 14,
  height = 11
)

grid.arrange(
  p_rank,
  p_heat,
  p_ideo,
  p_snp,
  ncol = 2,
  top = paste(
    "Genomic regions contributing to geographic structure -",
    PROFILE
  )
)

dev.off()

ggsave(
  file.path(
    OUTPUT_DIR,
    paste0(PROFILE, "_community_ranking.pdf")
  ),
  p_rank,
  width = 8,
  height = 5.5,
  dpi = 300,
  device = "pdf"
)

ggsave(
  file.path(
    OUTPUT_DIR,
    paste0(PROFILE, "_chromosome_contribution_ideogram.pdf")
  ),
  p_ideo,
  width = 9,
  height = 6,
  dpi = 300,
  device = "pdf"
)

cat("Community/chromosome contribution figures completed.\n")
