#!/usr/bin/env Rscript

# =============================================================
# compare_pggb_cactus_pca.R
#
# Creates a faceted comparison of PGGB and minigraph-cactus PCA
# coordinates. PCA axes are calculated independently for each method.
#
# Color = geographic region
# Label = sample identifier
#
# Usage:
#   Rscript compare_pggb_cactus_pca.R \
#       PGGB_COORDINATES.tsv \
#       CACTUS_COORDINATES.tsv \
#       OUTPUT_PDF
# =============================================================

required_packages <- c(
  "ggplot2",
  "ggrepel",
  "scales"
)

for (package_name in required_packages) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    install.packages(
      package_name,
      repos = "https://cloud.r-project.org",
      quiet = TRUE
    )
  }
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(scales)
})

arguments <- commandArgs(
  trailingOnly = TRUE
)

if (length(arguments) != 3) {
  stop(
    paste(
      "Usage:",
      "Rscript compare_pggb_cactus_pca.R",
      "PGGB_COORDINATES.tsv CACTUS_COORDINATES.tsv OUTPUT_PDF"
    )
  )
}

pggb_file <- arguments[1]
cactus_file <- arguments[2]
output_pdf <- arguments[3]

pggb <- read.delim(
  pggb_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cactus <- read.delim(
  cactus_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

pggb$method <- "PGGB"
cactus$method <- "minigraph-cactus"

required_columns <- c(
  "sample_id",
  "region",
  "PC1",
  "PC2"
)

for (table_name in c("pggb", "cactus")) {
  current_table <- get(table_name)
  missing <- setdiff(
    required_columns,
    colnames(current_table)
  )

  if (length(missing) > 0) {
    stop(
      sprintf(
        "%s coordinates are missing columns: %s",
        table_name,
        paste(missing, collapse = ", ")
      )
    )
  }
}

combined <- rbind(
  pggb[, c(required_columns, "method")],
  cactus[, c(required_columns, "method")]
)

regions <- sort(
  unique(combined$region)
)

region_colors <- setNames(
  scales::hue_pal()(
    length(regions)
  ),
  regions
)

comparison_plot <- ggplot(
  combined,
  aes(
    x = PC1,
    y = PC2,
    color = region,
    label = sample_id
  )
) +
  geom_point(
    size = 4.5,
    alpha = 0.9
  ) +
  geom_text_repel(
    size = 3.8,
    fontface = "bold",
    box.padding = 0.45,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = region_colors
  ) +
  facet_wrap(
    ~ method,
    scales = "free"
  ) +
  labs(
    title = "Population structure from all retained SNPs",
    subtitle = "PCA was calculated independently for each graph method",
    x = "PC1",
    y = "PC2",
    color = "Geographic region"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 13
    ),
    plot.subtitle = element_text(
      color = "grey40",
      size = 10
    ),
    strip.background = element_rect(
      fill = "grey95",
      color = "grey70"
    ),
    strip.text = element_text(
      face = "bold"
    ),
    panel.grid.major = element_line(
      color = "grey92",
      linewidth = 0.4
    ),
    legend.position = "bottom"
  )

ggsave(
  output_pdf,
  plot = comparison_plot,
  width = 11,
  height = 5.5,
  dpi = 300,
  device = "pdf"
)

cat(sprintf(
  "Comparison plot saved to: %s\n",
  output_pdf
))
