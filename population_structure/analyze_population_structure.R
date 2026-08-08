#!/usr/bin/env Rscript

# =============================================================
# analyze_population_structure.R
#
# Performs population-structure analyses from an imputed SNP matrix.
#
# Generated outputs:
#   1. PCA coordinates and variance explained.
#   2. PC1 versus PC2 and PC1 versus PC3 plots.
#   3. Scree plot.
#   4. Euclidean and Hamming distance matrices.
#   5. Distance heatmaps.
#   6. Hierarchical-clustering dendrograms.
#
# Plot aesthetics:
#   Color = geographic region
#   Label = sample identifier
#
# Usage:
#   Rscript analyze_population_structure.R \
#       SNP_MATRIX.tsv \
#       sample_metadata.tsv \
#       SOURCE_LABEL \
#       OUTPUT_DIR
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

# =========================
# COMMAND-LINE ARGUMENTS
# =========================

arguments <- commandArgs(
  trailingOnly = TRUE
)

if (length(arguments) != 4) {
  stop(
    paste(
      "Usage:",
      "Rscript analyze_population_structure.R",
      "SNP_MATRIX.tsv sample_metadata.tsv SOURCE_LABEL OUTPUT_DIR"
    )
  )
}

SNP_MATRIX_FILE <- arguments[1]
METADATA_FILE <- arguments[2]
SOURCE_LABEL <- arguments[3]
OUTPUT_DIR <- arguments[4]

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

start_time <- Sys.time()

cat("=============================================================\n")
cat("POPULATION-STRUCTURE ANALYSIS\n")
cat(sprintf("Start       : %s\n", format(start_time, "%Y-%m-%d %H:%M:%S")))
cat(sprintf("Source      : %s\n", SOURCE_LABEL))
cat(sprintf("SNP matrix  : %s\n", SNP_MATRIX_FILE))
cat(sprintf("Metadata    : %s\n", METADATA_FILE))
cat(sprintf("OUTPUT_DIR  : %s\n", OUTPUT_DIR))
cat("=============================================================\n\n")

# =========================
# LOAD DATA
# =========================

snp_matrix <- read.delim(
  SNP_MATRIX_FILE,
  row.names = 1,
  check.names = FALSE
)

metadata <- read.delim(
  METADATA_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_metadata_columns <- c(
  "sample_id",
  "region"
)

missing_columns <- setdiff(
  required_metadata_columns,
  colnames(metadata)
)

if (length(missing_columns) > 0) {
  stop(
    sprintf(
      "Metadata is missing required columns: %s",
      paste(missing_columns, collapse = ", ")
    )
  )
}

if (anyDuplicated(metadata$sample_id)) {
  stop("Metadata contains duplicated sample_id values.")
}

sample_ids <- colnames(
  snp_matrix
)

missing_metadata <- setdiff(
  sample_ids,
  metadata$sample_id
)

if (length(missing_metadata) > 0) {
  stop(
    sprintf(
      "Metadata is missing samples: %s",
      paste(missing_metadata, collapse = ", ")
    )
  )
}

metadata <- metadata[
  match(sample_ids, metadata$sample_id),
  ,
  drop = FALSE
]

genotype_matrix <- t(
  as.matrix(snp_matrix)
)

storage.mode(genotype_matrix) <- "numeric"

if (anyNA(genotype_matrix)) {
  stop(
    "The PCA input matrix contains missing values. Use the imputed matrix."
  )
}

if (nrow(genotype_matrix) < 3) {
  stop("At least three samples are required.")
}

if (ncol(genotype_matrix) < 2) {
  stop("At least two SNP features are required.")
}

cat(sprintf(
  "Samples      : %d\nSNP features : %d\n\n",
  nrow(genotype_matrix),
  ncol(genotype_matrix)
))

# =========================
# REGION COLORS
# =========================

regions <- sort(
  unique(metadata$region)
)

region_colors <- setNames(
  scales::hue_pal()(
    length(regions)
  ),
  regions
)

write.table(
  data.frame(
    region = names(region_colors),
    color = unname(region_colors)
  ),
  file.path(
    OUTPUT_DIR,
    sprintf("%s_region_colors.tsv", SOURCE_LABEL)
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# =========================
# PCA
# =========================

cat("[1/4] Running PCA...\n")

pca <- prcomp(
  genotype_matrix,
  center = TRUE,
  scale. = FALSE
)

variance_explained <- (
  pca$sdev^2
  / sum(pca$sdev^2)
  * 100
)

number_of_axes <- min(
  10,
  length(variance_explained)
)

scores <- as.data.frame(
  pca$x[
    ,
    seq_len(min(3, ncol(pca$x))),
    drop = FALSE
  ]
)

scores$sample_id <- rownames(
  scores
)

scores <- merge(
  scores,
  metadata,
  by = "sample_id",
  sort = FALSE
)

scores <- scores[
  match(sample_ids, scores$sample_id),
  ,
  drop = FALSE
]

coordinates_output <- file.path(
  OUTPUT_DIR,
  sprintf("%s_pca_coordinates.tsv", SOURCE_LABEL)
)

write.table(
  scores,
  coordinates_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

variance_table <- data.frame(
  PC = paste0(
    "PC",
    seq_len(number_of_axes)
  ),
  pct_variance = round(
    variance_explained[
      seq_len(number_of_axes)
    ],
    4
  ),
  cumulative_variance = round(
    cumsum(
      variance_explained[
        seq_len(number_of_axes)
      ]
    ),
    4
  )
)

write.table(
  variance_table,
  file.path(
    OUTPUT_DIR,
    sprintf("%s_pca_variance_explained.tsv", SOURCE_LABEL)
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(sprintf(
  "  PC1=%.2f%%  PC2=%.2f%%",
  variance_explained[1],
  variance_explained[2]
))

if (length(variance_explained) >= 3) {
  cat(sprintf(
    "  PC3=%.2f%%",
    variance_explained[3]
  ))
}

cat("\n")

# =========================
# PCA PLOTS
# =========================

plot_pca_axes <- function(
    scores,
    x_axis,
    y_axis,
    output_file
) {
  x_index <- as.integer(
    sub("PC", "", x_axis)
  )

  y_index <- as.integer(
    sub("PC", "", y_axis)
  )

  plot_object <- ggplot(
    scores,
    aes(
      x = .data[[x_axis]],
      y = .data[[y_axis]],
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
    labs(
      title = sprintf(
        "%s SNP PCA",
        SOURCE_LABEL
      ),
      subtitle = sprintf(
        "%d SNP features across %d genomes",
        ncol(genotype_matrix),
        nrow(genotype_matrix)
      ),
      x = sprintf(
        "%s (%.1f%%)",
        x_axis,
        variance_explained[x_index]
      ),
      y = sprintf(
        "%s (%.1f%%)",
        y_axis,
        variance_explained[y_index]
      ),
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
      panel.grid.major = element_line(
        color = "grey92",
        linewidth = 0.4
      ),
      legend.position = "right"
    )

  ggsave(
    output_file,
    plot = plot_object,
    width = 7,
    height = 5.5,
    dpi = 300,
    device = "pdf"
  )
}

plot_pca_axes(
  scores,
  "PC1",
  "PC2",
  file.path(
    OUTPUT_DIR,
    sprintf("%s_pca_PC1vsPC2.pdf", SOURCE_LABEL)
  )
)

if ("PC3" %in% colnames(scores)) {
  plot_pca_axes(
    scores,
    "PC1",
    "PC3",
    file.path(
      OUTPUT_DIR,
      sprintf("%s_pca_PC1vsPC3.pdf", SOURCE_LABEL)
    )
  )
}

scree_plot <- ggplot(
  variance_table,
  aes(
    x = factor(PC, levels = PC),
    y = pct_variance,
    group = 1
  )
) +
  geom_col(
    width = 0.65
  ) +
  geom_line(
    linewidth = 0.7
  ) +
  geom_point(
    size = 2.4
  ) +
  labs(
    title = sprintf(
      "Scree plot — %s",
      SOURCE_LABEL
    ),
    x = "Principal component",
    y = "Variance explained (%)"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

ggsave(
  file.path(
    OUTPUT_DIR,
    sprintf("%s_pca_screeplot.pdf", SOURCE_LABEL)
  ),
  plot = scree_plot,
  width = 6,
  height = 4,
  dpi = 300,
  device = "pdf"
)

# =========================
# DISTANCE MATRICES
# =========================

cat("[2/4] Calculating distance matrices...\n")

euclidean_distance <- as.matrix(
  dist(
    genotype_matrix,
    method = "euclidean"
  )
)

hamming_distance <- matrix(
  0,
  nrow = nrow(genotype_matrix),
  ncol = nrow(genotype_matrix),
  dimnames = list(
    rownames(genotype_matrix),
    rownames(genotype_matrix)
  )
)

for (i in seq_len(nrow(genotype_matrix))) {
  for (j in seq_len(nrow(genotype_matrix))) {
    hamming_distance[i, j] <- mean(
      genotype_matrix[i, ] != genotype_matrix[j, ]
    )
  }
}

write.table(
  euclidean_distance,
  file.path(
    OUTPUT_DIR,
    sprintf("%s_euclidean_distance.tsv", SOURCE_LABEL)
  ),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

write.table(
  hamming_distance,
  file.path(
    OUTPUT_DIR,
    sprintf("%s_hamming_distance.tsv", SOURCE_LABEL)
  ),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# =========================
# HEATMAPS AND CLUSTERING
# =========================

cat("[3/4] Generating heatmaps and dendrograms...\n")

plot_distance_heatmap <- function(
    distance_matrix,
    title,
    output_file
) {
  sample_order <- hclust(
    as.dist(distance_matrix),
    method = "average"
  )$order

  ordered_matrix <- distance_matrix[
    sample_order,
    sample_order,
    drop = FALSE
  ]

  long_table <- as.data.frame(
    as.table(ordered_matrix)
  )

  colnames(long_table) <- c(
    "sample_x",
    "sample_y",
    "distance"
  )

  heatmap_plot <- ggplot(
    long_table,
    aes(
      x = sample_x,
      y = sample_y,
      fill = distance
    )
  ) +
    geom_tile(
      color = "white",
      linewidth = 0.4
    ) +
    geom_text(
      aes(
        label = sprintf("%.3f", distance)
      ),
      size = 3
    ) +
    scale_fill_gradient(
      low = "white",
      high = "black"
    ) +
    labs(
      title = title,
      x = NULL,
      y = NULL,
      fill = "Distance"
    ) +
    theme_classic(
      base_size = 11
    ) +
    theme(
      plot.title = element_text(
        face = "bold"
      ),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      )
    )

  ggsave(
    output_file,
    plot = heatmap_plot,
    width = 6,
    height = 5.5,
    dpi = 300,
    device = "pdf"
  )
}

plot_distance_heatmap(
  euclidean_distance,
  sprintf(
    "Euclidean genetic distance — %s",
    SOURCE_LABEL
  ),
  file.path(
    OUTPUT_DIR,
    sprintf("%s_euclidean_distance_heatmap.pdf", SOURCE_LABEL)
  )
)

plot_distance_heatmap(
  hamming_distance,
  sprintf(
    "Hamming genetic distance — %s",
    SOURCE_LABEL
  ),
  file.path(
    OUTPUT_DIR,
    sprintf("%s_hamming_distance_heatmap.pdf", SOURCE_LABEL)
  )
)

plot_dendrogram <- function(
    distance_matrix,
    title,
    output_file
) {
  clustering <- hclust(
    as.dist(distance_matrix),
    method = "average"
  )

  pdf(
    output_file,
    width = 7,
    height = 5
  )

  plot(
    clustering,
    main = title,
    xlab = "",
    sub = "",
    ylab = "Distance",
    hang = -1
  )

  dev.off()
}

plot_dendrogram(
  euclidean_distance,
  sprintf(
    "UPGMA clustering — Euclidean distance — %s",
    SOURCE_LABEL
  ),
  file.path(
    OUTPUT_DIR,
    sprintf("%s_euclidean_UPGMA_dendrogram.pdf", SOURCE_LABEL)
  )
)

plot_dendrogram(
  hamming_distance,
  sprintf(
    "UPGMA clustering — Hamming distance — %s",
    SOURCE_LABEL
  ),
  file.path(
    OUTPUT_DIR,
    sprintf("%s_hamming_UPGMA_dendrogram.pdf", SOURCE_LABEL)
  )
)

# =========================
# FINAL REPORT
# =========================

cat("[4/4] Finalizing analysis...\n")

end_time <- Sys.time()

cat("\n=============================================================\n")
cat("POPULATION-STRUCTURE SUMMARY\n")
cat(sprintf("Source          : %s\n", SOURCE_LABEL))
cat(sprintf("Samples         : %d\n", nrow(genotype_matrix)))
cat(sprintf("SNP features    : %d\n", ncol(genotype_matrix)))
cat(sprintf("Regions         : %s\n", paste(regions, collapse = ", ")))
cat(sprintf("Results saved to: %s\n", OUTPUT_DIR))
cat(sprintf("Total runtime   : %s\n", end_time - start_time))
cat(sprintf("End             : %s\n", format(end_time, "%Y-%m-%d %H:%M:%S")))
cat("=============================================================\n")
