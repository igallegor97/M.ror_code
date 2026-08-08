#!/usr/bin/env Rscript

# =============================================================
# pca_global.R
#
# Global PCA combining all PGGB communities.
# Displays population structure among the five genomes.
#
# Usage:
#   Rscript pca_global.R
#
# Configuration can be provided through environment variables.
# =============================================================

# =========================
# REQUIRED PACKAGES
# =========================

required_packages <- c(
  "vcfR",
  "ggplot2",
  "ggrepel"
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
  library(vcfR)
  library(ggplot2)
  library(ggrepel)
})

# =========================
# CONFIGURATION
# =========================

BASE_DIR <- Sys.getenv(
  "PGGB_BASE",
  "/Storage/data1/isabella.gallego/MAESTRIA/PGGB_results_pacbio/pggb_partitioned_results"
)

OUTPUT_DIR <- Sys.getenv(
  "OUTPUT_DIR",
  "/Storage/data1/isabella.gallego/MAESTRIA/PGGB_results_pacbio/pca_results"
)

GRAPH_HASH <- Sys.getenv(
  "GRAPH_HASH",
  "bf3285f"
)

# COMMUNITIES may be provided as a comma-separated environment variable.
# If it is not defined, communities 0 through 9 are used.
communities_env <- Sys.getenv(
  "COMMUNITIES",
  ""
)

if (nchar(communities_env) > 0) {
  COMMUNITIES <- trimws(
    strsplit(communities_env, ",")[[1]]
  )
} else {
  COMMUNITIES <- paste0(
    "community.",
    0:9
  )
}

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

start_time <- Sys.time()

cat("=============================================================\n")
cat("GLOBAL PCA\n")
cat(sprintf(
  "Start       : %s\n",
  format(start_time, "%Y-%m-%d %H:%M:%S")
))
cat(sprintf(
  "BASE_DIR    : %s\n",
  BASE_DIR
))
cat(sprintf(
  "OUTPUT_DIR  : %s\n",
  OUTPUT_DIR
))
cat(sprintf(
  "GRAPH_HASH  : %s\n",
  GRAPH_HASH
))
cat(sprintf(
  "Communities : %s\n",
  paste(COMMUNITIES, collapse = ", ")
))
cat("=============================================================\n\n")

# =========================
# PLOT SETTINGS
# =========================

SAMPLE_COLORS <- c(
  B3   = "#0072B2",
  C26  = "#E69F00",
  CO8  = "#009E73",
  CO84 = "#CC79A7",
  E7   = "#56B4E9"
)

COMMUNITY_SHAPES <- setNames(
  0:9 %% 25,
  paste0("community.", 0:9)
)

# =========================
# BUILD GLOBAL GENOTYPE MATRIX
# =========================

# Strategy:
#   1. Read the VCF file from each community.
#   2. Convert genotypes into a binary presence/absence matrix.
#   3. Remove monomorphic variants.
#   4. Concatenate variants by rows.
#
# Sample columns must be consistent across communities.

cat("[1/3] Reading VCF files and building the global matrix...\n")

all_matrices <- list()
expected_samples <- NULL

for (community in COMMUNITIES) {
  community_dir <- file.path(
    BASE_DIR,
    sprintf(
      "all_pacbio_pansn.fasta.%s.%s",
      GRAPH_HASH,
      community
    )
  )

  vcf_matches <- list.files(
    community_dir,
    pattern = "variants\\.vcf(\\.gz)?$",
    full.names = TRUE
  )

  vcf_path <- if (length(vcf_matches) > 0) {
    vcf_matches[1]
  } else {
    NA_character_
  }

  if (is.na(vcf_path) || !file.exists(vcf_path)) {
    cat(sprintf(
      "  [SKIP] VCF not found for %s\n",
      community
    ))
    next
  }

  vcf <- tryCatch(
    read.vcfR(
      vcf_path,
      verbose = FALSE
    ),
    error = function(error) {
      cat(sprintf(
        "  [ERROR] %s: %s\n",
        community,
        error$message
      ))
      NULL
    }
  )

  if (is.null(vcf)) {
    next
  }

  gt <- extract.gt(
    vcf,
    element = "GT",
    as.numeric = FALSE
  )

  if (is.null(gt) || nrow(gt) == 0) {
    cat(sprintf(
      "  [WARNING] No genotype data found for %s\n",
      community
    ))
    next
  }

  # Verify sample consistency across communities
  sample_names <- colnames(gt)

  if (is.null(expected_samples)) {
    expected_samples <- sample_names
  } else if (!identical(
    sort(expected_samples),
    sort(sample_names)
  )) {
    cat(sprintf(
      "  [WARNING] %s contains a different sample set — keeping common samples\n",
      community
    ))

    common_samples <- intersect(
      expected_samples,
      sample_names
    )

    gt <- gt[
      ,
      common_samples,
      drop = FALSE
    ]
  }

  # Convert genotypes to a binary matrix:
  #   0 = missing or reference genotype
  #   1 = any non-reference genotype
  gt_binary <- apply(
    gt,
    2,
    function(genotype_column) {
      ifelse(
        is.na(genotype_column)
          | genotype_column == "."
          | genotype_column == "0",
        0L,
        1L
      )
    }
  )

  # Ensure that the result remains a matrix
  gt_binary <- as.matrix(gt_binary)

  # Remove monomorphic variants
  keep_variants <- apply(
    gt_binary,
    1,
    function(variant_row) {
      var(variant_row) > 0
    }
  )

  gt_binary <- gt_binary[
    keep_variants,
    ,
    drop = FALSE
  ]

  if (nrow(gt_binary) > 0) {
    # Add the community name to each row for traceability
    rownames(gt_binary) <- paste0(
      community,
      "_",
      seq_len(nrow(gt_binary))
    )

    all_matrices[[community]] <- gt_binary

    cat(sprintf(
      "  [OK] %s: %d polymorphic variants, %d samples\n",
      community,
      nrow(gt_binary),
      ncol(gt_binary)
    ))
  } else {
    cat(sprintf(
      "  [WARNING] No polymorphic variants retained for %s\n",
      community
    ))
  }
}

if (length(all_matrices) == 0) {
  stop(
    "[ERROR] No genotype matrices could be built. Check the input paths."
  )
}

# Identify samples shared by all community matrices
all_columns <- Reduce(
  intersect,
  lapply(
    all_matrices,
    colnames
  )
)

if (length(all_columns) == 0) {
  stop(
    "[ERROR] No samples are shared across all community matrices."
  )
}

cat(sprintf(
  "\n  Samples shared across all communities: %s\n",
  paste(all_columns, collapse = ", ")
))

# Ensure that every matrix has the same sample columns and column order
gt_global <- do.call(
  rbind,
  lapply(
    all_matrices,
    function(matrix_data) {
      matrix_data[
        ,
        all_columns,
        drop = FALSE
      ]
    }
  )
)

cat(sprintf(
  "\n  Global matrix: %d variants × %d samples\n\n",
  nrow(gt_global),
  ncol(gt_global)
))

# =========================
# GLOBAL PCA
# =========================

cat("[2/3] Running global PCA...\n")

# PCA requires samples as rows and variants as columns
pca_matrix <- t(gt_global)

pca <- prcomp(
  pca_matrix,
  center = TRUE,
  scale. = FALSE
)

variance_explained <- (
  summary(pca)$importance[2, ] * 100
)

cat(sprintf(
  "  PC1: %.2f%%  PC2: %.2f%%  PC3: %.2f%%\n",
  variance_explained[1],
  variance_explained[2],
  variance_explained[3]
))

# =========================
# SAVE PCA TABLES
# =========================

# PCA coordinates
scores <- as.data.frame(
  pca$x[, 1:3, drop = FALSE]
)

scores$sample <- rownames(scores)

# PGGB VCF sample columns are expected to already use clean names:
# C26, CO8, CO84, B3, and E7
scores$sample_clean <- scores$sample

coordinates_output <- file.path(
  OUTPUT_DIR,
  "global_pca_coordinates.tsv"
)

write.table(
  scores,
  coordinates_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(sprintf(
  "  PCA coordinates saved to: %s\n",
  coordinates_output
))

# Variance explained
number_of_pcs <- min(
  10,
  length(variance_explained)
)

variance_df <- data.frame(
  PC = names(
    variance_explained[
      seq_len(number_of_pcs)
    ]
  ),
  pct_var = round(
    variance_explained[
      seq_len(number_of_pcs)
    ],
    3
  ),
  cum_var = round(
    cumsum(
      variance_explained[
        seq_len(number_of_pcs)
      ]
    ),
    3
  )
)

variance_output <- file.path(
  OUTPUT_DIR,
  "global_pca_variance_explained.tsv"
)

write.table(
  variance_df,
  variance_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(sprintf(
  "  Variance explained saved to: %s\n",
  variance_output
))

# =========================
# GENERATE PLOTS
# =========================

cat("[3/3] Generating plots...\n")

pc1_label <- sprintf(
  "PC1 (%.1f%%)",
  variance_explained[1]
)

pc2_label <- sprintf(
  "PC2 (%.1f%%)",
  variance_explained[2]
)

pc3_label <- sprintf(
  "PC3 (%.1f%%)",
  variance_explained[3]
)

# Use predefined colors when all sample names are recognized.
# Otherwise, generate a fallback color palette.
sample_colors <- SAMPLE_COLORS[
  scores$sample_clean
]

if (any(is.na(sample_colors))) {
  unique_samples <- unique(
    scores$sample_clean
  )

  sample_colors <- setNames(
    colorRampPalette(
      c(
        "#0072B2",
        "#E69F00",
        "#009E73",
        "#CC79A7",
        "#56B4E9"
      )
    )(
      length(unique_samples)
    ),
    unique_samples
  )
}

# -------------------------------------------------------------
# Plot A: PC1 versus PC2
# -------------------------------------------------------------

plot_pc1_pc2 <- ggplot(
  scores,
  aes(
    x = PC1,
    y = PC2,
    color = sample_clean,
    label = sample_clean
  )
) +
  geom_point(
    size = 5,
    alpha = 0.9
  ) +
  geom_text_repel(
    size = 4,
    fontface = "bold",
    box.padding = 0.5,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = sample_colors
  ) +
  labs(
    title = "Global PCA — all PGGB communities",
    subtitle = "Genotype matrix combining variants from all communities",
    x = pc1_label,
    y = pc2_label,
    color = "Sample"
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
      size = 10,
      color = "grey40"
    ),
    panel.grid.major = element_line(
      color = "grey92",
      linewidth = 0.4
    ),
    legend.position = "right"
  )

pc1_pc2_output <- file.path(
  OUTPUT_DIR,
  "global_pca_PC1vsPC2.pdf"
)

ggsave(
  pc1_pc2_output,
  plot = plot_pc1_pc2,
  width = 7,
  height = 5.5,
  dpi = 300,
  device = "pdf"
)

cat(sprintf(
  "  PC1 versus PC2 plot: %s\n",
  pc1_pc2_output
))

# -------------------------------------------------------------
# Plot B: PC1 versus PC3
# -------------------------------------------------------------

plot_pc1_pc3 <- ggplot(
  scores,
  aes(
    x = PC1,
    y = PC3,
    color = sample_clean,
    label = sample_clean
  )
) +
  geom_point(
    size = 5,
    alpha = 0.9
  ) +
  geom_text_repel(
    size = 4,
    fontface = "bold",
    box.padding = 0.5,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = sample_colors
  ) +
  labs(
    title = "Global PCA — all PGGB communities",
    x = pc1_label,
    y = pc3_label,
    color = "Sample"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 13
    ),
    panel.grid.major = element_line(
      color = "grey92",
      linewidth = 0.4
    )
  )

pc1_pc3_output <- file.path(
  OUTPUT_DIR,
  "global_pca_PC1vsPC3.pdf"
)

ggsave(
  pc1_pc3_output,
  plot = plot_pc1_pc3,
  width = 7,
  height = 5.5,
  dpi = 300,
  device = "pdf"
)

cat(sprintf(
  "  PC1 versus PC3 plot: %s\n",
  pc1_pc3_output
))

# -------------------------------------------------------------
# Plot C: Scree plot
# -------------------------------------------------------------

scree_df <- data.frame(
  PC = factor(
    paste0(
      "PC",
      seq_len(number_of_pcs)
    ),
    levels = paste0(
      "PC",
      seq_len(number_of_pcs)
    )
  ),
  pct_var = variance_explained[
    seq_len(number_of_pcs)
  ]
)

scree_plot <- ggplot(
  scree_df,
  aes(
    x = PC,
    y = pct_var
  )
) +
  geom_col(
    fill = "#0072B2",
    width = 0.6
  ) +
  geom_line(
    aes(group = 1),
    color = "grey30",
    linewidth = 0.7
  ) +
  geom_point(
    size = 2.5,
    color = "grey20"
  ) +
  labs(
    title = "Scree plot — Global PCA",
    x = "Principal component",
    y = "Variance explained (%)"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 13
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

scree_output <- file.path(
  OUTPUT_DIR,
  "global_pca_screeplot.pdf"
)

ggsave(
  scree_output,
  plot = scree_plot,
  width = 6,
  height = 4,
  dpi = 300,
  device = "pdf"
)

cat(sprintf(
  "  Scree plot: %s\n",
  scree_output
))

# =========================
# FINAL REPORT
# =========================

end_time <- Sys.time()
elapsed_time <- end_time - start_time

cat("\n=============================================================\n")
cat("GLOBAL PCA SUMMARY\n")
cat(sprintf(
  "Variants used    : %d\n",
  nrow(gt_global)
))
cat(sprintf(
  "Samples          : %s\n",
  paste(all_columns, collapse = ", ")
))
cat(sprintf(
  "Communities used : %d\n",
  length(all_matrices)
))
cat(sprintf(
  "PC1 variance     : %.2f%%\n",
  variance_explained[1]
))
cat(sprintf(
  "PC2 variance     : %.2f%%\n",
  variance_explained[2]
))
cat(sprintf(
  "PC3 variance     : %.2f%%\n",
  variance_explained[3]
))

cat("\nFirst principal components:\n")
print(
  head(
    variance_df,
    5
  ),
  row.names = FALSE
)

cat(sprintf(
  "\nResults saved to : %s\n",
  OUTPUT_DIR
))
cat(sprintf(
  "Total runtime    : %s\n",
  elapsed_time
))
cat(sprintf(
  "End              : %s\n",
  format(end_time, "%Y-%m-%d %H:%M:%S")
))
cat("=============================================================\n")
