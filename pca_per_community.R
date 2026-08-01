#!/usr/bin/env Rscript

# =============================================================
# pca_per_community.R
#
# PCA per PGGB community.
# Each community produces one PCA using the five genomes as samples.
#
# Usage:
#   Rscript pca_per_community.R
#
# Dependencies:
#   vcfR, ggplot2, ggrepel, and scales
#   Missing packages are installed automatically.
# =============================================================

# =========================
# REQUIRED PACKAGES
# =========================

required_packages <- c(
  "vcfR",
  "ggplot2",
  "ggrepel",
  "scales"
)

for (package_name in required_packages) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    message(
      paste0(
        "[INFO] Installing package: ",
        package_name
      )
    )

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
  library(scales)
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

# Number of principal components to save
N_PCS <- 3

# COMMUNITIES may be provided as a comma-separated environment variable.
# If it is not defined, communities 0 through 9 are used.
communities_env <- Sys.getenv(
  "COMMUNITIES",
  ""
)

if (nchar(communities_env) > 0) {
  COMMUNITIES <- trimws(
    strsplit(
      communities_env,
      ","
    )[[1]]
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
cat("PCA PER COMMUNITY\n")
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

# Accessible color palette based on Wong (2011)
SAMPLE_COLORS <- c(
  B3   = "#0072B2",
  C26  = "#E69F00",
  CO8  = "#009E73",
  CO84 = "#CC79A7",
  E7   = "#56B4E9"
)

# =========================
# FUNCTIONS
# =========================

#' Read a VCF file and build a binary genotype matrix.
#'
#' Rows represent variants and columns represent samples.
#'
#' Genotype encoding:
#'   0 = missing or reference genotype
#'   1 = at least one alternative allele is present
#'
#' @param vcf_path Path to a VCF or compressed VCF file.
#'
#' @return A binary genotype matrix or NULL.
build_genotype_matrix <- function(vcf_path) {
  vcf <- read.vcfR(
    vcf_path,
    verbose = FALSE
  )

  genotype_matrix <- extract.gt(
    vcf,
    element = "GT",
    as.numeric = FALSE
  )

  if (
    is.null(genotype_matrix)
    || nrow(genotype_matrix) == 0
  ) {
    warning(
      paste(
        "No genotype data found in:",
        vcf_path
      )
    )

    return(NULL)
  }

  # Convert genotypes to binary values.
  #
  # This preserves the behavior of the original script:
  # any value other than ".", NA, or exactly "0" is treated as
  # containing an alternative allele.
  genotype_binary <- apply(
    genotype_matrix,
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

  # Ensure the result remains a matrix
  genotype_binary <- as.matrix(
    genotype_binary
  )

  # Remove monomorphic variants because they do not contribute to PCA
  keep_variants <- apply(
    genotype_binary,
    1,
    function(variant_row) {
      var(variant_row) > 0
    }
  )

  genotype_binary <- genotype_binary[
    keep_variants,
    ,
    drop = FALSE
  ]

  cat(sprintf(
    paste0(
      "    Total variants: %d | ",
      "polymorphic variants used in PCA: %d\n"
    ),
    nrow(genotype_matrix),
    sum(keep_variants)
  ))

  return(genotype_binary)
}


#' Run PCA on a genotype matrix.
#'
#' The matrix is transposed so that samples are rows and variants are
#' columns.
#'
#' @param genotype_matrix Binary genotype matrix.
#'
#' @return A prcomp object.
run_pca <- function(genotype_matrix) {
  pca_matrix <- t(
    genotype_matrix
  )

  pca <- prcomp(
    pca_matrix,
    center = TRUE,
    scale. = FALSE
  )

  return(pca)
}


#' Convert PCA results into a coordinate table.
#'
#' @param pca A prcomp object.
#' @param community Community identifier.
#' @param n_pcs Maximum number of principal components to retain.
#'
#' @return A list containing PCA scores and variance explained.
pca_to_dataframe <- function(
    pca,
    community,
    n_pcs = 3
) {
  available_pcs <- min(
    n_pcs,
    ncol(pca$x)
  )

  scores <- as.data.frame(
    pca$x[
      ,
      seq_len(available_pcs),
      drop = FALSE
    ]
  )

  scores$sample <- rownames(
    scores
  )

  scores$community <- community

  variance_explained <- (
    summary(pca)$importance[
      2,
      seq_len(available_pcs)
    ] * 100
  )

  return(
    list(
      scores = scores,
      variance_explained = variance_explained
    )
  )
}


#' Generate a publication-style PCA plot.
#'
#' @param scores PCA coordinate table.
#' @param variance_explained Percentage of variance explained by each PC.
#' @param title Plot title.
#' @param output_path Output PDF path.
#' @param color_map Named vector containing sample colors.
#'
#' @return A ggplot object.
plot_pca <- function(
    scores,
    variance_explained,
    title,
    output_path,
    color_map = SAMPLE_COLORS
) {
  pc1_label <- sprintf(
    "PC1 (%.1f%%)",
    variance_explained[1]
  )

  pc2_label <- sprintf(
    "PC2 (%.1f%%)",
    variance_explained[2]
  )

  # Use the predefined color map when all sample names are recognized.
  # Otherwise, generate an automatic palette.
  samples_present <- unique(
    scores$sample
  )

  colors_to_use <- color_map[
    samples_present
  ]

  if (any(is.na(colors_to_use))) {
    colors_to_use <- setNames(
      scales::hue_pal()(
        length(samples_present)
      ),
      samples_present
    )
  }

  pca_plot <- ggplot(
    scores,
    aes(
      x = PC1,
      y = PC2,
      color = sample,
      label = sample
    )
  ) +
    geom_point(
      size = 4,
      alpha = 0.9
    ) +
    geom_text_repel(
      size = 3.5,
      fontface = "bold",
      box.padding = 0.4,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = colors_to_use
    ) +
    labs(
      title = title,
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
      legend.position = "right",
      panel.grid.major = element_line(
        color = "grey92",
        linewidth = 0.4
      )
    )

  ggsave(
    output_path,
    plot = pca_plot,
    width = 6,
    height = 5,
    dpi = 300,
    device = "pdf"
  )

  cat(sprintf(
    "    Plot saved to: %s\n",
    output_path
  ))

  return(pca_plot)
}

# =========================
# PROCESS EACH COMMUNITY
# =========================

all_scores <- list()
all_variance <- list()
failed_communities <- character()

for (community in COMMUNITIES) {
  cat(sprintf(
    "\n[%s]\n",
    community
  ))

  community_dir <- file.path(
    BASE_DIR,
    sprintf(
      "all_pacbio_pansn.fasta.%s.%s",
      GRAPH_HASH,
      community
    )
  )

  if (!dir.exists(community_dir)) {
    cat(sprintf(
      "  [SKIP] Community directory not found: %s\n",
      community_dir
    ))

    failed_communities <- c(
      failed_communities,
      community
    )

    next
  }

  vcf_matches <- list.files(
    community_dir,
    pattern = "variants\\.vcf(\\.gz)?$",
    full.names = TRUE
  )

  if (length(vcf_matches) == 0) {
    cat(sprintf(
      "  [SKIP] No variants.vcf or variants.vcf.gz file found in: %s\n",
      community_dir
    ))

    failed_communities <- c(
      failed_communities,
      community
    )

    next
  }

  vcf_path <- vcf_matches[1]

  cat(sprintf(
    "  VCF: %s\n",
    vcf_path
  ))

  # Build genotype matrix
  genotype_matrix <- tryCatch(
    build_genotype_matrix(
      vcf_path
    ),
    error = function(error) {
      cat(sprintf(
        "  [ERROR] %s\n",
        error$message
      ))

      NULL
    }
  )

  if (
    is.null(genotype_matrix)
    || nrow(genotype_matrix) < 2
  ) {
    cat(
      "  [SKIP] Insufficient polymorphic variants for PCA\n"
    )

    failed_communities <- c(
      failed_communities,
      community
    )

    next
  }

  # Run PCA
  pca <- tryCatch(
    run_pca(
      genotype_matrix
    ),
    error = function(error) {
      cat(sprintf(
        "  [ERROR] PCA failed: %s\n",
        error$message
      ))

      NULL
    }
  )

  if (is.null(pca)) {
    failed_communities <- c(
      failed_communities,
      community
    )

    next
  }

  result <- pca_to_dataframe(
    pca,
    community,
    n_pcs = N_PCS
  )

  scores <- result$scores
  variance_explained <- result$variance_explained

  if (length(variance_explained) < 2) {
    cat(
      "  [SKIP] PCA produced fewer than two principal components\n"
    )

    failed_communities <- c(
      failed_communities,
      community
    )

    next
  }

  cat(sprintf(
    "    Variance explained: PC1=%.1f%% PC2=%.1f%%\n",
    variance_explained[1],
    variance_explained[2]
  ))

  # Save PCA coordinates
  coordinates_output <- file.path(
    OUTPUT_DIR,
    sprintf(
      "%s_pca_coordinates.tsv",
      community
    )
  )

  write.table(
    scores,
    coordinates_output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  cat(sprintf(
    "    Coordinates saved to: %s\n",
    coordinates_output
  ))

  # Generate PCA plot
  plot_output <- file.path(
    OUTPUT_DIR,
    sprintf(
      "%s_pca.pdf",
      community
    )
  )

  plot_pca(
    scores = scores,
    variance_explained = variance_explained,
    title = sprintf(
      "PCA — %s (PGGB)",
      community
    ),
    output_path = plot_output
  )

  # Save variance explained
  variance_output <- file.path(
    OUTPUT_DIR,
    sprintf(
      "%s_variance_explained.tsv",
      community
    )
  )

  variance_df <- data.frame(
    PC = names(
      variance_explained
    ),
    pct_var = round(
      variance_explained,
      3
    )
  )

  write.table(
    variance_df,
    variance_output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  cat(sprintf(
    "    Variance explained saved to: %s\n",
    variance_output
  ))

  all_scores[[community]] <- scores
  all_variance[[community]] <- variance_explained
}

# =========================
# MULTI-COMMUNITY PANEL
# =========================

if (length(all_scores) > 1) {
  cat("\n[Multi-community panel]\n")

  # Combine previously calculated PCA coordinates.
  #
  # Each facet uses free scales because PCA coordinates are calculated
  # independently for each community and are therefore not directly
  # comparable in absolute scale.
  combined_scores <- do.call(
    rbind,
    all_scores
  )

  panel_plot <- ggplot(
    combined_scores,
    aes(
      x = PC1,
      y = PC2,
      color = sample,
      label = sample
    )
  ) +
    geom_point(
      size = 2.5,
      alpha = 0.85
    ) +
    geom_text_repel(
      size = 2.8,
      show.legend = FALSE,
      box.padding = 0.3,
      max.overlaps = 5
    ) +
    scale_color_manual(
      values = SAMPLE_COLORS,
      na.value = "grey50"
    ) +
    facet_wrap(
      ~ community,
      scales = "free"
    ) +
    labs(
      title = "PCA per PGGB community",
      subtitle = "PCA was calculated independently for each community",
      x = "PC1",
      y = "PC2",
      color = "Sample"
    ) +
    theme_classic(
      base_size = 10
    ) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 12
      ),
      plot.subtitle = element_text(
        size = 9,
        color = "grey40"
      ),
      strip.background = element_rect(
        fill = "grey95",
        color = "grey70"
      ),
      strip.text = element_text(
        face = "bold",
        size = 9
      ),
      panel.grid.major = element_line(
        color = "grey92",
        linewidth = 0.3
      ),
      legend.position = "bottom"
    )

  panel_output <- file.path(
    OUTPUT_DIR,
    "all_communities_pca_panel.pdf"
  )

  ggsave(
    panel_output,
    plot = panel_plot,
    width = 14,
    height = 10,
    dpi = 300,
    device = "pdf"
  )

  cat(sprintf(
    "  Panel saved to: %s\n",
    panel_output
  ))
}

# =========================
# FINAL REPORT
# =========================

end_time <- Sys.time()
elapsed_time <- end_time - start_time

cat("\n=============================================================\n")
cat("PCA PER COMMUNITY SUMMARY\n")
cat(sprintf(
  "Communities requested : %d\n",
  length(COMMUNITIES)
))
cat(sprintf(
  "Communities processed : %d\n",
  length(all_scores)
))
cat(sprintf(
  "Communities failed    : %d",
  length(failed_communities)
))

if (length(failed_communities) > 0) {
  cat(sprintf(
    "  %s\n",
    paste(
      failed_communities,
      collapse = ", "
    )
  ))
} else {
  cat("\n")
}

cat(sprintf(
  "Results saved to      : %s\n",
  OUTPUT_DIR
))
cat(sprintf(
  "Total runtime         : %s\n",
  elapsed_time
))
cat(sprintf(
  "End                   : %s\n",
  format(
    end_time,
    "%Y-%m-%d %H:%M:%S"
  )
))
cat("=============================================================\n")
