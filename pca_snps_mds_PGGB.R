#!/usr/bin/env Rscript

# =============================================================
# pca_snps_mds.R
#
# Analysis 1: Global PCA using SNPs only.
# Analysis 2: MDS using Jaccard distance across all variant types.
#
# Both analyses use PGGB VCF files, with one VCF per community.
# The global genotype matrices are built by concatenating variants
# from all communities.
#
# Usage:
#   Rscript pca_snps_mds.R
#
# Dependencies:
#   vcfR, ggplot2, and ggrepel
#   Missing packages are installed automatically.
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
    message(
      paste0(
        "[INFO] Installing package: ",
        package_name
      )
    )

    install.packages(
      package_name,
      lib = Sys.getenv("R_LIBS_USER"),
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

COMMUNITIES <- if (nchar(communities_env) > 0) {
  trimws(
    strsplit(
      communities_env,
      ","
    )[[1]]
  )
} else {
  paste0(
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
cat("SNP PCA AND JACCARD MDS\n")
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
# BASE FUNCTIONS
# =========================

#' Find the variants VCF file inside a community directory.
#'
#' @param community_dir Path to a PGGB community directory.
#'
#' @return Path to the VCF file, or NA if no file is found.
find_vcf <- function(community_dir) {
  if (!dir.exists(community_dir)) {
    return(NA_character_)
  }

  vcf_matches <- list.files(
    community_dir,
    pattern = "variants\\.vcf(\\.gz)?$",
    full.names = TRUE
  )

  if (length(vcf_matches) == 0) {
    return(NA_character_)
  }

  return(vcf_matches[1])
}


#' Classify a variant according to REF and ALT allele lengths.
#'
#' Only the first alternative allele is considered when the record
#' contains multiple ALT alleles.
#'
#' @param ref Reference allele.
#' @param alt Alternative allele or comma-separated alleles.
#'
#' @return Variant type: SNP, MNP, INS, DEL, or COMPLEX.
classify_variant <- function(ref, alt) {
  first_alt <- sub(
    ",.*",
    "",
    alt
  )

  if (
    nchar(ref) == 1
    && nchar(first_alt) == 1
  ) {
    return("SNP")
  }

  if (nchar(ref) == nchar(first_alt)) {
    return("MNP")
  }

  if (nchar(ref) < nchar(first_alt)) {
    return("INS")
  }

  if (nchar(ref) > nchar(first_alt)) {
    return("DEL")
  }

  return("COMPLEX")
}


#' Remove monomorphic rows from a genotype matrix.
#'
#' @param genotype_matrix Binary genotype matrix.
#'
#' @return Matrix containing only polymorphic rows.
filter_polymorphic_variants <- function(genotype_matrix) {
  if (
    is.null(genotype_matrix)
    || nrow(genotype_matrix) == 0
  ) {
    return(
      genotype_matrix[
        FALSE,
        ,
        drop = FALSE
      ]
    )
  }

  keep_variants <- apply(
    genotype_matrix,
    1,
    function(variant_row) {
      var(variant_row) > 0
    }
  )

  genotype_matrix[
    keep_variants,
    ,
    drop = FALSE
  ]
}


#' Read a VCF and generate binary genotype matrices by variant type.
#'
#' Rows represent variants and columns represent samples.
#'
#' Genotype encoding:
#'   0 = missing or reference genotype
#'   1 = at least one alternative allele is present
#'
#' Returns matrices for:
#'   all, SNP, INS, DEL, and MNP
#'
#' @param vcf_path Path to a VCF or compressed VCF file.
#'
#' @return Named list of genotype matrices, or NULL.
read_vcf_by_type <- function(vcf_path) {
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
    return(NULL)
  }

  # Convert genotypes to binary values.
  #
  # This preserves the behavior of the original script:
  # values other than NA, ".", or exactly "0" are classified as
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

  genotype_binary <- as.matrix(
    genotype_binary
  )

  reference_alleles <- getREF(vcf)
  alternative_alleles <- getALT(vcf)

  variant_types <- mapply(
    classify_variant,
    reference_alleles,
    alternative_alleles,
    USE.NAMES = FALSE
  )

  all_matrix <- filter_polymorphic_variants(
    genotype_binary
  )

  snp_matrix <- filter_polymorphic_variants(
    genotype_binary[
      variant_types == "SNP",
      ,
      drop = FALSE
    ]
  )

  insertion_matrix <- filter_polymorphic_variants(
    genotype_binary[
      variant_types == "INS",
      ,
      drop = FALSE
    ]
  )

  deletion_matrix <- filter_polymorphic_variants(
    genotype_binary[
      variant_types == "DEL",
      ,
      drop = FALSE
    ]
  )

  mnp_matrix <- filter_polymorphic_variants(
    genotype_binary[
      variant_types == "MNP",
      ,
      drop = FALSE
    ]
  )

  return(
    list(
      all = all_matrix,
      SNP = snp_matrix,
      INS = insertion_matrix,
      DEL = deletion_matrix,
      MNP = mnp_matrix
    )
  )
}


#' Combine genotype matrices by rows.
#'
#' @param matrix_list List of genotype matrices.
#' @param sample_columns Sample columns to retain and order.
#'
#' @return Combined genotype matrix, or NULL.
combine_matrices <- function(
    matrix_list,
    sample_columns
) {
  if (
    length(matrix_list) == 0
    || length(sample_columns) == 0
  ) {
    return(NULL)
  }

  do.call(
    rbind,
    lapply(
      matrix_list,
      function(genotype_matrix) {
        genotype_matrix[
          ,
          sample_columns,
          drop = FALSE
        ]
      }
    )
  )
}


#' Generate a PCA or MDS ordination plot.
#'
#' @param scores_df Data frame containing coordinates and sample names.
#' @param variance_explained Named vector with variance percentages.
#' @param x_column Name of the x-axis coordinate.
#' @param y_column Name of the y-axis coordinate.
#' @param title Plot title.
#' @param subtitle Optional plot subtitle.
#' @param output_path Output PDF path.
#' @param color_map Named vector of sample colors.
#' @param width Plot width in inches.
#' @param height Plot height in inches.
#'
#' @return Invisibly returns the ggplot object.
plot_ordination <- function(
    scores_df,
    variance_explained = NULL,
    x_column,
    y_column,
    title,
    subtitle = NULL,
    output_path,
    color_map = SAMPLE_COLORS,
    width = 7,
    height = 5.5
) {
  x_label <- if (!is.null(variance_explained)) {
    sprintf(
      "%s (%.1f%%)",
      x_column,
      variance_explained[x_column]
    )
  } else {
    x_column
  }

  y_label <- if (!is.null(variance_explained)) {
    sprintf(
      "%s (%.1f%%)",
      y_column,
      variance_explained[y_column]
    )
  } else {
    y_column
  }

  samples_present <- unique(
    scores_df$sample
  )

  colors_to_use <- color_map[
    samples_present
  ]

  if (any(is.na(colors_to_use))) {
    colors_to_use <- setNames(
      colorRampPalette(
        unname(SAMPLE_COLORS)
      )(
        length(samples_present)
      ),
      samples_present
    )
  }

  ordination_plot <- ggplot(
    scores_df,
    aes(
      x = .data[[x_column]],
      y = .data[[y_column]],
      color = sample,
      label = sample
    )
  ) +
    geom_point(
      size = 5,
      alpha = 0.9
    ) +
    geom_text_repel(
      size = 3.8,
      fontface = "bold",
      box.padding = 0.5,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = colors_to_use
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = x_label,
      y = y_label,
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

  ggsave(
    output_path,
    plot = ordination_plot,
    width = width,
    height = height,
    dpi = 300,
    device = "pdf"
  )

  cat(sprintf(
    "  Saved to: %s\n",
    output_path
  ))

  return(
    invisible(ordination_plot)
  )
}

# =========================
# BUILD GLOBAL MATRICES
# =========================

cat("[1/4] Reading VCF files and building global matrices...\n")

matrices_by_type <- list(
  all = list(),
  SNP = list()
)

sample_sets <- list()
processed_communities <- character()
failed_communities <- character()

for (community in COMMUNITIES) {
  community_dir <- file.path(
    BASE_DIR,
    sprintf(
      "all_pacbio_pansn.fasta.%s.%s",
      GRAPH_HASH,
      community
    )
  )

  vcf_path <- find_vcf(
    community_dir
  )

  if (
    is.na(vcf_path)
    || !file.exists(vcf_path)
  ) {
    cat(sprintf(
      "  [SKIP] %s: VCF not found\n",
      community
    ))

    failed_communities <- c(
      failed_communities,
      community
    )

    next
  }

  matrices <- tryCatch(
    read_vcf_by_type(
      vcf_path
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

  if (is.null(matrices)) {
    failed_communities <- c(
      failed_communities,
      community
    )

    next
  }

  sample_sets[[community]] <- colnames(
    matrices$all
  )

  for (variant_type in c("all", "SNP")) {
    genotype_matrix <- matrices[
      [variant_type]
    ]

    if (
      is.null(genotype_matrix)
      || nrow(genotype_matrix) == 0
    ) {
      next
    }

    rownames(genotype_matrix) <- paste0(
      community,
      "_",
      seq_len(nrow(genotype_matrix))
    )

    matrices_by_type[
      [variant_type]
    ][[community]] <- genotype_matrix
  }

  processed_communities <- c(
    processed_communities,
    community
  )

  cat(sprintf(
    paste0(
      "  [OK] %s: all=%d  SNP=%d  ",
      "INS=%d  DEL=%d  MNP=%d\n"
    ),
    community,
    nrow(matrices$all),
    nrow(matrices$SNP),
    nrow(matrices$INS),
    nrow(matrices$DEL),
    nrow(matrices$MNP)
  ))
}

if (length(processed_communities) == 0) {
  stop(
    "[ERROR] No communities were successfully processed."
  )
}

# Identify samples shared by every successfully processed community
common_samples <- Reduce(
  intersect,
  sample_sets[
    processed_communities
  ]
)

if (length(common_samples) == 0) {
  stop(
    "[ERROR] No samples are shared across the processed communities."
  )
}

cat(sprintf(
  "\n  Samples shared across communities: %s\n",
  paste(common_samples, collapse = ", ")
))

global_all_matrix <- combine_matrices(
  matrices_by_type$all,
  common_samples
)

global_snp_matrix <- combine_matrices(
  matrices_by_type$SNP,
  common_samples
)

if (!is.null(global_all_matrix)) {
  cat(sprintf(
    paste0(
      "\n  Global matrix, all variant types: ",
      "%d variants × %d samples\n"
    ),
    nrow(global_all_matrix),
    ncol(global_all_matrix)
  ))
} else {
  cat(
    "\n  [WARNING] Global all-variant matrix could not be built\n"
  )
}

if (!is.null(global_snp_matrix)) {
  cat(sprintf(
    paste0(
      "  Global matrix, SNPs only: ",
      "%d variants × %d samples\n\n"
    ),
    nrow(global_snp_matrix),
    ncol(global_snp_matrix)
  ))
} else {
  cat(
    "  [WARNING] Global SNP matrix could not be built\n\n"
  )
}

# =========================
# PCA USING SNPs ONLY
# =========================

cat("[2/4] Running PCA using SNPs only...\n")

if (
  !is.null(global_snp_matrix)
  && nrow(global_snp_matrix) >= 2
  && ncol(global_snp_matrix) >= 2
) {
  snp_pca <- prcomp(
    t(global_snp_matrix),
    center = TRUE,
    scale. = FALSE
  )

  snp_variance_explained <- (
    summary(snp_pca)$importance[2, ] * 100
  )

  names(snp_variance_explained) <- paste0(
    "PC",
    seq_along(snp_variance_explained)
  )

  number_of_pca_axes <- min(
    3,
    ncol(snp_pca$x)
  )

  cat(sprintf(
    "  PC1=%.2f%%",
    snp_variance_explained[1]
  ))

  if (number_of_pca_axes >= 2) {
    cat(sprintf(
      "  PC2=%.2f%%",
      snp_variance_explained[2]
    ))
  }

  if (number_of_pca_axes >= 3) {
    cat(sprintf(
      "  PC3=%.2f%%",
      snp_variance_explained[3]
    ))
  }

  cat("\n")

  snp_scores <- as.data.frame(
    snp_pca$x[
      ,
      seq_len(number_of_pca_axes),
      drop = FALSE
    ]
  )

  snp_scores$sample <- rownames(
    snp_scores
  )

  # Plot PC1 versus PC2
  if (number_of_pca_axes >= 2) {
    plot_ordination(
      scores_df = snp_scores,
      variance_explained = snp_variance_explained,
      x_column = "PC1",
      y_column = "PC2",
      title = "PCA — SNPs only across all communities",
      subtitle = sprintf(
        "%d SNP loci · %d samples",
        nrow(global_snp_matrix),
        ncol(global_snp_matrix)
      ),
      output_path = file.path(
        OUTPUT_DIR,
        "snp_pca_PC1vsPC2.pdf"
      )
    )
  }

  # Plot PC1 versus PC3
  if (number_of_pca_axes >= 3) {
    plot_ordination(
      scores_df = snp_scores,
      variance_explained = snp_variance_explained,
      x_column = "PC1",
      y_column = "PC3",
      title = "PCA — SNPs only across all communities",
      output_path = file.path(
        OUTPUT_DIR,
        "snp_pca_PC1vsPC3.pdf"
      )
    )
  }

  # Scree plot
  number_of_scree_axes <- min(
    10,
    length(snp_variance_explained)
  )

  scree_df <- data.frame(
    PC = factor(
      paste0(
        "PC",
        seq_len(number_of_scree_axes)
      ),
      levels = paste0(
        "PC",
        seq_len(number_of_scree_axes)
      )
    ),
    pct_var = snp_variance_explained[
      seq_len(number_of_scree_axes)
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
      title = "Scree plot — SNP PCA",
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

  scree_output <- file.path(
    OUTPUT_DIR,
    "snp_pca_screeplot.pdf"
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
    "  Scree plot saved to: %s\n",
    scree_output
  ))

  # Save PCA coordinates
  snp_coordinates_output <- file.path(
    OUTPUT_DIR,
    "snp_pca_coordinates.tsv"
  )

  write.table(
    snp_scores,
    snp_coordinates_output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  cat(sprintf(
    "  PCA coordinates saved to: %s\n",
    snp_coordinates_output
  ))

  # Save variance explained
  snp_variance_df <- data.frame(
    PC = names(
      snp_variance_explained
    ),
    pct_var = round(
      snp_variance_explained,
      3
    ),
    cum_var = round(
      cumsum(snp_variance_explained),
      3
    )
  )

  snp_variance_output <- file.path(
    OUTPUT_DIR,
    "snp_pca_variance_explained.tsv"
  )

  write.table(
    snp_variance_df[
      seq_len(number_of_scree_axes),
      ,
      drop = FALSE
    ],
    snp_variance_output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  cat(sprintf(
    "  Variance explained saved to: %s\n",
    snp_variance_output
  ))

} else {
  cat(
    "  [SKIP] Insufficient polymorphic SNPs for PCA\n"
  )
}

# =========================
# JACCARD MDS
# =========================

cat(
  "\n[3/4] Running MDS with Jaccard distance across all variant types...\n"
)

if (
  !is.null(global_all_matrix)
  && nrow(global_all_matrix) >= 2
  && ncol(global_all_matrix) >= 2
) {
  number_of_samples <- ncol(
    global_all_matrix
  )

  sample_names <- colnames(
    global_all_matrix
  )

  jaccard_matrix <- matrix(
    0,
    number_of_samples,
    number_of_samples,
    dimnames = list(
      sample_names,
      sample_names
    )
  )

  # Jaccard distance between two binary sample vectors:
  #
  # distance(A, B) = 1 - |A intersection B| / |A union B|
  for (i in seq_len(number_of_samples)) {
    for (j in seq_len(number_of_samples)) {
      if (i == j) {
        next
      }

      sample_a <- global_all_matrix[
        ,
        i
      ]

      sample_b <- global_all_matrix[
        ,
        j
      ]

      union_size <- sum(
        sample_a | sample_b
      )

      intersection_size <- sum(
        sample_a & sample_b
      )

      jaccard_matrix[i, j] <- if (union_size > 0) {
        1 - intersection_size / union_size
      } else {
        0
      }
    }
  }

  cat("  Jaccard distance matrix:\n")
  print(
    round(
      jaccard_matrix,
      4
    )
  )

  # Save the Jaccard distance matrix
  jaccard_output <- file.path(
    OUTPUT_DIR,
    "jaccard_distance_matrix.tsv"
  )

  write.table(
    as.data.frame(jaccard_matrix),
    jaccard_output,
    sep = "\t",
    quote = FALSE,
    row.names = TRUE
  )

  cat(sprintf(
    "  Jaccard matrix saved to: %s\n",
    jaccard_output
  ))

  # Classical MDS, also known as principal coordinates analysis
  maximum_mds_axes <- min(
    3,
    number_of_samples - 1
  )

  mds <- cmdscale(
    as.dist(jaccard_matrix),
    k = maximum_mds_axes,
    eig = TRUE
  )

  mds_scores <- as.data.frame(
    mds$points
  )

  colnames(mds_scores) <- paste0(
    "MDS",
    seq_len(ncol(mds_scores))
  )

  mds_scores$sample <- rownames(
    mds_scores
  )

  # Percentage represented by each MDS axis, calculated from
  # positive eigenvalues
  positive_eigenvalues <- mds$eig[
    mds$eig > 0
  ]

  mds_variance_explained <- (
    positive_eigenvalues
    / sum(positive_eigenvalues)
    * 100
  )

  names(mds_variance_explained) <- paste0(
    "MDS",
    seq_along(mds_variance_explained)
  )

  cat(sprintf(
    "  MDS1=%.2f%%",
    mds_variance_explained[1]
  ))

  if (length(mds_variance_explained) >= 2) {
    cat(sprintf(
      "  MDS2=%.2f%%",
      mds_variance_explained[2]
    ))
  }

  if (length(mds_variance_explained) >= 3) {
    cat(sprintf(
      "  MDS3=%.2f%%",
      mds_variance_explained[3]
    ))
  }

  cat("\n")

  # Plot MDS1 versus MDS2
  if (ncol(mds$points) >= 2) {
    plot_ordination(
      scores_df = mds_scores,
      variance_explained = mds_variance_explained,
      x_column = "MDS1",
      y_column = "MDS2",
      title = "MDS — Jaccard distance across all variant types",
      subtitle = sprintf(
        "%d variant loci · %d samples",
        nrow(global_all_matrix),
        ncol(global_all_matrix)
      ),
      output_path = file.path(
        OUTPUT_DIR,
        "mds_jaccard_MDS1vsMDS2.pdf"
      )
    )
  }

  # Plot MDS1 versus MDS3
  if (ncol(mds$points) >= 3) {
    plot_ordination(
      scores_df = mds_scores,
      variance_explained = mds_variance_explained,
      x_column = "MDS1",
      y_column = "MDS3",
      title = "MDS — Jaccard distance across all variant types",
      output_path = file.path(
        OUTPUT_DIR,
        "mds_jaccard_MDS1vsMDS3.pdf"
      )
    )
  }

  # -----------------------------------------------------------
  # PCA versus MDS comparison panel
  # -----------------------------------------------------------

  if (
    exists("snp_scores")
    && "PC1" %in% colnames(snp_scores)
    && "PC2" %in% colnames(snp_scores)
    && "MDS1" %in% colnames(mds_scores)
    && "MDS2" %in% colnames(mds_scores)
  ) {
    pca_panel_scores <- data.frame(
      Dim1 = snp_scores$PC1,
      Dim2 = snp_scores$PC2,
      sample = snp_scores$sample,
      method = "PCA (SNPs)"
    )

    mds_panel_scores <- data.frame(
      Dim1 = mds_scores$MDS1,
      Dim2 = mds_scores$MDS2,
      sample = mds_scores$sample,
      method = "MDS (Jaccard, all variants)"
    )

    combined_scores <- rbind(
      pca_panel_scores,
      mds_panel_scores
    )

    comparison_panel <- ggplot(
      combined_scores,
      aes(
        x = Dim1,
        y = Dim2,
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
        values = SAMPLE_COLORS,
        na.value = "grey50"
      ) +
      facet_wrap(
        ~ method,
        scales = "free"
      ) +
      labs(
        title = "PCA using SNPs versus MDS using Jaccard distance",
        x = "Dimension 1",
        y = "Dimension 2",
        color = "Sample"
      ) +
      theme_classic(
        base_size = 11
      ) +
      theme(
        plot.title = element_text(
          face = "bold",
          size = 12
        ),
        strip.background = element_rect(
          fill = "grey95",
          color = "grey70"
        ),
        strip.text = element_text(
          face = "bold",
          size = 10
        ),
        panel.grid.major = element_line(
          color = "grey92",
          linewidth = 0.3
        ),
        legend.position = "bottom"
      )

    comparison_output <- file.path(
      OUTPUT_DIR,
      "pca_snps_vs_mds_panel.pdf"
    )

    ggsave(
      comparison_output,
      plot = comparison_panel,
      width = 11,
      height = 5,
      dpi = 300,
      device = "pdf"
    )

    cat(sprintf(
      "  Comparison panel saved to: %s\n",
      comparison_output
    ))
  }

  # Save MDS coordinates
  mds_coordinates_output <- file.path(
    OUTPUT_DIR,
    "mds_jaccard_coordinates.tsv"
  )

  write.table(
    mds_scores,
    mds_coordinates_output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  cat(sprintf(
    "  MDS coordinates saved to: %s\n",
    mds_coordinates_output
  ))

} else {
  cat(
    "  [SKIP] Global all-variant matrix is empty or insufficient\n"
  )
}

# =========================
# FINAL REPORT
# =========================

end_time <- Sys.time()
elapsed_time <- end_time - start_time

cat("\n[4/4] Final report\n")
cat("=============================================================\n")
cat("SNP PCA AND JACCARD MDS SUMMARY\n")
cat(sprintf(
  "Communities requested : %d\n",
  length(COMMUNITIES)
))
cat(sprintf(
  "Communities processed : %d\n",
  length(processed_communities)
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

if (!is.null(global_all_matrix)) {
  cat(sprintf(
    "All-variant loci      : %d\n",
    nrow(global_all_matrix)
  ))
}

if (!is.null(global_snp_matrix)) {
  cat(sprintf(
    "SNP loci              : %d\n",
    nrow(global_snp_matrix)
  ))
}

cat(sprintf(
  "Samples               : %s\n",
  paste(common_samples, collapse = ", ")
))
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
  format(end_time, "%Y-%m-%d %H:%M:%S")
))

cat("\nGenerated files:\n")

generated_files <- list.files(
  OUTPUT_DIR,
  pattern = "snp_pca|mds_jaccard|pca_snps_vs_mds",
  full.names = TRUE
)

if (length(generated_files) > 0) {
  for (generated_file in generated_files) {
    cat(sprintf(
      "  %s\n",
      basename(generated_file)
    ))
  }
} else {
  cat("  No matching output files were generated.\n")
}

cat("=============================================================\n")
