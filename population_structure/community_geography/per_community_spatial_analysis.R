#!/usr/bin/env Rscript

# =============================================================================
# per_community_spatial_analysis.R
#
# Computes population-structure and geography-versus-genetics statistics
# independently for every PGGB community represented in a final SNP matrix.
#
# For each community:
#   - retained SNP count and sample count
#   - normalized Manhattan genetic-distance matrix
#   - PCA coordinates and variance explained
#   - Mantel Pearson and Spearman tests
#   - MRM: genetic distance ~ log10(geographic distance + 1)
#          + sequencing-technology distance
#   - descriptive distance-decay slope and R2
#
# The complete and high-precision coordinate subsets are evaluated.
#
# Usage:
# Rscript per_community_spatial_analysis.R \
#   SNP_MATRIX FEATURE_METADATA GEO_MATRIX COORDINATES METADATA \
#   PROFILE OUTPUT_DIR MIN_SNPS PERMUTATIONS
# =============================================================================

required_packages <- c("vegan", "ecodist")

for (package_name in required_packages) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    stop(sprintf("Required R package not installed: %s", package_name))
  }
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 9) {
  stop(
    paste(
      "Usage: Rscript per_community_spatial_analysis.R",
      "MATRIX FEATURES GEO COORDINATES METADATA PROFILE OUTPUT",
      "MIN_SNPS PERMUTATIONS"
    )
  )
}

MATRIX_FILE <- args[1]
FEATURE_FILE <- args[2]
GEO_FILE <- args[3]
COORDINATE_FILE <- args[4]
METADATA_FILE <- args[5]
PROFILE <- args[6]
OUTPUT_DIR <- args[7]
MIN_SNPS <- as.integer(args[8])
PERMUTATIONS <- as.integer(args[9])

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "distance_matrices"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "pca"), recursive = TRUE, showWarnings = FALSE)

matrix_df <- read.delim(
  MATRIX_FILE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

feature_df <- read.delim(
  FEATURE_FILE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

geographic_distance <- as.matrix(
  read.delim(
    GEO_FILE,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
)

coordinates <- read.delim(
  COORDINATE_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

metadata <- read.delim(
  METADATA_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!"community" %in% colnames(feature_df)) {
  stop("Feature metadata must contain a column named 'community'.")
}

common_features <- intersect(rownames(matrix_df), rownames(feature_df))
matrix_df <- matrix_df[common_features, , drop = FALSE]
feature_df <- feature_df[common_features, , drop = FALSE]

spatial_samples <- intersect(colnames(matrix_df), rownames(geographic_distance))
matrix_df <- matrix_df[, spatial_samples, drop = FALSE]
geographic_distance <- geographic_distance[spatial_samples, spatial_samples, drop = FALSE]

technology <- metadata$technology[match(spatial_samples, metadata$sample_id)]
technology_distance <- outer(
  technology,
  technology,
  FUN = function(value_1, value_2) as.numeric(value_1 != value_2)
)
diag(technology_distance) <- 0
rownames(technology_distance) <- spatial_samples
colnames(technology_distance) <- spatial_samples

high_precision_samples <- coordinates$sample_id[
  coordinates$coordinate_precision %in% c("locality", "municipality")
]
high_precision_samples <- intersect(high_precision_samples, spatial_samples)

subset_definitions <- list(
  all_valid_coordinates = spatial_samples,
  high_precision = high_precision_samples
)

normalize_mrm_coefficients <- function(fit) {
  coefficient_matrix <- as.data.frame(fit$coef, check.names = FALSE)
  if (ncol(coefficient_matrix) < 2) {
    stop("Unexpected ecodist::MRM coefficient structure.")
  }

  data.frame(
    term = rownames(coefficient_matrix),
    estimate = as.numeric(coefficient_matrix[[1]]),
    p_value = as.numeric(coefficient_matrix[[2]]),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

normalize_mrm_model <- function(fit) {
  r_squared_values <- as.numeric(fit$r.squared)
  f_test_values <- as.numeric(fit$F.test)

  data.frame(
    r_squared = ifelse(length(r_squared_values) >= 1, r_squared_values[1], NA_real_),
    r_squared_p_value = ifelse(length(r_squared_values) >= 2, r_squared_values[2], NA_real_),
    f_statistic = ifelse(length(f_test_values) >= 1, f_test_values[1], NA_real_),
    f_test_p_value = ifelse(length(f_test_values) >= 2, f_test_values[2], NA_real_)
  )
}

communities <- sort(unique(feature_df$community))
summary_rows <- list()
mantel_rows <- list()
mrm_coefficient_rows <- list()
mrm_model_rows <- list()
regression_rows <- list()
pca_variance_rows <- list()

for (community_name in communities) {
  community_features <- rownames(feature_df)[feature_df$community == community_name]
  community_features <- intersect(community_features, rownames(matrix_df))

  retained_snps <- length(community_features)

  if (retained_snps < MIN_SNPS) {
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      profile = PROFILE,
      community = community_name,
      retained_snps = retained_snps,
      variable_snps = 0,
      samples = ncol(matrix_df),
      analysis_status = sprintf("skipped: fewer than %d retained SNPs", MIN_SNPS)
    )
    next
  }

  community_matrix <- t(as.matrix(matrix_df[community_features, , drop = FALSE]))
  variable_columns <- apply(community_matrix, 2, function(column) var(column) > 0)
  community_matrix <- community_matrix[, variable_columns, drop = FALSE]
  variable_snps <- ncol(community_matrix)

  if (variable_snps < MIN_SNPS) {
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      profile = PROFILE,
      community = community_name,
      retained_snps = retained_snps,
      variable_snps = variable_snps,
      samples = nrow(community_matrix),
      analysis_status = sprintf("skipped: fewer than %d variable SNPs", MIN_SNPS)
    )
    next
  }

  genetic_distance <- as.matrix(
    dist(community_matrix, method = "manhattan")
  ) / variable_snps

  write.table(
    genetic_distance,
    file.path(
      OUTPUT_DIR,
      "distance_matrices",
      paste0(PROFILE, "_", community_name, "_genetic_distance.tsv")
    ),
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )

  pca <- prcomp(
    community_matrix,
    center = TRUE,
    scale. = FALSE
  )

  variance_explained <- summary(pca)$importance[2, ] * 100
  number_pcs <- min(5, ncol(pca$x))

  pca_scores <- as.data.frame(pca$x[, seq_len(number_pcs), drop = FALSE])
  pca_scores$sample_id <- rownames(pca_scores)
  pca_scores$profile <- PROFILE
  pca_scores$community <- community_name

  pca_scores <- merge(
    pca_scores,
    metadata,
    by = "sample_id",
    all.x = TRUE,
    sort = FALSE
  )

  write.table(
    pca_scores,
    file.path(
      OUTPUT_DIR,
      "pca",
      paste0(PROFILE, "_", community_name, "_pca_coordinates.tsv")
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  for (pc_index in seq_len(number_pcs)) {
    pca_variance_rows[[length(pca_variance_rows) + 1]] <- data.frame(
      profile = PROFILE,
      community = community_name,
      PC = paste0("PC", pc_index),
      pct_variance = variance_explained[pc_index],
      cumulative_variance = sum(variance_explained[seq_len(pc_index)])
    )
  }

  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    profile = PROFILE,
    community = community_name,
    retained_snps = retained_snps,
    variable_snps = variable_snps,
    samples = nrow(community_matrix),
    analysis_status = "completed"
  )

  for (subset_name in names(subset_definitions)) {
    subset_samples <- subset_definitions[[subset_name]]
    subset_samples <- intersect(subset_samples, rownames(genetic_distance))

    if (length(subset_samples) < 4) {
      next
    }

    genetic_subset <- genetic_distance[
      subset_samples,
      subset_samples,
      drop = FALSE
    ]

    geographic_subset <- geographic_distance[
      subset_samples,
      subset_samples,
      drop = FALSE
    ]

    technology_subset <- technology_distance[
      subset_samples,
      subset_samples,
      drop = FALSE
    ]

    for (method_name in c("pearson", "spearman")) {
      mantel_result <- vegan::mantel(
        as.dist(genetic_subset),
        as.dist(geographic_subset),
        method = method_name,
        permutations = PERMUTATIONS,
        na.rm = TRUE
      )

      mantel_rows[[length(mantel_rows) + 1]] <- data.frame(
        profile = PROFILE,
        community = community_name,
        subset = subset_name,
        method = method_name,
        n_samples = length(subset_samples),
        retained_snps = retained_snps,
        variable_snps = variable_snps,
        statistic = unname(mantel_result$statistic),
        p_value = mantel_result$signif,
        permutations = mantel_result$permutations
      )
    }

    mrm_fit <- ecodist::MRM(
      as.dist(genetic_subset) ~
        as.dist(log10(geographic_subset + 1)) +
        as.dist(technology_subset),
      nperm = PERMUTATIONS
    )

    mrm_coefficients <- normalize_mrm_coefficients(mrm_fit)
    mrm_coefficients$profile <- PROFILE
    mrm_coefficients$community <- community_name
    mrm_coefficients$subset <- subset_name
    mrm_coefficients$n_samples <- length(subset_samples)
    mrm_coefficients$retained_snps <- retained_snps
    mrm_coefficients$variable_snps <- variable_snps
    mrm_coefficients$permutations <- PERMUTATIONS

    mrm_coefficient_rows[[length(mrm_coefficient_rows) + 1]] <- mrm_coefficients

    mrm_model <- normalize_mrm_model(mrm_fit)
    mrm_model$profile <- PROFILE
    mrm_model$community <- community_name
    mrm_model$subset <- subset_name
    mrm_model$n_samples <- length(subset_samples)
    mrm_model$retained_snps <- retained_snps
    mrm_model$variable_snps <- variable_snps
    mrm_model$permutations <- PERMUTATIONS

    mrm_model_rows[[length(mrm_model_rows) + 1]] <- mrm_model

    pair_positions <- upper.tri(genetic_subset)
    genetic_vector <- genetic_subset[pair_positions]
    geographic_vector <- geographic_subset[pair_positions]

    regression_fit <- lm(
      genetic_vector ~ log10(geographic_vector + 1)
    )

    regression_summary <- summary(regression_fit)

    regression_rows[[length(regression_rows) + 1]] <- data.frame(
      profile = PROFILE,
      community = community_name,
      subset = subset_name,
      n_samples = length(subset_samples),
      n_pairs = length(genetic_vector),
      retained_snps = retained_snps,
      variable_snps = variable_snps,
      intercept = unname(coef(regression_fit)[1]),
      slope_log10_km = unname(coef(regression_fit)[2]),
      r_squared = regression_summary$r.squared,
      adjusted_r_squared = regression_summary$adj.r.squared,
      note = "Descriptive effect size only; pairwise rows are non-independent."
    )
  }
}

write.table(
  do.call(rbind, summary_rows),
  file.path(OUTPUT_DIR, paste0(PROFILE, "_community_analysis_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  do.call(rbind, mantel_rows),
  file.path(OUTPUT_DIR, paste0(PROFILE, "_community_mantel.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  do.call(rbind, mrm_coefficient_rows),
  file.path(OUTPUT_DIR, paste0(PROFILE, "_community_mrm_coefficients.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  do.call(rbind, mrm_model_rows),
  file.path(OUTPUT_DIR, paste0(PROFILE, "_community_mrm_models.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  do.call(rbind, regression_rows),
  file.path(OUTPUT_DIR, paste0(PROFILE, "_community_regression.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  do.call(rbind, pca_variance_rows),
  file.path(OUTPUT_DIR, paste0(PROFILE, "_community_pca_variance.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(sprintf("Per-community spatial analysis completed for profile %s.\n", PROFILE))
