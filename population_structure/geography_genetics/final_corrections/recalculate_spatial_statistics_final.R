#!/usr/bin/env Rscript

# =============================================================================
# recalculate_spatial_statistics_final.R
#
# Final corrected spatial-genetic statistics for the PGGB24 analysis.
#
# Corrections implemented:
#   1. Ordered geographic-distance classes.
#   2. Formal model residuals instead of z-score subtraction.
#   3. Standardized MRM output with coefficients, p-values, R2 and model p-value.
#   4. Explicit treatment of zero-km pairs.
#   5. Sensitivity analysis for co-located municipal coordinates.
#   6. High-precision coordinate subset.
#
# Usage:
#   Rscript recalculate_spatial_statistics_final.R \
#       validated_coordinates.tsv metadata.tsv \
#       all11_distance_matrix.tsv conservative9_distance_matrix.tsv \
#       output_directory
# =============================================================================

required_packages <- c("geosphere", "vegan", "ecodist")

for (package_name in required_packages) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    stop(sprintf("Required R package not installed: %s", package_name))
  }
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
  stop(
    paste(
      "Usage: Rscript recalculate_spatial_statistics_final.R",
      "COORDINATES METADATA ALL11_DISTANCE CONS9_DISTANCE OUTPUT_DIR"
    )
  )
}

COORDINATE_FILE <- args[1]
METADATA_FILE <- args[2]
ALL11_FILE <- args[3]
CONSERVATIVE9_FILE <- args[4]
OUTPUT_DIR <- args[5]

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

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

read_distance_matrix <- function(path) {
  matrix_df <- read.delim(
    path,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  as.matrix(matrix_df)
}

all11 <- read_distance_matrix(ALL11_FILE)
conservative9 <- read_distance_matrix(CONSERVATIVE9_FILE)

valid <- coordinates[
  !is.na(coordinates$latitude) &
  !is.na(coordinates$longitude) &
  coordinates$latitude != "" &
  coordinates$longitude != "",
  ,
  drop = FALSE
]

valid$latitude <- as.numeric(valid$latitude)
valid$longitude <- as.numeric(valid$longitude)

samples <- valid$sample_id

if (
  !all(samples %in% rownames(all11)) ||
  !all(samples %in% rownames(conservative9))
) {
  stop("At least one spatial sample is absent from a genetic-distance matrix.")
}

all11 <- all11[samples, samples, drop = FALSE]
conservative9 <- conservative9[samples, samples, drop = FALSE]

coordinate_matrix <- as.matrix(valid[, c("longitude", "latitude")])

geographic_distance <- geosphere::distm(
  coordinate_matrix,
  fun = geosphere::distGeo
) / 1000

rownames(geographic_distance) <- samples
colnames(geographic_distance) <- samples

write.table(
  geographic_distance,
  file.path(OUTPUT_DIR, "geographic_distance_matrix_km_final.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

upper_triangle_table <- function(matrix_object) {
  positions <- which(upper.tri(matrix_object), arr.ind = TRUE)

  data.frame(
    sample_1 = rownames(matrix_object)[positions[, 1]],
    sample_2 = colnames(matrix_object)[positions[, 2]],
    value = matrix_object[positions],
    stringsAsFactors = FALSE
  )
}

geographic_long <- upper_triangle_table(geographic_distance)
all11_long <- upper_triangle_table(all11)
conservative9_long <- upper_triangle_table(conservative9)

pairs <- geographic_long
colnames(pairs)[3] <- "geographic_distance_km"
pairs$genetic_distance_all11 <- all11_long$value
pairs$genetic_distance_conservative9 <- conservative9_long$value

metadata_1 <- metadata[
  match(pairs$sample_1, metadata$sample_id),
  ,
  drop = FALSE
]

metadata_2 <- metadata[
  match(pairs$sample_2, metadata$sample_id),
  ,
  drop = FALSE
]

coordinate_1 <- valid[
  match(pairs$sample_1, valid$sample_id),
  ,
  drop = FALSE
]

coordinate_2 <- valid[
  match(pairs$sample_2, valid$sample_id),
  ,
  drop = FALSE
]

pairs$country_1 <- metadata_1$country
pairs$country_2 <- metadata_2$country
pairs$region_1 <- metadata_1$region_state
pairs$region_2 <- metadata_2$region_state
pairs$technology_1 <- metadata_1$technology
pairs$technology_2 <- metadata_2$technology
pairs$coordinate_precision_1 <- coordinate_1$coordinate_precision
pairs$coordinate_precision_2 <- coordinate_2$coordinate_precision

pairs$same_country <- ifelse(
  pairs$country_1 == pairs$country_2,
  "same country",
  "different countries"
)

pairs$same_region <- ifelse(
  pairs$country_1 == pairs$country_2 &
  pairs$region_1 == pairs$region_2,
  "same region",
  "different regions"
)

pairs$technology_pair <- ifelse(
  pairs$technology_1 == pairs$technology_2,
  paste0(pairs$technology_1, "-", pairs$technology_2),
  "mixed"
)

pairs$coordinate_precision_pair <- ifelse(
  pairs$coordinate_precision_1 == "state" |
  pairs$coordinate_precision_2 == "state",
  "includes state-level coordinate",
  "locality/municipality only"
)

pairs$zero_km_pair <- pairs$geographic_distance_km < 1e-9
pairs$zero_km_interpretation <- ifelse(
  pairs$zero_km_pair,
  "co-located at available municipal/locality resolution",
  "spatially separated at available resolution"
)

# -----------------------------------------------------------------------------
# Ordered geographic-distance classes
# -----------------------------------------------------------------------------

maximum_distance <- ceiling(max(pairs$geographic_distance_km, na.rm = TRUE))

distance_breaks <- c(
  0,
  246,
  780,
  1099,
  1735,
  maximum_distance
)

distance_labels <- c(
  "[0, 246]",
  "(246, 780]",
  "(780, 1099]",
  "(1099, 1735]",
  sprintf("(1735, %d]", maximum_distance)
)

pairs$distance_class <- cut(
  pairs$geographic_distance_km,
  breaks = distance_breaks,
  labels = distance_labels,
  include.lowest = TRUE,
  right = TRUE,
  ordered_result = TRUE
)

# -----------------------------------------------------------------------------
# Formal residuals from fitted distance-decay models
#
# These are used only to rank mismatched pairs. Ordinary LM p-values are not
# used for inference because pairwise rows are non-independent.
# -----------------------------------------------------------------------------

fit_all11 <- lm(
  genetic_distance_all11 ~ log10(geographic_distance_km + 1),
  data = pairs
)

fit_conservative9 <- lm(
  genetic_distance_conservative9 ~ log10(geographic_distance_km + 1),
  data = pairs
)

pairs$model_fitted_all11 <- fitted(fit_all11)
pairs$model_fitted_conservative9 <- fitted(fit_conservative9)

pairs$model_residual_all11 <- residuals(fit_all11)
pairs$model_residual_conservative9 <- residuals(fit_conservative9)

pairs$standardized_residual_all11 <- as.numeric(
  scale(pairs$model_residual_all11)
)

pairs$standardized_residual_conservative9 <- as.numeric(
  scale(pairs$model_residual_conservative9)
)

pairs$residual_interpretation_all11 <- ifelse(
  pairs$standardized_residual_all11 < 0,
  "genetically closer than fitted distance-decay model",
  "genetically farther than fitted distance-decay model"
)

pairs$residual_interpretation_conservative9 <- ifelse(
  pairs$standardized_residual_conservative9 < 0,
  "genetically closer than fitted distance-decay model",
  "genetically farther than fitted distance-decay model"
)

write.table(
  pairs,
  file.path(OUTPUT_DIR, "geographic_genetic_pair_table_final.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Descriptive regressions
# -----------------------------------------------------------------------------

extract_lm_summary <- function(fit, profile, subset_name, n_pairs) {
  fit_summary <- summary(fit)

  data.frame(
    profile = profile,
    subset = subset_name,
    n_pairs = n_pairs,
    intercept = unname(coef(fit)[1]),
    slope_log10_km = unname(coef(fit)[2]),
    r_squared = fit_summary$r.squared,
    adjusted_r_squared = fit_summary$adj.r.squared,
    note = paste(
      "Effect-size description only;",
      "pairwise rows are non-independent and LM p-values are not interpreted."
    )
  )
}

regression_results <- rbind(
  extract_lm_summary(
    fit_all11,
    "all11",
    "all_valid_pairs",
    nrow(pairs)
  ),
  extract_lm_summary(
    fit_conservative9,
    "conservative9",
    "all_valid_pairs",
    nrow(pairs)
  )
)

nonzero_pairs <- pairs[!pairs$zero_km_pair, , drop = FALSE]

fit_all11_nonzero <- lm(
  genetic_distance_all11 ~ log10(geographic_distance_km + 1),
  data = nonzero_pairs
)

fit_conservative9_nonzero <- lm(
  genetic_distance_conservative9 ~ log10(geographic_distance_km + 1),
  data = nonzero_pairs
)

regression_results <- rbind(
  regression_results,
  extract_lm_summary(
    fit_all11_nonzero,
    "all11",
    "excluding_zero_km_pairs",
    nrow(nonzero_pairs)
  ),
  extract_lm_summary(
    fit_conservative9_nonzero,
    "conservative9",
    "excluding_zero_km_pairs",
    nrow(nonzero_pairs)
  )
)

write.table(
  regression_results,
  file.path(OUTPUT_DIR, "descriptive_regression_results_final.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Mantel tests
# -----------------------------------------------------------------------------

subset_definitions <- list(
  all_valid_coordinates = samples,
  high_precision = valid$sample_id[
    valid$coordinate_precision %in% c("locality", "municipality")
  ],
  Colombia_only = valid$sample_id[valid$country == "Colombia"],
  Illumina_only = valid$sample_id[
    metadata$technology[
      match(valid$sample_id, metadata$sample_id)
    ] == "Illumina"
  ],
  excluding_Colombia = valid$sample_id[valid$country != "Colombia"]
)

mantel_rows <- list()
subset_rows <- list()

for (subset_name in names(subset_definitions)) {
  subset_samples <- unique(subset_definitions[[subset_name]])
  subset_samples <- subset_samples[subset_samples %in% samples]

  subset_rows[[length(subset_rows) + 1]] <- data.frame(
    subset = subset_name,
    n_samples = length(subset_samples),
    samples = paste(subset_samples, collapse = ",")
  )

  if (length(subset_samples) < 4) {
    next
  }

  geographic_subset <- geographic_distance[
    subset_samples,
    subset_samples,
    drop = FALSE
  ]

  for (profile in c("all11", "conservative9")) {
    genetic_full <- if (profile == "all11") {
      all11
    } else {
      conservative9
    }

    genetic_subset <- genetic_full[
      subset_samples,
      subset_samples,
      drop = FALSE
    ]

    for (method_name in c("pearson", "spearman")) {
      mantel_result <- vegan::mantel(
        as.dist(genetic_subset),
        as.dist(geographic_subset),
        method = method_name,
        permutations = 9999,
        na.rm = TRUE
      )

      mantel_rows[[length(mantel_rows) + 1]] <- data.frame(
        profile = profile,
        subset = subset_name,
        method = method_name,
        n_samples = length(subset_samples),
        statistic = unname(mantel_result$statistic),
        p_value = mantel_result$signif,
        permutations = mantel_result$permutations
      )
    }
  }
}

mantel_results <- do.call(rbind, mantel_rows)
subset_summary <- do.call(rbind, subset_rows)

write.table(
  subset_summary,
  file.path(OUTPUT_DIR, "spatial_subset_summary_final.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  mantel_results,
  file.path(OUTPUT_DIR, "mantel_results_final.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Co-location sensitivity
#
# For every duplicated coordinate group, retain one representative at a time.
# With three duplicated pairs, this creates 2^3 = 8 representative datasets.
# This avoids selectively deleting individual matrix cells.
# -----------------------------------------------------------------------------

coordinate_key <- paste(
  sprintf("%.8f", valid$latitude),
  sprintf("%.8f", valid$longitude),
  sep = "|"
)

coordinate_groups <- split(valid$sample_id, coordinate_key)
duplicated_groups <- coordinate_groups[lengths(coordinate_groups) > 1]

expand_representatives <- function(groups) {
  if (length(groups) == 0) {
    return(data.frame(combination_id = 1))
  }

  combinations <- expand.grid(
    groups,
    stringsAsFactors = FALSE,
    KEEP.OUT.ATTRS = FALSE
  )

  colnames(combinations) <- paste0(
    "duplicate_group_",
    seq_len(ncol(combinations))
  )

  combinations$combination_id <- seq_len(nrow(combinations))
  combinations
}

representative_combinations <- expand_representatives(duplicated_groups)

colocation_mantel_rows <- list()

for (row_index in seq_len(nrow(representative_combinations))) {
  selected_representatives <- if (length(duplicated_groups) > 0) {
    as.character(
      representative_combinations[
        row_index,
        grep("^duplicate_group_", colnames(representative_combinations)),
        drop = TRUE
      ]
    )
  } else {
    character(0)
  }

  duplicated_samples <- unlist(duplicated_groups, use.names = FALSE)

  retained_samples <- c(
    setdiff(samples, duplicated_samples),
    selected_representatives
  )

  retained_samples <- unique(retained_samples)

  if (length(retained_samples) < 4) {
    next
  }

  for (profile in c("all11", "conservative9")) {
    genetic_full <- if (profile == "all11") {
      all11
    } else {
      conservative9
    }

    genetic_subset <- genetic_full[
      retained_samples,
      retained_samples,
      drop = FALSE
    ]

    geographic_subset <- geographic_distance[
      retained_samples,
      retained_samples,
      drop = FALSE
    ]

    for (method_name in c("pearson", "spearman")) {
      mantel_result <- vegan::mantel(
        as.dist(genetic_subset),
        as.dist(geographic_subset),
        method = method_name,
        permutations = 9999,
        na.rm = TRUE
      )

      colocation_mantel_rows[[length(colocation_mantel_rows) + 1]] <- data.frame(
        combination_id = representative_combinations$combination_id[row_index],
        profile = profile,
        method = method_name,
        n_samples = length(retained_samples),
        selected_representatives = paste(
          selected_representatives,
          collapse = ","
        ),
        statistic = unname(mantel_result$statistic),
        p_value = mantel_result$signif,
        permutations = mantel_result$permutations
      )
    }
  }
}

colocation_mantel <- do.call(rbind, colocation_mantel_rows)

write.table(
  colocation_mantel,
  file.path(OUTPUT_DIR, "colocation_representative_mantel_results.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

summarize_colocation <- function(dataframe) {
  result_rows <- list()

  grouping_keys <- unique(
    dataframe[, c("profile", "method"), drop = FALSE]
  )

  for (index in seq_len(nrow(grouping_keys))) {
    profile_name <- grouping_keys$profile[index]
    method_name <- grouping_keys$method[index]

    subset <- dataframe[
      dataframe$profile == profile_name &
      dataframe$method == method_name,
      ,
      drop = FALSE
    ]

    result_rows[[length(result_rows) + 1]] <- data.frame(
      profile = profile_name,
      method = method_name,
      combinations = nrow(subset),
      statistic_median = median(subset$statistic),
      statistic_q025 = unname(quantile(subset$statistic, 0.025)),
      statistic_q975 = unname(quantile(subset$statistic, 0.975)),
      p_value_median = median(subset$p_value),
      significant_fraction_0.05 = mean(subset$p_value <= 0.05)
    )
  }

  do.call(rbind, result_rows)
}

colocation_summary <- summarize_colocation(colocation_mantel)

write.table(
  colocation_summary,
  file.path(OUTPUT_DIR, "colocation_representative_mantel_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# MRM models
#
# MRM uses permutation tests for coefficients and model R2.
# Models are run for all valid coordinates and the high-precision subset.
# -----------------------------------------------------------------------------

technology_vector <- metadata$technology[
  match(samples, metadata$sample_id)
]

technology_distance <- outer(
  technology_vector,
  technology_vector,
  FUN = function(value_1, value_2) {
    as.numeric(value_1 != value_2)
  }
)

diag(technology_distance) <- 0
rownames(technology_distance) <- samples
colnames(technology_distance) <- samples

normalize_mrm_coefficients <- function(
  fit,
  profile,
  subset_name,
  n_samples
) {
  coefficient_matrix <- as.data.frame(
    fit$coef,
    check.names = FALSE
  )

  if (ncol(coefficient_matrix) < 2) {
    stop("Unexpected ecodist::MRM coefficient structure.")
  }

  data.frame(
    profile = profile,
    subset = subset_name,
    n_samples = n_samples,
    term = rownames(coefficient_matrix),
    estimate = as.numeric(coefficient_matrix[[1]]),
    p_value = as.numeric(coefficient_matrix[[2]]),
    permutations = 9999,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

normalize_mrm_model <- function(
  fit,
  profile,
  subset_name,
  n_samples
) {
  r_squared_values <- as.numeric(fit$r.squared)
  f_test_values <- as.numeric(fit$F.test)

  data.frame(
    profile = profile,
    subset = subset_name,
    n_samples = n_samples,
    r_squared = ifelse(
      length(r_squared_values) >= 1,
      r_squared_values[1],
      NA_real_
    ),
    r_squared_p_value = ifelse(
      length(r_squared_values) >= 2,
      r_squared_values[2],
      NA_real_
    ),
    f_statistic = ifelse(
      length(f_test_values) >= 1,
      f_test_values[1],
      NA_real_
    ),
    f_test_p_value = ifelse(
      length(f_test_values) >= 2,
      f_test_values[2],
      NA_real_
    ),
    permutations = 9999
  )
}

mrm_coefficient_rows <- list()
mrm_model_rows <- list()

mrm_subsets <- list(
  all_valid_coordinates = samples,
  high_precision = valid$sample_id[
    valid$coordinate_precision %in% c("locality", "municipality")
  ]
)

for (subset_name in names(mrm_subsets)) {
  subset_samples <- mrm_subsets[[subset_name]]
  subset_samples <- subset_samples[subset_samples %in% samples]

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

  for (profile in c("all11", "conservative9")) {
    genetic_full <- if (profile == "all11") {
      all11
    } else {
      conservative9
    }

    genetic_subset <- genetic_full[
      subset_samples,
      subset_samples,
      drop = FALSE
    ]

    fit <- ecodist::MRM(
      as.dist(genetic_subset) ~
        as.dist(log10(geographic_subset + 1)) +
        as.dist(technology_subset),
      nperm = 9999
    )

    mrm_coefficient_rows[[length(mrm_coefficient_rows) + 1]] <-
      normalize_mrm_coefficients(
        fit,
        profile,
        subset_name,
        length(subset_samples)
      )

    mrm_model_rows[[length(mrm_model_rows) + 1]] <-
      normalize_mrm_model(
        fit,
        profile,
        subset_name,
        length(subset_samples)
      )
  }
}

mrm_coefficients <- do.call(rbind, mrm_coefficient_rows)
mrm_models <- do.call(rbind, mrm_model_rows)

write.table(
  mrm_coefficients,
  file.path(OUTPUT_DIR, "mrm_coefficients_final.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  mrm_models,
  file.path(OUTPUT_DIR, "mrm_model_summary_final.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Distance-class summary and residual-ranked pairs
# -----------------------------------------------------------------------------

summarize_distance_classes <- function(dataframe, value_column, profile) {
  rows <- list()

  for (distance_level in levels(dataframe$distance_class)) {
    values <- dataframe[
      dataframe$distance_class == distance_level,
      value_column
    ]

    values <- values[!is.na(values)]

    rows[[length(rows) + 1]] <- data.frame(
      profile = profile,
      distance_class = distance_level,
      n_pairs = length(values),
      median = median(values),
      mean = mean(values),
      sd = sd(values),
      q025 = unname(quantile(values, 0.025)),
      q975 = unname(quantile(values, 0.975))
    )
  }

  do.call(rbind, rows)
}

distance_class_summary <- rbind(
  summarize_distance_classes(
    pairs,
    "genetic_distance_all11",
    "all11"
  ),
  summarize_distance_classes(
    pairs,
    "genetic_distance_conservative9",
    "conservative9"
  )
)

write.table(
  distance_class_summary,
  file.path(OUTPUT_DIR, "distance_class_summary_final.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

all11_residual_ranking <- pairs[
  order(pairs$standardized_residual_all11),
  ,
  drop = FALSE
]

conservative9_residual_ranking <- pairs[
  order(pairs$standardized_residual_conservative9),
  ,
  drop = FALSE
]

write.table(
  all11_residual_ranking,
  file.path(OUTPUT_DIR, "model_residual_pairs_all11.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  conservative9_residual_ranking,
  file.path(OUTPUT_DIR, "model_residual_pairs_conservative9.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Final corrected spatial-genetic statistics completed.\n")
