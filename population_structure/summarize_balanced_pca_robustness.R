#!/usr/bin/env Rscript

# =============================================================
# summarize_balanced_pca_robustness.R
#
# Summarizes balanced-PCA robustness replicates for PGGB and Cactus.
#
# Each replicate is aligned to the complete-SNP PCA using orthogonal
# Procrustes transformation. The script then calculates:
#
#   - Procrustes residual sum of squares
#   - Pearson and Spearman correlation of pairwise sample distances
#   - Mean and standard deviation of aligned PC coordinates
#   - 95% empirical coordinate intervals
#   - Replicate-cloud plots colored by geographic region
#   - Metric distribution plots
#
# Usage:
#   Rscript summarize_balanced_pca_robustness.R \
#       ROBUSTNESS_DIR \
#       FULL_PGGB_COORDINATES.tsv \
#       FULL_CACTUS_COORDINATES.tsv \
#       OUTPUT_DIR
# =============================================================

required_packages <- c(
  "ggplot2",
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
  library(scales)
})

arguments <- commandArgs(
  trailingOnly = TRUE
)

if (length(arguments) != 4) {
  stop(
    paste(
      "Usage:",
      "Rscript summarize_balanced_pca_robustness.R",
      "ROBUSTNESS_DIR FULL_PGGB_COORDINATES.tsv",
      "FULL_CACTUS_COORDINATES.tsv OUTPUT_DIR"
    )
  )
}

ROBUSTNESS_DIR <- arguments[1]
FULL_PGGB_FILE <- arguments[2]
FULL_CACTUS_FILE <- arguments[3]
OUTPUT_DIR <- arguments[4]

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# =========================
# FUNCTIONS
# =========================

read_coordinates <- function(path) {
  table <- read.delim(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required <- c(
    "sample_id",
    "PC1",
    "PC2"
  )

  missing <- setdiff(
    required,
    colnames(table)
  )

  if (length(missing) > 0) {
    stop(
      sprintf(
        "Coordinate file is missing columns %s: %s",
        paste(missing, collapse = ", "),
        path
      )
    )
  }

  table
}


center_matrix <- function(matrix_data) {
  scale(
    matrix_data,
    center = TRUE,
    scale = FALSE
  )
}


procrustes_align <- function(target, configuration) {
  target_centered <- center_matrix(
    target
  )

  configuration_centered <- center_matrix(
    configuration
  )

  cross_product <- t(
    configuration_centered
  ) %*% target_centered

  decomposition <- svd(
    cross_product
  )

  rotation <- decomposition$u %*% t(
    decomposition$v
  )

  rotated <- configuration_centered %*% rotation

  denominator <- sum(
    rotated^2
  )

  scale_factor <- if (denominator > 0) {
    sum(rotated * target_centered) / denominator
  } else {
    1
  }

  aligned <- rotated * scale_factor

  residual_ss <- sum(
    (target_centered - aligned)^2
  )

  list(
    aligned = aligned,
    residual_ss = residual_ss,
    scale = scale_factor
  )
}


distance_vector <- function(matrix_data) {
  as.vector(
    dist(matrix_data)
  )
}


summarize_method <- function(
    method,
    full_coordinates_file,
    robustness_dir,
    output_dir
) {
  full_table <- read_coordinates(
    full_coordinates_file
  )

  full_table <- full_table[
    order(full_table$sample_id),
    ,
    drop = FALSE
  ]

  target <- as.matrix(
    full_table[, c("PC1", "PC2")]
  )

  rownames(target) <- full_table$sample_id

  coordinate_files <- Sys.glob(
    file.path(
      robustness_dir,
      "seed_*",
      method,
      "pca_coordinates.tsv"
    )
  )

  if (length(coordinate_files) == 0) {
    stop(
      sprintf(
        "No robustness coordinate files found for %s",
        method
      )
    )
  }

  aligned_rows <- list()
  metric_rows <- list()

  target_distances <- distance_vector(
    target
  )

  for (coordinate_file in coordinate_files) {
    replicate <- read_coordinates(
      coordinate_file
    )

    replicate <- replicate[
      order(replicate$sample_id),
      ,
      drop = FALSE
    ]

    if (!identical(
      replicate$sample_id,
      full_table$sample_id
    )) {
      stop(
        sprintf(
          "Sample mismatch in: %s",
          coordinate_file
        )
      )
    }

    configuration <- as.matrix(
      replicate[, c("PC1", "PC2")]
    )

    alignment <- procrustes_align(
      target,
      configuration
    )

    seed <- if ("seed" %in% colnames(replicate)) {
      replicate$seed[1]
    } else {
      as.integer(
        sub(
          "seed_([0-9]+).*",
          "\\1",
          coordinate_file
        )
      )
    }

    aligned_table <- data.frame(
      method = method,
      seed = seed,
      sample_id = replicate$sample_id,
      aligned_PC1 = alignment$aligned[, 1],
      aligned_PC2 = alignment$aligned[, 2],
      region = if ("region" %in% colnames(replicate)) {
        replicate$region
      } else {
        "Unknown"
      },
      stringsAsFactors = FALSE
    )

    aligned_rows[[length(aligned_rows) + 1]] <- aligned_table

    replicate_distances <- distance_vector(
      configuration
    )

    metric_rows[[length(metric_rows) + 1]] <- data.frame(
      method = method,
      seed = seed,
      procrustes_residual_ss = alignment$residual_ss,
      distance_pearson = cor(
        target_distances,
        replicate_distances,
        method = "pearson"
      ),
      distance_spearman = cor(
        target_distances,
        replicate_distances,
        method = "spearman"
      ),
      stringsAsFactors = FALSE
    )
  }

  aligned_all <- do.call(
    rbind,
    aligned_rows
  )

  metrics <- do.call(
    rbind,
    metric_rows
  )

  coordinate_summary <- do.call(
    rbind,
    lapply(
      split(
        aligned_all,
        aligned_all$sample_id
      ),
      function(sample_table) {
        data.frame(
          method = method,
          sample_id = sample_table$sample_id[1],
          region = sample_table$region[1],
          replicates = nrow(sample_table),
          mean_PC1 = mean(sample_table$aligned_PC1),
          sd_PC1 = sd(sample_table$aligned_PC1),
          q025_PC1 = quantile(
            sample_table$aligned_PC1,
            0.025
          ),
          q975_PC1 = quantile(
            sample_table$aligned_PC1,
            0.975
          ),
          mean_PC2 = mean(sample_table$aligned_PC2),
          sd_PC2 = sd(sample_table$aligned_PC2),
          q025_PC2 = quantile(
            sample_table$aligned_PC2,
            0.025
          ),
          q975_PC2 = quantile(
            sample_table$aligned_PC2,
            0.975
          ),
          stringsAsFactors = FALSE
        )
      }
    )
  )

  method_output <- file.path(
    output_dir,
    method
  )

  dir.create(
    method_output,
    recursive = TRUE,
    showWarnings = FALSE
  )

  write.table(
    aligned_all,
    file.path(
      method_output,
      "aligned_replicate_coordinates.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    metrics,
    file.path(
      method_output,
      "robustness_metrics.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    coordinate_summary,
    file.path(
      method_output,
      "sample_coordinate_stability.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  regions <- sort(
    unique(aligned_all$region)
  )

  region_colors <- setNames(
    scales::hue_pal()(
      length(regions)
    ),
    regions
  )

  cloud_plot <- ggplot(
    aligned_all,
    aes(
      x = aligned_PC1,
      y = aligned_PC2,
      color = region,
      group = sample_id
    )
  ) +
    geom_point(
      alpha = 0.18,
      size = 1.8
    ) +
    geom_point(
      data = coordinate_summary,
      aes(
        x = mean_PC1,
        y = mean_PC2,
        color = region
      ),
      size = 4.5,
      shape = 21,
      fill = "white",
      stroke = 1.2
    ) +
    geom_text(
      data = coordinate_summary,
      aes(
        x = mean_PC1,
        y = mean_PC2,
        label = sample_id
      ),
      color = "black",
      fontface = "bold",
      size = 3.5,
      nudge_y = 0.03,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = region_colors
    ) +
    labs(
      title = sprintf(
        "%s balanced-PCA robustness",
        method
      ),
      subtitle = sprintf(
        "%d balanced resampling replicates aligned to the complete-SNP PCA",
        length(unique(aligned_all$seed))
      ),
      x = "Aligned PC1",
      y = "Aligned PC2",
      color = "Geographic region"
    ) +
    theme_classic(
      base_size = 12
    ) +
    theme(
      plot.title = element_text(
        face = "bold"
      ),
      panel.grid.major = element_line(
        color = "grey92",
        linewidth = 0.4
      ),
      legend.position = "right"
    )

  ggsave(
    file.path(
      method_output,
      "balanced_pca_replicate_cloud.pdf"
    ),
    plot = cloud_plot,
    width = 7,
    height = 5.5,
    dpi = 300,
    device = "pdf"
  )

  metric_long <- rbind(
    data.frame(
      seed = metrics$seed,
      metric = "Distance Pearson correlation",
      value = metrics$distance_pearson
    ),
    data.frame(
      seed = metrics$seed,
      metric = "Distance Spearman correlation",
      value = metrics$distance_spearman
    )
  )

  metric_plot <- ggplot(
    metric_long,
    aes(
      x = metric,
      y = value
    )
  ) +
    geom_boxplot(
      width = 0.55,
      outlier.shape = NA
    ) +
    geom_jitter(
      width = 0.12,
      alpha = 0.35,
      size = 1.5
    ) +
    coord_cartesian(
      ylim = c(-1, 1)
    ) +
    labs(
      title = sprintf(
        "%s robustness relative to the complete-SNP PCA",
        method
      ),
      x = NULL,
      y = "Correlation of pairwise sample distances"
    ) +
    theme_classic(
      base_size = 12
    ) +
    theme(
      plot.title = element_text(
        face = "bold"
      ),
      axis.text.x = element_text(
        angle = 15,
        hjust = 1
      )
    )

  ggsave(
    file.path(
      method_output,
      "distance_correlation_distribution.pdf"
    ),
    plot = metric_plot,
    width = 7,
    height = 4.5,
    dpi = 300,
    device = "pdf"
  )

  metric_summary <- data.frame(
    method = method,
    replicates = nrow(metrics),
    procrustes_residual_median = median(
      metrics$procrustes_residual_ss
    ),
    procrustes_residual_q025 = quantile(
      metrics$procrustes_residual_ss,
      0.025
    ),
    procrustes_residual_q975 = quantile(
      metrics$procrustes_residual_ss,
      0.975
    ),
    distance_pearson_median = median(
      metrics$distance_pearson
    ),
    distance_pearson_q025 = quantile(
      metrics$distance_pearson,
      0.025
    ),
    distance_pearson_q975 = quantile(
      metrics$distance_pearson,
      0.975
    ),
    distance_spearman_median = median(
      metrics$distance_spearman
    ),
    distance_spearman_q025 = quantile(
      metrics$distance_spearman,
      0.025
    ),
    distance_spearman_q975 = quantile(
      metrics$distance_spearman,
      0.975
    )
  )

  write.table(
    metric_summary,
    file.path(
      method_output,
      "robustness_metric_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  list(
    aligned = aligned_all,
    metrics = metrics,
    coordinate_summary = coordinate_summary,
    metric_summary = metric_summary
  )
}

# =========================
# MAIN WORKFLOW
# =========================

start_time <- Sys.time()

cat("=============================================================\n")
cat("BALANCED PCA ROBUSTNESS SUMMARY\n")
cat(sprintf("Start          : %s\n", format(start_time, "%Y-%m-%d %H:%M:%S")))
cat(sprintf("Robustness dir : %s\n", ROBUSTNESS_DIR))
cat(sprintf("Output dir     : %s\n", OUTPUT_DIR))
cat("=============================================================\n\n")

pggb_results <- summarize_method(
  method = "PGGB",
  full_coordinates_file = FULL_PGGB_FILE,
  robustness_dir = ROBUSTNESS_DIR,
  output_dir = OUTPUT_DIR
)

cactus_results <- summarize_method(
  method = "Cactus",
  full_coordinates_file = FULL_CACTUS_FILE,
  robustness_dir = ROBUSTNESS_DIR,
  output_dir = OUTPUT_DIR
)

combined_metric_summary <- rbind(
  pggb_results$metric_summary,
  cactus_results$metric_summary
)

write.table(
  combined_metric_summary,
  file.path(
    OUTPUT_DIR,
    "PGGB_Cactus_robustness_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

combined_metrics <- rbind(
  pggb_results$metrics,
  cactus_results$metrics
)

comparison_plot <- ggplot(
  combined_metrics,
  aes(
    x = method,
    y = distance_spearman,
    fill = method
  )
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    alpha = 0.35,
    size = 1.5
  ) +
  coord_cartesian(
    ylim = c(-1, 1)
  ) +
  labs(
    title = "Balanced-PCA robustness by graph method",
    x = NULL,
    y = "Spearman correlation with complete-PCA distances"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    legend.position = "none"
  )

ggsave(
  file.path(
    OUTPUT_DIR,
    "PGGB_Cactus_robustness_comparison.pdf"
  ),
  plot = comparison_plot,
  width = 6,
  height = 4.5,
  dpi = 300,
  device = "pdf"
)

end_time <- Sys.time()

cat("\n")
print(
  combined_metric_summary,
  row.names = FALSE
)

cat("\n=============================================================\n")
cat("ROBUSTNESS SUMMARY COMPLETED\n")
cat(sprintf("Total runtime : %s\n", end_time - start_time))
cat(sprintf("End           : %s\n", format(end_time, "%Y-%m-%d %H:%M:%S")))
cat("=============================================================\n")
