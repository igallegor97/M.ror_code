#!/usr/bin/env Rscript

# =============================================================================
# make_spatial_genetic_figures_publication.R
#
# Publication-ready figures for the PGGB24 geography-versus-genetics analysis.
#
# Visual corrections:
#   1. Country legends retain their fill colors.
#   2. Coordinate-precision symbols are larger and visually distinct.
#   3. Sequencing technology is represented by point outline type, not a
#      competing color scale.
#   4. geom_smooth() no longer inherits the shape aesthetic.
#   5. IBD log1p axes use explicit, readable breaks.
#   6. Zero-km pairs use a restrained cross symbol.
#   7. Network-map axes are labelled Longitude and Latitude.
#   8. all11 and conservative9 use identical scales, colors, and terminology.
#   9. Map legends use explicit override.aes values.
#
# Usage:
#   Rscript make_spatial_genetic_figures_publication.R \
#       validated_coordinates.tsv metadata.tsv \
#       geographic_genetic_pair_table_final.tsv \
#       mantel_results_final.tsv mrm_coefficients_final.tsv \
#       output_directory
# =============================================================================

required_packages <- c(
  "ggplot2",
  "ggrepel",
  "sf",
  "rnaturalearth",
  "rnaturalearthdata",
  "ggspatial",
  "scales"
)

for (package_name in required_packages) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    stop(sprintf("Required R package not installed: %s", package_name))
  }
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(sf)
  library(rnaturalearth)
  library(ggspatial)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 6) {
  stop(
    paste(
      "Usage: Rscript make_spatial_genetic_figures_publication.R",
      "COORDINATES METADATA PAIRS MANTEL MRM_COEFFICIENTS OUTPUT_DIR"
    )
  )
}

COORDINATE_FILE <- args[1]
METADATA_FILE <- args[2]
PAIR_FILE <- args[3]
MANTEL_FILE <- args[4]
MRM_FILE <- args[5]
OUTPUT_DIR <- args[6]

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

pairs <- read.delim(
  PAIR_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

mantel_results <- read.delim(
  MANTEL_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

mrm_coefficients <- read.delim(
  MRM_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# -----------------------------------------------------------------------------
# Data preparation
# -----------------------------------------------------------------------------

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

metadata_columns_to_add <- metadata[, c(
  "sample_id",
  setdiff(colnames(metadata), colnames(valid))
), drop = FALSE]

valid <- merge(
  valid,
  metadata_columns_to_add,
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

valid <- valid[
  match(
    coordinates$sample_id[
      coordinates$sample_id %in% valid$sample_id
    ],
    valid$sample_id
  ),
  ,
  drop = FALSE
]

valid$coordinate_precision <- factor(
  valid$coordinate_precision,
  levels = c("locality", "municipality", "state")
)

valid$technology <- factor(
  valid$technology,
  levels = c("Illumina", "PacBio")
)

pairs$zero_km_pair <- factor(
  pairs$zero_km_pair,
  levels = c(FALSE, TRUE),
  labels = c(
    "Positive geographic distance",
    "0 km at available spatial resolution"
  )
)

pairs$same_country <- factor(
  pairs$same_country,
  levels = c("same country", "different countries"),
  labels = c("Same country", "Different countries")
)

points_sf <- st_as_sf(
  valid,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)

high_precision <- valid[
  valid$coordinate_precision %in% c("locality", "municipality"),
  ,
  drop = FALSE
]

high_precision_sf <- st_as_sf(
  high_precision,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)

world <- ne_countries(
  scale = "medium",
  returnclass = "sf"
)

latin_america <- world[
  world$continent %in% c("South America", "North America"),
]

# -----------------------------------------------------------------------------
# Publication style and scales
# -----------------------------------------------------------------------------

country_palette <- c(
  "Bolivia" = "#D55E00",
  "Colombia" = "#C9A900",
  "Costa Rica" = "#009E73",
  "Ecuador" = "#00A6D6",
  "Mexico" = "#0072B2",
  "Peru" = "#CC79A7"
)

precision_shapes <- c(
  "locality" = 21,
  "municipality" = 22,
  "state" = 24
)

# Technology is shown by outline line type:
# solid black = Illumina; dark grey = PacBio.
technology_outline <- c(
  "Illumina" = "black",
  "PacBio" = "grey35"
)

pair_palette <- c(
  "Same country" = "#0072B2",
  "Different countries" = "#D55E00"
)

profile_palette <- c(
  "All 11 selected communities" = "#0072B2",
  "Conservative 9 selected communities" = "#E69F00"
)

shared_theme <- theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14,
      margin = margin(b = 4)
    ),
    plot.subtitle = element_text(
      size = 10.5,
      color = "grey30",
      margin = margin(b = 8)
    ),
    axis.title = element_text(size = 11.5),
    axis.text = element_text(size = 10),
    legend.title = element_text(face = "bold", size = 10.5),
    legend.text = element_text(size = 9.5),
    legend.position = "right",
    legend.box = "vertical",
    plot.margin = margin(10, 12, 10, 10)
  )

map_guides <- guides(
  fill = guide_legend(
    order = 1,
    override.aes = list(
      shape = 21,
      size = 4.5,
      color = "black",
      stroke = 0.8,
      alpha = 1
    )
  ),
  shape = guide_legend(
    order = 2,
    override.aes = list(
      size = 4.8,
      fill = "white",
      color = "black",
      stroke = 1.1,
      alpha = 1
    )
  ),
  color = guide_legend(
    order = 3,
    override.aes = list(
      shape = 21,
      size = 4.5,
      fill = "white",
      stroke = 1.3,
      alpha = 1
    )
  )
)

# -----------------------------------------------------------------------------
# Sampling maps
# -----------------------------------------------------------------------------

make_sampling_map <- function(
  point_data,
  point_sf,
  title_text,
  subtitle_text,
  x_limits,
  y_limits,
  file_name,
  width,
  height
) {
  plot_object <- ggplot() +
    geom_sf(
      data = latin_america,
      fill = "grey97",
      color = "grey65",
      linewidth = 0.28
    ) +
    geom_sf(
      data = point_sf,
      aes(
        fill = country,
        shape = coordinate_precision,
        color = technology
      ),
      size = 4.4,
      stroke = 1.15,
      alpha = 0.98
    ) +
    geom_text_repel(
      data = point_data,
      aes(
        x = longitude,
        y = latitude,
        label = sample_id
      ),
      size = 3,
      max.overlaps = Inf,
      show.legend = FALSE,
      box.padding = 0.35,
      point.padding = 0.25,
      min.segment.length = 0
    ) +
    scale_fill_manual(
      values = country_palette,
      name = "Country",
      drop = FALSE
    ) +
    scale_shape_manual(
      values = precision_shapes,
      name = "Coordinate precision",
      drop = TRUE
    ) +
    scale_color_manual(
      values = technology_outline,
      name = "Sequencing technology",
      drop = FALSE
    ) +
    map_guides +
    annotation_north_arrow(
      location = "bl",
      which_north = "true",
      pad_x = grid::unit(0.15, "in"),
      pad_y = grid::unit(0.18, "in"),
      style = north_arrow_fancy_orienteering
    ) +
    coord_sf(
      xlim = x_limits,
      ylim = y_limits,
      expand = FALSE
    ) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = "Longitude",
      y = "Latitude"
    ) +
    shared_theme

  ggsave(
    file.path(OUTPUT_DIR, file_name),
    plot_object,
    width = width,
    height = height,
    dpi = 300,
    device = "pdf"
  )

  ggsave(
    file.path(
      OUTPUT_DIR,
      sub("\\.pdf$", ".png", file_name)
    ),
    plot_object,
    width = width,
    height = height,
    dpi = 300
  )
}

make_sampling_map(
  point_data = valid,
  point_sf = points_sf,
  title_text = "Geographic distribution of Moniliophthora roreri isolates",
  subtitle_text = paste(
    "Fill indicates country; shape indicates coordinate precision;",
    "outline indicates sequencing technology."
  ),
  x_limits = c(-105, -55),
  y_limits = c(-25, 25),
  file_name = "Fig1_sampling_map_Latin_America_publication.pdf",
  width = 10.5,
  height = 7.2
)

make_sampling_map(
  point_data = high_precision,
  point_sf = high_precision_sf,
  title_text = "High-precision sampling locations",
  subtitle_text = "Locality- and municipality-level coordinates only.",
  x_limits = c(-90, -65),
  y_limits = c(-20, 15),
  file_name = "FigS1_sampling_map_high_precision_publication.pdf",
  width = 9.2,
  height = 7.2
)

# Regional panel: include only samples from Colombia and Ecuador.
# This prevents labels from out-of-frame samples being pushed onto the
# panel boundary by geom_text_repel().
regional_valid <- droplevels(
  valid[
    valid$country %in% c("Colombia", "Ecuador"),
    ,
    drop = FALSE
  ]
)

regional_points_sf <- st_as_sf(
  regional_valid,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)

make_sampling_map(
  point_data = regional_valid,
  point_sf = regional_points_sf,
  title_text = "Regional sampling view: Colombia and Ecuador",
  subtitle_text = paste0(
    "Fill indicates country; shape indicates coordinate precision;\n",
    "outline indicates sequencing technology."
  ),
  x_limits = c(-80.5, -72),
  y_limits = c(-5.5, 13),
  file_name = "Fig2_sampling_map_Colombia_Ecuador_publication.pdf",
  width = 8.8,
  height = 8.2
)

# -----------------------------------------------------------------------------
# Isolation-by-distance plots
# -----------------------------------------------------------------------------

shared_x_limits <- range(
  pairs$geographic_distance_km,
  na.rm = TRUE
)

shared_y_limits <- range(
  c(
    pairs$genetic_distance_all11,
    pairs$genetic_distance_conservative9
  ),
  na.rm = TRUE
)

shared_y_padding <- diff(shared_y_limits) * 0.04

shared_y_limits <- c(
  max(0, shared_y_limits[1] - shared_y_padding),
  shared_y_limits[2] + shared_y_padding
)

ibd_breaks <- c(0, 50, 100, 250, 500, 1000, 2000, 4000)

ibd_breaks <- ibd_breaks[
  ibd_breaks >= shared_x_limits[1] &
  ibd_breaks <= shared_x_limits[2]
]

extract_annotation <- function(profile) {
  mantel_row <- mantel_results[
    mantel_results$profile == profile &
    mantel_results$subset == "all_valid_coordinates" &
    mantel_results$method == "pearson",
    ,
    drop = FALSE
  ]

  geographic_terms <- mrm_coefficients[
    mrm_coefficients$profile == profile &
    mrm_coefficients$subset == "all_valid_coordinates" &
    grepl(
      "geographic|log10",
      mrm_coefficients$term,
      ignore.case = TRUE
    ),
    ,
    drop = FALSE
  ]

  mantel_text <- if (nrow(mantel_row) > 0) {
    sprintf(
      "Mantel r = %.3f; p = %.4f",
      mantel_row$statistic[1],
      mantel_row$p_value[1]
    )
  } else {
    "Mantel result unavailable"
  }

  mrm_text <- if (nrow(geographic_terms) > 0) {
    sprintf(
      "MRM geographic beta = %.4f; p = %.4f",
      geographic_terms$estimate[1],
      geographic_terms$p_value[1]
    )
  } else {
    "MRM geographic coefficient unavailable"
  }

  paste(mantel_text, mrm_text, sep = "\n")
}

make_ibd_plot <- function(
  response_column,
  profile,
  profile_label,
  file_name
) {
  annotation_text <- extract_annotation(profile)

  plot_object <- ggplot() +
    geom_point(
      data = pairs,
      aes(
        x = geographic_distance_km,
        y = .data[[response_column]],
        color = same_country,
        shape = zero_km_pair
      ),
      alpha = 0.72,
      size = 2.5,
      stroke = 0.85
    ) +
    geom_smooth(
      data = pairs,
      mapping = aes(
        x = geographic_distance_km,
        y = .data[[response_column]],
        group = 1
      ),
      inherit.aes = FALSE,
      method = "lm",
      formula = y ~ log10(x + 1),
      se = TRUE,
      color = "black",
      fill = "grey75",
      alpha = 0.45,
      linewidth = 0.85
    ) +
    scale_x_continuous(
      trans = "log1p",
      breaks = ibd_breaks,
      labels = label_number(
        accuracy = 1,
        big.mark = ","
      ),
      limits = shared_x_limits,
      expand = expansion(mult = c(0.025, 0.05))
    ) +
    scale_y_continuous(
      limits = shared_y_limits,
      breaks = pretty_breaks(n = 6),
      expand = expansion(mult = c(0.01, 0.04))
    ) +
    scale_color_manual(
      values = pair_palette,
      name = "Pair category",
      drop = FALSE
    ) +
    scale_shape_manual(
      values = c(
        "Positive geographic distance" = 16,
        "0 km at available spatial resolution" = 4
      ),
      name = "Spatial resolution",
      drop = FALSE
    ) +
    guides(
      color = guide_legend(
        order = 1,
        override.aes = list(
          shape = 16,
          size = 3,
          alpha = 1
        )
      ),
      shape = guide_legend(
        order = 2,
        override.aes = list(
          color = "black",
          size = 3.4,
          stroke = 1
        )
      )
    ) +
    annotate(
      "label",
      x = Inf,
      y = Inf,
      label = annotation_text,
      hjust = 1.03,
      vjust = 1.15,
      size = 3.25,
      label.size = 0.2,
      fill = alpha("white", 0.88)
    ) +
    labs(
      title = paste(
        "Isolation by distance -",
        profile_label
      ),
      subtitle = sprintf(
        paste(
          "%d isolate pairs; %d zero-km pairs.",
          "Zero km indicates co-location at the available",
          "municipal/locality resolution."
        ),
        nrow(pairs),
        sum(
          pairs$zero_km_pair ==
          "0 km at available spatial resolution"
        )
      ),
      x = "Geographic distance (km; log1p scale)",
      y = "Genetic distance"
    ) +
    shared_theme +
    theme(
      legend.position = "right",
      panel.grid.major.y = element_line(
        color = "grey92",
        linewidth = 0.35
      )
    )

  ggsave(
    file.path(OUTPUT_DIR, file_name),
    plot_object,
    width = 8.5,
    height = 6.3,
    dpi = 300,
    device = "pdf"
  )

  ggsave(
    file.path(
      OUTPUT_DIR,
      sub("\\.pdf$", ".png", file_name)
    ),
    plot_object,
    width = 8.5,
    height = 6.3,
    dpi = 300
  )
}

make_ibd_plot(
  "genetic_distance_all11",
  "all11",
  "all 11 selected communities",
  "Fig3A_IBD_all11_publication.pdf"
)

make_ibd_plot(
  "genetic_distance_conservative9",
  "conservative9",
  "conservative 9 selected communities",
  "Fig3B_IBD_conservative9_publication.pdf"
)

# -----------------------------------------------------------------------------
# Network maps
# -----------------------------------------------------------------------------

coordinate_lookup <- valid[
  ,
  c("sample_id", "longitude", "latitude"),
  drop = FALSE
]

add_coordinates <- function(dataframe) {
  coordinate_1 <- coordinate_lookup[
    match(dataframe$sample_1, coordinate_lookup$sample_id),
    ,
    drop = FALSE
  ]

  coordinate_2 <- coordinate_lookup[
    match(dataframe$sample_2, coordinate_lookup$sample_id),
    ,
    drop = FALSE
  ]

  dataframe$longitude_1 <- coordinate_1$longitude
  dataframe$latitude_1 <- coordinate_1$latitude
  dataframe$longitude_2 <- coordinate_2$longitude
  dataframe$latitude_2 <- coordinate_2$latitude

  dataframe$identical_endpoints <- (
    abs(dataframe$longitude_1 - dataframe$longitude_2) < 1e-12 &
    abs(dataframe$latitude_1 - dataframe$latitude_2) < 1e-12
  )

  dataframe
}

nearest_neighbor_edges <- function(response_column) {
  selected_edges <- list()

  for (sample_id in valid$sample_id) {
    sample_pairs <- pairs[
      pairs$sample_1 == sample_id |
      pairs$sample_2 == sample_id,
      ,
      drop = FALSE
    ]

    sample_pairs <- sample_pairs[
      order(sample_pairs[[response_column]]),
      ,
      drop = FALSE
    ]

    if (nrow(sample_pairs) > 0) {
      selected_edges[[length(selected_edges) + 1]] <-
        sample_pairs[1, , drop = FALSE]
    }
  }

  add_coordinates(
    unique(do.call(rbind, selected_edges))
  )
}

top_similarity_edges <- function(response_column) {
  threshold <- quantile(
    pairs[[response_column]],
    0.10,
    na.rm = TRUE
  )

  add_coordinates(
    pairs[
      pairs[[response_column]] <= threshold,
      ,
      drop = FALSE
    ]
  )
}

plot_network <- function(
  edges,
  response_column,
  title_text,
  subtitle_text,
  file_name
) {
  drawable_edges <- edges[
    !edges$identical_endpoints,
    ,
    drop = FALSE
  ]

  collocated_edges <- edges[
    edges$identical_endpoints,
    ,
    drop = FALSE
  ]

  similarity <- 1 - drawable_edges[[response_column]]

  if (length(similarity) > 0) {
    drawable_edges$similarity_scaled <- rescale(
      similarity,
      to = c(0.4, 1)
    )
  }

  plot_object <- ggplot() +
    geom_sf(
      data = latin_america,
      fill = "grey97",
      color = "grey70",
      linewidth = 0.25
    )

  if (nrow(drawable_edges) > 0) {
    plot_object <- plot_object +
      geom_curve(
        data = drawable_edges,
        aes(
          x = longitude_1,
          y = latitude_1,
          xend = longitude_2,
          yend = latitude_2,
          linewidth = similarity_scaled,
          alpha = similarity_scaled
        ),
        curvature = 0.10,
        color = "grey20",
        lineend = "round"
      )
  }

  if (nrow(collocated_edges) > 0) {
    plot_object <- plot_object +
      geom_point(
        data = collocated_edges,
        aes(
          x = longitude_1,
          y = latitude_1
        ),
        shape = 1,
        size = 7,
        stroke = 1.1,
        color = "black"
      )
  }

  plot_object <- plot_object +
    geom_sf(
      data = points_sf,
      aes(
        fill = country,
        shape = coordinate_precision,
        color = technology
      ),
      size = 4.1,
      stroke = 1.05,
      alpha = 0.98
    ) +
    geom_text_repel(
      data = valid,
      aes(
        x = longitude,
        y = latitude,
        label = sample_id
      ),
      size = 2.8,
      max.overlaps = Inf,
      show.legend = FALSE,
      box.padding = 0.3,
      point.padding = 0.2,
      min.segment.length = 0
    ) +
    scale_fill_manual(
      values = country_palette,
      name = "Country",
      drop = FALSE
    ) +
    scale_shape_manual(
      values = precision_shapes,
      name = "Coordinate precision",
      drop = TRUE
    ) +
    scale_color_manual(
      values = technology_outline,
      name = "Sequencing technology",
      drop = FALSE
    ) +
    scale_linewidth_continuous(
      name = "Relative genetic similarity",
      range = c(0.5, 3.2),
      breaks = c(0.4, 0.6, 0.8, 1.0)
    ) +
    scale_alpha_continuous(
      name = "Relative genetic similarity",
      range = c(0.3, 0.85),
      guide = "none"
    ) +
    guides(
      fill = guide_legend(
        order = 1,
        override.aes = list(
          shape = 21,
          size = 4.3,
          color = "black",
          stroke = 0.8,
          alpha = 1
        )
      ),
      shape = guide_legend(
        order = 2,
        override.aes = list(
          size = 4.7,
          fill = "white",
          color = "black",
          stroke = 1.1,
          alpha = 1
        )
      ),
      color = guide_legend(
        order = 3,
        override.aes = list(
          shape = 21,
          size = 4.3,
          fill = "white",
          stroke = 1.2,
          alpha = 1
        )
      ),
      linewidth = guide_legend(order = 4)
    ) +
    coord_sf(
      xlim = c(-105, -55),
      ylim = c(-25, 25),
      expand = FALSE
    ) +
    labs(
      title = title_text,
      subtitle = paste(
        subtitle_text,
        "Open circles denote linked pairs with identical available coordinates."
      ),
      x = "Longitude",
      y = "Latitude"
    ) +
    shared_theme

  ggsave(
    file.path(OUTPUT_DIR, file_name),
    plot_object,
    width = 10.5,
    height = 7.2,
    dpi = 300,
    device = "pdf"
  )

  ggsave(
    file.path(
      OUTPUT_DIR,
      sub("\\.pdf$", ".png", file_name)
    ),
    plot_object,
    width = 10.5,
    height = 7.2,
    dpi = 300
  )
}

nearest_all11 <- nearest_neighbor_edges(
  "genetic_distance_all11"
)

nearest_conservative9 <- nearest_neighbor_edges(
  "genetic_distance_conservative9"
)

top_all11 <- top_similarity_edges(
  "genetic_distance_all11"
)

top_conservative9 <- top_similarity_edges(
  "genetic_distance_conservative9"
)

write.table(
  nearest_all11,
  file.path(
    OUTPUT_DIR,
    "nearest_neighbor_edges_all11_publication.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  nearest_conservative9,
  file.path(
    OUTPUT_DIR,
    "nearest_neighbor_edges_conservative9_publication.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  top_all11,
  file.path(
    OUTPUT_DIR,
    "top10pct_similarity_edges_all11_publication.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  top_conservative9,
  file.path(
    OUTPUT_DIR,
    "top10pct_similarity_edges_conservative9_publication.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

plot_network(
  nearest_all11,
  "genetic_distance_all11",
  "Nearest genetic neighbor - all 11 selected communities",
  "One nearest genetic neighbor is retained for each isolate.",
  "Fig4A_nearest_neighbor_map_all11_publication.pdf"
)

plot_network(
  nearest_conservative9,
  "genetic_distance_conservative9",
  "Nearest genetic neighbor - conservative 9 selected communities",
  "One nearest genetic neighbor is retained for each isolate.",
  "Fig4B_nearest_neighbor_map_conservative9_publication.pdf"
)

plot_network(
  top_all11,
  "genetic_distance_all11",
  "Top 10% most genetically similar pairs - all 11 selected communities",
  "Edges represent the lowest decile of genetic distances.",
  "Fig5A_top_similarity_network_all11_publication.pdf"
)

plot_network(
  top_conservative9,
  "genetic_distance_conservative9",
  "Top 10% most genetically similar pairs - conservative 9 selected communities",
  "Edges represent the lowest decile of genetic distances.",
  "Fig5B_top_similarity_network_conservative9_publication.pdf"
)

# -----------------------------------------------------------------------------
# Residual mismatch plots
# -----------------------------------------------------------------------------

plot_model_residuals <- function(
  residual_column,
  interpretation_column,
  profile_label,
  file_name
) {
  ordered_pairs <- pairs[
    order(pairs[[residual_column]]),
    ,
    drop = FALSE
  ]

  selected <- rbind(
    head(ordered_pairs, 5),
    tail(ordered_pairs, 5)
  )

  selected$pair <- paste(
    selected$sample_1,
    selected$sample_2,
    sep = " - "
  )

  selected$pair <- factor(
    selected$pair,
    levels = selected$pair[
      order(selected[[residual_column]])
    ]
  )

  residual_palette <- c(
    "genetically closer than fitted distance-decay model" =
      "#0072B2",
    "genetically farther than fitted distance-decay model" =
      "#D55E00"
  )

  residual_labels <- c(
    "genetically closer than fitted distance-decay model" =
      "Closer than fitted model",
    "genetically farther than fitted distance-decay model" =
      "Farther than fitted model"
  )

  plot_object <- ggplot(
    selected,
    aes(
      x = pair,
      y = .data[[residual_column]],
      fill = .data[[interpretation_column]]
    )
  ) +
    geom_col(
      width = 0.72,
      color = "grey20",
      linewidth = 0.25
    ) +
    coord_flip(clip = "off") +
    geom_hline(
      yintercept = 0,
      linewidth = 0.55
    ) +
    scale_fill_manual(
      values = residual_palette,
      labels = residual_labels,
      name = "Residual interpretation"
    ) +
    guides(
      fill = guide_legend(
        nrow = 1,
        byrow = TRUE,
        title.position = "top"
      )
    ) +
    labs(
      title = paste(
        "Genetic-geographic model residuals -",
        profile_label
      ),
      subtitle = paste0(
        "Negative residuals indicate lower genetic distance than fitted;\n",
        "positive residuals indicate higher genetic distance than fitted."
      ),
      x = NULL,
      y = "Standardized residual"
    ) +
    shared_theme +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      plot.subtitle = element_text(
        size = 10.2,
        lineheight = 1.05,
        margin = margin(b = 9)
      ),
      plot.margin = margin(12, 18, 12, 12)
    )

  ggsave(
    file.path(OUTPUT_DIR, file_name),
    plot_object,
    width = 10.2,
    height = 6.5,
    dpi = 300,
    device = "pdf"
  )

  ggsave(
    file.path(
      OUTPUT_DIR,
      sub("\\.pdf$", ".png", file_name)
    ),
    plot_object,
    width = 10.2,
    height = 6.5,
    dpi = 300
  )
}

plot_model_residuals(
  "standardized_residual_all11",
  "residual_interpretation_all11",
  "all 11 selected communities",
  "Fig7A_model_residual_pairs_all11_publication.pdf"
)

plot_model_residuals(
  "standardized_residual_conservative9",
  "residual_interpretation_conservative9",
  "conservative 9 selected communities",
  "Fig7B_model_residual_pairs_conservative9_publication.pdf"
)

# -----------------------------------------------------------------------------
# Ordered distance-class plot
# -----------------------------------------------------------------------------

last_distance_level <- grep(
  "^\\(1735,",
  unique(as.character(pairs$distance_class)),
  value = TRUE
)[1]

pairs$distance_class <- factor(
  pairs$distance_class,
  levels = c(
    "[0, 246]",
    "(246, 780]",
    "(780, 1099]",
    "(1099, 1735]",
    last_distance_level
  ),
  ordered = TRUE
)

distance_long <- rbind(
  data.frame(
    distance_class = pairs$distance_class,
    profile = "All 11 selected communities",
    genetic_distance = pairs$genetic_distance_all11
  ),
  data.frame(
    distance_class = pairs$distance_class,
    profile = "Conservative 9 selected communities",
    genetic_distance =
      pairs$genetic_distance_conservative9
  )
)

distance_class_plot <- ggplot(
  distance_long,
  aes(
    x = distance_class,
    y = genetic_distance,
    fill = profile
  )
) +
  geom_boxplot(
    outlier.alpha = 0.38,
    outlier.size = 1.5,
    width = 0.66,
    position = position_dodge(width = 0.76),
    color = "grey20",
    linewidth = 0.35
  ) +
  scale_fill_manual(
    values = profile_palette,
    name = "Profile"
  ) +
  scale_y_continuous(
    limits = shared_y_limits,
    breaks = pretty_breaks(n = 6)
  ) +
  labs(
    title = "Genetic distance across ordered geographic-distance classes",
    subtitle = paste(
      "Intervals are ordered numerically.",
      "The first class includes zero-km pairs at available spatial resolution."
    ),
    x = "Geographic-distance class (km)",
    y = "Genetic distance"
  ) +
  shared_theme +
  theme(
    axis.text.x = element_text(
      angle = 25,
      hjust = 1
    ),
    legend.position = "bottom"
  )

ggsave(
  file.path(
    OUTPUT_DIR,
    "Fig8_distance_class_boxplots_publication.pdf"
  ),
  distance_class_plot,
  width = 9.7,
  height = 5.9,
  dpi = 300,
  device = "pdf"
)

ggsave(
  file.path(
    OUTPUT_DIR,
    "Fig8_distance_class_boxplots_publication.png"
  ),
  distance_class_plot,
  width = 9.7,
  height = 5.9,
  dpi = 300
)

cat("Final polished spatial-genetic figures completed.\n")
