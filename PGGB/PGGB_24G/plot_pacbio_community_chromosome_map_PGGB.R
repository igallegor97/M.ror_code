#!/usr/bin/env Rscript

# =============================================================
# plot_pacbio_community_chromosome_map.R
#
# Creates a publication-ready two-panel figure showing how the
# chromosome-scale PacBio sequences were assigned to PGGB
# communities after partition-before-pggb.
#
# Panel A:
#   Syntenic groups (Group1-Group11) -> PGGB communities
#
# Panel B:
#   PacBio ungrouped sequences -> PGGB communities
#
# Links are colored by PacBio genome:
#   B3, C26, CO8, CO84, E7
#
# The figure is generated directly from:
#   pacbio_chromosome_to_communities.tsv
#   community_summary.tsv
#
# Outputs:
#   pacbio_chromosome_community_map.pdf
#   pacbio_chromosome_community_map.png
#   panel_A_grouped_chromosomes.pdf
#   panel_B_ungrouped_sequences.pdf
#   figure_edge_data.tsv
#
# Usage:
#   Rscript plot_pacbio_community_chromosome_map.R
#
# Configuration is read from environment variables:
#   MAP_TSV
#   COMMUNITY_SUMMARY
#   OUTPUT_DIR
# =============================================================

required_packages <- c("ggplot2")

for (package_name in required_packages) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    stop(
      sprintf(
        "Required R package is not installed: %s",
        package_name
      )
    )
  }
}

suppressPackageStartupMessages({
  library(ggplot2)
})

# =========================
# CONFIGURATION
# =========================

MAP_TSV <- Sys.getenv(
  "MAP_TSV",
  paste0(
    "/Storage/data1/isabella.gallego/MAESTRIA/",
    "PGGB_results_24G/community_sequence_tables/",
    "pacbio_chromosome_to_communities.tsv"
  )
)

COMMUNITY_SUMMARY <- Sys.getenv(
  "COMMUNITY_SUMMARY",
  paste0(
    "/Storage/data1/isabella.gallego/MAESTRIA/",
    "PGGB_results_24G/community_sequence_tables/",
    "community_summary.tsv"
  )
)

OUTPUT_DIR <- Sys.getenv(
  "OUTPUT_DIR",
  paste0(
    "/Storage/data1/isabella.gallego/MAESTRIA/",
    "PGGB_results_24G/community_sequence_figures"
  )
)

SAMPLE_COLORS <- c(
  B3   = "#0072B2",
  C26  = "#E69F00",
  CO8  = "#009E73",
  CO84 = "#CC79A7",
  E7   = "#56B4E9"
)

SAMPLE_ORDER <- c(
  "B3",
  "C26",
  "CO8",
  "CO84",
  "E7"
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# =========================
# PRE-RUN CHECKS
# =========================

if (!file.exists(MAP_TSV)) {
  stop(
    sprintf(
      "MAP_TSV not found: %s",
      MAP_TSV
    )
  )
}

if (!file.exists(COMMUNITY_SUMMARY)) {
  stop(
    sprintf(
      "COMMUNITY_SUMMARY not found: %s",
      COMMUNITY_SUMMARY
    )
  )
}

start_time <- Sys.time()

cat("=============================================================\n")
cat("PACBIO CHROMOSOME-TO-COMMUNITY FIGURE\n")
cat(sprintf("Start             : %s\n", format(start_time, "%Y-%m-%d %H:%M:%S")))
cat(sprintf("MAP_TSV           : %s\n", MAP_TSV))
cat(sprintf("COMMUNITY_SUMMARY : %s\n", COMMUNITY_SUMMARY))
cat(sprintf("OUTPUT_DIR        : %s\n", OUTPUT_DIR))
cat("=============================================================\n\n")

# =========================
# LOAD TABLES
# =========================

mapping <- read.delim(
  MAP_TSV,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = ""
)

community_summary <- read.delim(
  COMMUNITY_SUMMARY,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = ""
)

required_map_columns <- c(
  "sequence_class",
  "canonical_label",
  "sample_sequence_details"
)

missing_map_columns <- setdiff(
  required_map_columns,
  colnames(mapping)
)

if (length(missing_map_columns) > 0) {
  stop(
    sprintf(
      "Mapping table is missing columns: %s",
      paste(missing_map_columns, collapse = ", ")
    )
  )
}

required_summary_columns <- c(
  "community",
  "n_sequences"
)

missing_summary_columns <- setdiff(
  required_summary_columns,
  colnames(community_summary)
)

if (length(missing_summary_columns) > 0) {
  stop(
    sprintf(
      "Community summary is missing columns: %s",
      paste(missing_summary_columns, collapse = ", ")
    )
  )
}

# =========================
# PARSE SAMPLE-LEVEL EDGES
# =========================

parse_detail_string <- function(
    detail_string,
    sequence_class,
    canonical_label
) {
  if (is.na(detail_string) || detail_string == "") {
    return(NULL)
  }

  details <- strsplit(
    detail_string,
    ";",
    fixed = TRUE
  )[[1]]

  rows <- lapply(
    details,
    function(detail) {
      match_result <- regexec(
        "^([^:]+):([^@]+)@(community\\.[0-9]+)$",
        detail
      )

      captured <- regmatches(
        detail,
        match_result
      )[[1]]

      if (length(captured) != 4) {
        stop(
          sprintf(
            "Could not parse sample_sequence_details entry: %s",
            detail
          )
        )
      }

      data.frame(
        sequence_class = sequence_class,
        canonical_label = canonical_label,
        sample = captured[2],
        contig = captured[3],
        community = captured[4],
        stringsAsFactors = FALSE
      )
    }
  )

  do.call(
    rbind,
    rows
  )
}

edge_tables <- lapply(
  seq_len(nrow(mapping)),
  function(index) {
    parse_detail_string(
      mapping$sample_sequence_details[index],
      mapping$sequence_class[index],
      mapping$canonical_label[index]
    )
  }
)

edges <- do.call(
  rbind,
  edge_tables
)

edges <- edges[
  edges$sample %in% SAMPLE_ORDER,
  ,
  drop = FALSE
]

edges$sample <- factor(
  edges$sample,
  levels = SAMPLE_ORDER
)

# Add community sizes for labels.
edges <- merge(
  edges,
  community_summary[, c("community", "n_sequences")],
  by = "community",
  all.x = TRUE,
  sort = FALSE
)

edges$community_label <- sprintf(
  "%s\n%s sequences",
  edges$community,
  format(
    edges$n_sequences,
    big.mark = ",",
    scientific = FALSE
  )
)

write.table(
  edges,
  file.path(
    OUTPUT_DIR,
    "figure_edge_data.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# =========================
# HELPER FUNCTIONS
# =========================

extract_number <- function(value) {
  as.integer(
    sub(
      ".*?([0-9]+)$",
      "\\1",
      value
    )
  )
}

community_number <- function(value) {
  as.integer(
    sub(
      "community\\.",
      "",
      value
    )
  )
}

build_flow_plot <- function(
    plot_data,
    panel_title,
    panel_subtitle,
    left_title,
    right_title,
    highlight_group10 = FALSE
) {
  if (nrow(plot_data) == 0) {
    stop(
      sprintf(
        "No rows available for panel: %s",
        panel_title
      )
    )
  }

  # Left-node ordering follows the biological numeric order.
  left_labels <- unique(
    plot_data$canonical_label
  )

  left_labels <- left_labels[
    order(
      extract_number(left_labels)
    )
  ]

  # Right-node ordering follows the community number.
  community_lookup <- unique(
    plot_data[, c("community", "community_label")]
  )

  community_lookup <- community_lookup[
    order(
      community_number(
        community_lookup$community
      )
    ),
    ,
    drop = FALSE
  ]

  community_labels <- community_lookup$community_label

  left_y <- setNames(
    rev(seq_along(left_labels)),
    left_labels
  )

  right_y <- setNames(
    rev(seq_along(community_labels)),
    community_labels
  )

  plot_data$left_y <- unname(
    left_y[
      plot_data$canonical_label
    ]
  )

  plot_data$right_y <- unname(
    right_y[
      plot_data$community_label
    ]
  )

  # Small offsets allow all five sample-specific links to remain visible.
  sample_offsets <- setNames(
    seq(
      -0.16,
      0.16,
      length.out = length(SAMPLE_ORDER)
    ),
    SAMPLE_ORDER
  )

  plot_data$offset <- unname(
    sample_offsets[
      as.character(
        plot_data$sample
      )
    ]
  )

  plot_data$y_start <- (
    plot_data$left_y
    + plot_data$offset
  )

  plot_data$y_end <- (
    plot_data$right_y
    + plot_data$offset
  )

  left_nodes <- data.frame(
    canonical_label = left_labels,
    y = unname(
      left_y[
        left_labels
      ]
    ),
    stringsAsFactors = FALSE
  )

  right_nodes <- data.frame(
    community_label = community_labels,
    y = unname(
      right_y[
        community_labels
      ]
    ),
    stringsAsFactors = FALSE
  )

  plot_object <- ggplot() +
    geom_curve(
      data = plot_data,
      aes(
        x = 0.08,
        y = y_start,
        xend = 0.92,
        yend = y_end,
        color = sample,
        group = interaction(
          canonical_label,
          sample,
          community
        )
      ),
      curvature = 0.12,
      linewidth = 1.15,
      alpha = 0.82,
      lineend = "round"
    ) +
    geom_rect(
      data = left_nodes,
      aes(
        xmin = -0.02,
        xmax = 0.08,
        ymin = y - 0.31,
        ymax = y + 0.31
      ),
      fill = "white",
      color = "grey25",
      linewidth = 0.55
    ) +
    geom_rect(
      data = right_nodes,
      aes(
        xmin = 0.92,
        xmax = 1.02,
        ymin = y - 0.31,
        ymax = y + 0.31
      ),
      fill = "white",
      color = "grey25",
      linewidth = 0.55
    ) +
    geom_text(
      data = left_nodes,
      aes(
        x = -0.035,
        y = y,
        label = canonical_label
      ),
      hjust = 1,
      size = 3.7,
      fontface = "bold"
    ) +
    geom_text(
      data = right_nodes,
      aes(
        x = 1.035,
        y = y,
        label = community_label
      ),
      hjust = 0,
      size = 3.15,
      lineheight = 0.9
    ) +
    annotate(
      "text",
      x = -0.035,
      y = max(left_nodes$y) + 0.78,
      label = left_title,
      hjust = 1,
      fontface = "bold",
      size = 4
    ) +
    annotate(
      "text",
      x = 1.035,
      y = max(right_nodes$y) + 0.78,
      label = right_title,
      hjust = 0,
      fontface = "bold",
      size = 4
    ) +
    scale_color_manual(
      values = SAMPLE_COLORS,
      drop = FALSE
    ) +
    coord_cartesian(
      xlim = c(-0.42, 1.42),
      ylim = c(
        0.3,
        max(
          max(left_nodes$y),
          max(right_nodes$y)
        ) + 1.1
      ),
      clip = "off"
    ) +
    labs(
      title = panel_title,
      subtitle = panel_subtitle,
      color = "PacBio genome"
    ) +
    theme_void(
      base_size = 11
    ) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 13,
        hjust = 0
      ),
      plot.subtitle = element_text(
        size = 10,
        color = "grey35",
        margin = margin(
          b = 10
        )
      ),
      legend.position = "bottom",
      legend.title = element_text(
        face = "bold"
      ),
      plot.margin = margin(
        t = 10,
        r = 110,
        b = 10,
        l = 95
      )
    )

  if (
    highlight_group10
    && "Group10" %in% left_nodes$canonical_label
  ) {
    group10_y <- left_nodes$y[
      left_nodes$canonical_label == "Group10"
    ]

    plot_object <- plot_object +
      annotate(
        "rect",
        xmin = -0.055,
        xmax = 0.105,
        ymin = group10_y - 0.40,
        ymax = group10_y + 0.40,
        fill = NA,
        color = "#D55E00",
        linewidth = 1.1,
        linetype = "dashed"
      ) +
      annotate(
        "text",
        x = 0.19,
        y = group10_y + 0.44,
        label = "Split assignment",
        color = "#D55E00",
        fontface = "bold",
        hjust = 0,
        size = 3.2
      )
  }

  plot_object
}

# =========================
# PANEL A: GROUPED CHROMOSOMES
# =========================

grouped_edges <- edges[
  edges$sequence_class == "grouped",
  ,
  drop = FALSE
]

panel_a <- build_flow_plot(
  plot_data = grouped_edges,
  panel_title = "A  Chromosome-scale syntenic groups",
  panel_subtitle = paste(
    "Ten groups map completely to one community;",
    "Group10 is split between community.0 and community.1."
  ),
  left_title = "PacBio syntenic group",
  right_title = "PGGB community",
  highlight_group10 = TRUE
)

# =========================
# PANEL B: UNGROUPED SEQUENCES
# =========================

ungrouped_edges <- edges[
  edges$sequence_class == "ungrouped",
  ,
  drop = FALSE
]

panel_b <- build_flow_plot(
  plot_data = ungrouped_edges,
  panel_title = "B  PacBio ungrouped sequences",
  panel_subtitle = paste(
    "Ungrouped labels are assembly-specific identifiers;",
    "identical numbers do not imply homology."
  ),
  left_title = "Ungrouped sequence label",
  right_title = "PGGB community",
  highlight_group10 = FALSE
)

# =========================
# SAVE INDIVIDUAL PANELS
# =========================

ggsave(
  file.path(
    OUTPUT_DIR,
    "panel_A_grouped_chromosomes.pdf"
  ),
  plot = panel_a,
  width = 10.5,
  height = 8.2,
  dpi = 300,
  device = "pdf"
)

ggsave(
  file.path(
    OUTPUT_DIR,
    "panel_B_ungrouped_sequences.pdf"
  ),
  plot = panel_b,
  width = 10.5,
  height = 7.4,
  dpi = 300,
  device = "pdf"
)

# =========================
# COMBINED FIGURE
# =========================

combined_pdf <- file.path(
  OUTPUT_DIR,
  "pacbio_chromosome_community_map.pdf"
)

pdf(
  combined_pdf,
  width = 11,
  height = 15.5,
  onefile = TRUE
)

grid::grid.newpage()

layout <- grid::grid.layout(
  nrow = 2,
  ncol = 1,
  heights = grid::unit(
    c(0.54, 0.46),
    "npc"
  )
)

grid::pushViewport(
  grid::viewport(
    layout = layout
  )
)

print(
  panel_a,
  vp = grid::viewport(
    layout.pos.row = 1,
    layout.pos.col = 1
  )
)

print(
  panel_b,
  vp = grid::viewport(
    layout.pos.row = 2,
    layout.pos.col = 1
  )
)

grid::popViewport()
dev.off()

combined_png <- file.path(
  OUTPUT_DIR,
  "pacbio_chromosome_community_map.png"
)

png(
  combined_png,
  width = 3300,
  height = 4650,
  res = 300
)

grid::grid.newpage()

layout <- grid::grid.layout(
  nrow = 2,
  ncol = 1,
  heights = grid::unit(
    c(0.54, 0.46),
    "npc"
  )
)

grid::pushViewport(
  grid::viewport(
    layout = layout
  )
)

print(
  panel_a,
  vp = grid::viewport(
    layout.pos.row = 1,
    layout.pos.col = 1
  )
)

print(
  panel_b,
  vp = grid::viewport(
    layout.pos.row = 2,
    layout.pos.col = 1
  )
)

grid::popViewport()
dev.off()

# =========================
# FIGURE SUMMARY
# =========================

grouped_summary <- unique(
  grouped_edges[
    ,
    c(
      "canonical_label",
      "community"
    )
  ]
)

groups_in_one_community <- sum(
  table(
    grouped_summary$canonical_label
  ) == 1
)

groups_in_multiple_communities <- sum(
  table(
    grouped_summary$canonical_label
  ) > 1
)

summary_table <- data.frame(
  metric = c(
    "grouped_syntenic_groups",
    "groups_in_one_community",
    "groups_in_multiple_communities",
    "ungrouped_sequence_records",
    "ungrouped_communities"
  ),
  value = c(
    length(
      unique(
        grouped_edges$canonical_label
      )
    ),
    groups_in_one_community,
    groups_in_multiple_communities,
    nrow(
      ungrouped_edges
    ),
    length(
      unique(
        ungrouped_edges$community
      )
    )
  )
)

write.table(
  summary_table,
  file.path(
    OUTPUT_DIR,
    "figure_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

end_time <- Sys.time()

cat("\n=============================================================\n")
cat("FIGURE COMPLETED\n")
cat(sprintf("Grouped groups             : %d\n", summary_table$value[1]))
cat(sprintf("Groups in one community    : %d\n", summary_table$value[2]))
cat(sprintf("Groups split               : %d\n", summary_table$value[3]))
cat(sprintf("Ungrouped records          : %d\n", summary_table$value[4]))
cat(sprintf("PDF                        : %s\n", combined_pdf))
cat(sprintf("PNG                        : %s\n", combined_png))
cat(sprintf("Total runtime              : %s\n", end_time - start_time))
cat(sprintf("End                        : %s\n", format(end_time, "%Y-%m-%d %H:%M:%S")))
cat("=============================================================\n")
