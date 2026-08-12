#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop("Missing value after ", flag)
  args[[i + 1L]]
}

master_path <- get_arg("--master")
output_dir <- get_arg("--output-dir")
top_n <- as.integer(get_arg("--top-n", "20"))
short_threshold <- as.integer(get_arg("--short-threshold", "100"))
png_dpi <- as.integer(get_arg("--png-dpi", "300"))
go_obo <- get_arg("--go-obo", "")

if (is.null(master_path) || is.null(output_dir)) {
  stop("Usage: 09_plot_functional_annotation_PGGB.R --master TABLE --output-dir DIR [--top-n 20] [--go-obo go-basic.obo]")
}
if (!file.exists(master_path) || file.info(master_path)$size == 0) stop("Master table is missing or empty: ", master_path)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
plot_dir <- file.path(output_dir, "plots")
data_dir <- file.path(output_dir, "plot_data")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

samples_expected <- c("MrorB3", "MrorC26", "MrorCO8", "MrorCO84", "MrorE7")
required <- c("sample_id", "protein_id", "protein_length", "short_protein_flag",
              "eggnog_description", "GO_terms", "KEGG_ko", "pfam_domains", "cazy_families",
              "cazy_number_of_tools", "signalp_positive", "secretome_candidate",
              "effector_candidate", "effectorp_primary_class", "effectorp_dual_localized")

dt <- fread(master_path, sep = "\t", na.strings = c("", "-"), quote = "", showProgress = TRUE)
missing <- setdiff(required, names(dt))
if (length(missing)) stop("Final master is missing columns: ", paste(missing, collapse = ", "))
if (anyDuplicated(dt[, .(sample_id, protein_id)])) stop("Duplicated sample_id/protein_id keys detected")
if (!all(samples_expected %in% unique(dt$sample_id))) stop("One or more expected proteomes are absent")

dt[, sample_id := factor(sample_id, levels = samples_expected)]
dt[, protein_length := as.numeric(protein_length)]
if (anyNA(dt$protein_length) || any(dt$protein_length < 1)) stop("Invalid protein lengths detected")

as_bool <- function(x) toupper(trimws(as.character(x))) == "TRUE"
for (column in c("short_protein_flag", "signalp_positive", "secretome_candidate",
                 "effector_candidate", "effectorp_dual_localized")) {
  set(dt, j = column, value = as_bool(dt[[column]]))
}

if (any(dt$short_protein_flag != (dt$protein_length < short_threshold))) {
  stop("short_protein_flag disagrees with the configured threshold")
}
if (any(dt$secretome_candidate & !dt$signalp_positive)) stop("Secretome includes SignalP-negative proteins")
if (any(dt$effector_candidate & !dt$secretome_candidate)) stop("Effector set includes non-secretome proteins")

theme_pggb <- function(base_size = 10) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          plot.title = element_text(face = "bold", size = rel(1.25)),
          plot.subtitle = element_text(color = "grey35"),
          axis.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom", legend.title = element_text(face = "bold"),
          plot.caption = element_text(color = "grey40", hjust = 0),
          plot.margin = margin(8, 12, 8, 8))
}
theme_set(theme_pggb())

pal <- c(
  blue = "#0072B2", sky = "#56B4E9", green = "#009E73", orange = "#E69F00",
  vermillion = "#D55E00", purple = "#CC79A7", yellow = "#F0E442", grey = "#8A8A8A"
)

save_plot <- function(plot, stem, width, height) {
  for (extension in c("pdf", "svg", "png")) {
    path <- file.path(plot_dir, paste0(stem, ".", extension))
    if (extension == "png") {
      ggsave(path, plot, width = width, height = height, units = "in", dpi = png_dpi, bg = "white")
    } else if (extension == "svg") {
      ggsave(path, plot, width = width, height = height, units = "in", device = svglite::svglite, bg = "white")
    } else {
      ggsave(path, plot, width = width, height = height, units = "in", device = "pdf", bg = "white")
    }
  }
}

write_plot_data <- function(x, filename) fwrite(x, file.path(data_dir, filename), sep = "\t", na = "")

split_annotations <- function(data, column, separator = "[;,]") {
  x <- data[!is.na(get(column)) & nzchar(get(column)), .(sample_id, protein_id, value = get(column))]
  if (!nrow(x)) return(data.table(sample_id = factor(levels = samples_expected), protein_id = character(), term = character()))
  x[, term := strsplit(value, separator, perl = TRUE)]
  x <- x[, .(term = trimws(unlist(term))), by = .(sample_id, protein_id)]
  x[nzchar(term) & term != "-"]
}

top_term_plot <- function(long, title, subtitle, x_label, stem, clean = identity) {
  if (!nrow(long)) stop("No valid annotations were found for ", stem)
  counts <- unique(long[, .(sample_id, protein_id, term)])[, .N, by = .(sample_id, term)]
  totals <- counts[, .(total = sum(N)), by = term][order(-total, term)]
  totals <- head(totals, top_n)
  selected <- counts[term %in% totals$term]
  selected[, term_label := clean(term)]
  order_labels <- totals[, clean(term)]
  selected[, term_label := factor(term_label, levels = rev(order_labels))]
  write_plot_data(selected, paste0(stem, "_counts.tsv"))
  p <- ggplot(selected, aes(N, term_label, fill = sample_id)) +
    geom_col(position = position_dodge2(width = 0.85, preserve = "single"), width = 0.72) +
    scale_fill_manual(values = c(pal["blue"], pal["orange"], pal["green"], pal["purple"], pal["vermillion"]), drop = FALSE) +
    scale_x_continuous(labels = comma, expand = expansion(mult = c(0, 0.04))) +
    labs(title = title, subtitle = subtitle, x = x_label, y = NULL, fill = "Proteome",
         caption = paste0("Top ", min(top_n, nrow(totals)), " terms ranked by the sum across proteomes; proteins are counted once per term."))
  save_plot(p, stem, 10, max(6, 0.30 * nrow(totals) + 2.4))
  p
}

# 01: Pfam families/domains.
pfam <- split_annotations(dt, "pfam_domains")
p_pfam <- top_term_plot(pfam, "Most represented Pfam domains",
                         "Independent protein counts for each proteome", "Proteins", "figure_01_top_pfam")

# 02: consensus CAZy families.
cazy <- split_annotations(dt, "cazy_families", "[;,|]")
cazy[, term := sub("\\(.*$", "", term)]
cazy <- cazy[grepl("^(GH|GT|PL|CE|AA|CBM)[0-9]+", term)]
p_cazy <- top_term_plot(cazy, "Most represented consensus CAZy families",
                         "Only proteins supported by at least two dbCAN tools", "Proteins", "figure_02_top_cazy")

# 03: GO terms, with optional OBO names.
go <- split_annotations(dt, "GO_terms", "[,;]")
go <- go[grepl("^GO:[0-9]{7}$", term)]
go_names <- data.table(term = character(), name = character())
if (nzchar(go_obo)) {
  if (!file.exists(go_obo)) stop("GO OBO file does not exist: ", go_obo)
  lines <- readLines(go_obo, warn = FALSE)
  ids <- grep("^id: GO:[0-9]{7}$", lines)
  names_at <- grep("^name: ", lines)
  name_position <- findInterval(ids, names_at) + 1L
  valid <- name_position <= length(names_at)
  next_id <- c(ids[-1L], length(lines) + 1L)
  valid <- valid & names_at[pmin(name_position, length(names_at))] < next_id
  go_names <- data.table(
    term = sub("^id: ", "", lines[ids[valid]]),
    name = sub("^name: ", "", lines[names_at[name_position[valid]]])
  )
}
go_clean <- function(x) {
  if (!nrow(go_names)) return(x)
  labels <- go_names$name[match(x, go_names$term)]
  ifelse(is.na(labels), x, paste0(labels, " (", x, ")"))
}
p_go <- top_term_plot(go, "Most represented GO terms",
                       if (nrow(go_names)) "GO names obtained from the supplied go-basic.obo" else "GO identifiers shown because no ontology file was supplied",
                       "Proteins", "figure_03_top_go", go_clean)

# 04: KEGG orthology identifiers.
ko <- split_annotations(dt, "KEGG_ko", "[,;]")
ko[, term := sub("^ko:", "", term, ignore.case = TRUE)]
ko <- ko[grepl("^K[0-9]{5}$", term)]
p_ko <- top_term_plot(ko, "Most represented KEGG Orthology identifiers",
                       "Independent protein counts for each proteome", "Proteins", "figure_04_top_kegg_ko")

# 05: length distributions and short-protein counts.
length_plot_data <- copy(dt[, .(sample_id, protein_id, protein_length, short_protein_flag)])
write_plot_data(length_plot_data, "figure_05_protein_lengths.tsv")
p_length_a <- ggplot(length_plot_data, aes(protein_length, colour = sample_id, fill = sample_id)) +
  geom_density(alpha = 0.10, linewidth = 0.75, adjust = 1.1) +
  geom_vline(xintercept = short_threshold, linetype = 2, colour = pal["vermillion"]) +
  coord_cartesian(xlim = c(0, quantile(length_plot_data$protein_length, 0.99))) +
  scale_colour_manual(values = c(pal["blue"], pal["orange"], pal["green"], pal["purple"], pal["vermillion"])) +
  scale_fill_manual(values = c(pal["blue"], pal["orange"], pal["green"], pal["purple"], pal["vermillion"])) +
  labs(title = "Protein-length distributions", subtitle = "Displayed through the 99th percentile so extreme proteins do not compress the distribution",
       x = "Protein length (amino acids)", y = "Density", colour = "Proteome", fill = "Proteome")
short_counts <- dt[, .(total = .N, short = sum(short_protein_flag)), by = sample_id]
short_counts[, percent := 100 * short / total]
p_length_b <- ggplot(short_counts, aes(sample_id, percent, fill = sample_id)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = paste0(comma(short), "\n", number(percent, accuracy = 0.1), "%")), vjust = -0.25, size = 3) +
  scale_fill_manual(values = c(pal["blue"], pal["orange"], pal["green"], pal["purple"], pal["vermillion"])) +
  scale_y_continuous(labels = label_percent(scale = 1), expand = expansion(mult = c(0, 0.14))) +
  labs(title = paste0("Short proteins (<", short_threshold, " aa)"), x = NULL, y = "Percentage of proteome")
p_length <- p_length_a / p_length_b + plot_layout(heights = c(1.35, 1)) + plot_annotation(tag_levels = "A")
save_plot(p_length, "figure_05_protein_length_and_short_proteins", 9, 8)

# Shared per-sample summary.
summary_dt <- dt[, .(
  proteins = .N,
  eggnog_description = sum(!is.na(eggnog_description) & nzchar(eggnog_description)),
  GO = sum(!is.na(GO_terms) & nzchar(GO_terms)),
  KEGG_KO = sum(!is.na(KEGG_ko) & nzchar(KEGG_ko)),
  Pfam = sum(!is.na(pfam_domains) & nzchar(pfam_domains)),
  CAZy_consensus = sum(!is.na(cazy_families) & nzchar(cazy_families)),
  SignalP_positive = sum(signalp_positive),
  Secretome = sum(secretome_candidate),
  Effectors = sum(effector_candidate)
), by = sample_id]
write_plot_data(summary_dt, "functional_coverage_summary.tsv")

# 06: direct secretome/effector comparison.
comparison <- melt(summary_dt, id.vars = c("sample_id", "proteins"),
                   measure.vars = c("SignalP_positive", "Secretome", "Effectors"),
                   variable.name = "stage", value.name = "count")
comparison[, percent := 100 * count / proteins]
comparison[, stage := factor(stage, levels = c("SignalP_positive", "Secretome", "Effectors"),
                             labels = c("SignalP positive", "Final secretome", "Effector candidates"))]
write_plot_data(comparison, "figure_06_secretome_effector_comparison.tsv")
p_compare <- ggplot(comparison, aes(sample_id, percent, fill = stage)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72) +
  geom_text(aes(label = comma(count)), position = position_dodge(width = 0.8), vjust = -0.25, size = 2.8) +
  scale_fill_manual(values = c(pal["sky"], pal["green"], pal["purple"])) +
  scale_y_continuous(labels = label_percent(scale = 1), expand = expansion(mult = c(0, 0.13))) +
  labs(title = "Secretome and effector predictions by proteome", subtitle = "Labels give protein counts; bar height gives percentage of the complete proteome",
       x = NULL, y = "Percentage of proteome", fill = "Analysis stage")
save_plot(p_compare, "figure_06_secretome_effector_comparison", 9, 5.5)

# 07: presence/coverage heatmap.
coverage <- melt(summary_dt, id.vars = c("sample_id", "proteins"),
                 measure.vars = c("eggnog_description", "GO", "KEGG_KO", "Pfam", "CAZy_consensus", "SignalP_positive", "Secretome", "Effectors"),
                 variable.name = "annotation", value.name = "count")
coverage[, percent := 100 * count / proteins]
coverage[, annotation := factor(annotation, levels = rev(c("eggnog_description", "GO", "KEGG_KO", "Pfam", "CAZy_consensus", "SignalP_positive", "Secretome", "Effectors")),
                              labels = rev(c("eggNOG description", "GO", "KEGG KO", "Pfam", "CAZy consensus", "SignalP positive", "Final secretome", "Effectors")))]
write_plot_data(coverage, "figure_07_functional_coverage_heatmap.tsv")
p_heat <- ggplot(coverage, aes(sample_id, annotation, fill = percent)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = paste0(number(percent, accuracy = 0.1), "%\n", comma(count))), size = 3) +
  scale_fill_gradientn(colours = c("#F7FBFF", pal["sky"], pal["blue"]), limits = c(0, 100), labels = label_percent(scale = 1)) +
  coord_equal() +
  labs(title = "Functional annotation presence and coverage", subtitle = "Percentage and protein count in every cell",
       x = NULL, y = NULL, fill = "Proteome coverage")
save_plot(p_heat, "figure_07_functional_coverage_heatmap", 8.5, 6.2)

# 08: sequential filtering flow; only genuinely nested sets are used.
flow <- data.table(
  stage = factor(c("All proteins", "SignalP positive", "Final secretome", "Effector candidates"),
                 levels = c("All proteins", "SignalP positive", "Final secretome", "Effector candidates")),
  count = c(nrow(dt), sum(dt$signalp_positive), sum(dt$secretome_candidate), sum(dt$effector_candidate))
)
flow[, percent_all := 100 * count / count[1]]
write_plot_data(flow, "figure_08_secretome_effector_flow.tsv")
p_flow <- ggplot(flow, aes(stage, count, fill = stage)) +
  geom_col(width = 0.68, show.legend = FALSE) +
  geom_text(aes(label = paste0(comma(count), "\n", number(percent_all, accuracy = 0.1), "% of all proteins")), vjust = -0.25, size = 3.2) +
  scale_fill_manual(values = c(pal["blue"], pal["sky"], pal["green"], pal["purple"])) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.14))) +
  labs(title = "Sequential secretome-to-effector filtering", subtitle = "Combined counts from all five proteomes",
       x = NULL, y = "Proteins")
save_plot(p_flow, "figure_08_secretome_effector_flow", 8.5, 5.5)

# 09: EffectorP primary composition; dual localization is deliberately not stacked.
eff <- dt[secretome_candidate == TRUE, .(
  Non_effector = sum(!effector_candidate),
  Apoplastic_primary = sum(effector_candidate & effectorp_primary_class == "Apoplastic"),
  Cytoplasmic_primary = sum(effector_candidate & effectorp_primary_class == "Cytoplasmic"),
  Dual_localized = sum(effector_candidate & effectorp_dual_localized)
), by = sample_id]
eff_long <- melt(eff, id.vars = c("sample_id", "Dual_localized"),
                 measure.vars = c("Non_effector", "Apoplastic_primary", "Cytoplasmic_primary"),
                 variable.name = "class", value.name = "count")
eff_long[, class := factor(class, levels = c("Non_effector", "Apoplastic_primary", "Cytoplasmic_primary"),
                           labels = c("Non-effector", "Apoplastic-primary", "Cytoplasmic-primary"))]
write_plot_data(eff_long, "figure_09_effectorp_composition.tsv")
p_eff <- ggplot(eff_long, aes(sample_id, count, fill = class)) +
  geom_col(width = 0.72) +
  geom_text(data = eff, aes(sample_id, Non_effector + Apoplastic_primary + Cytoplasmic_primary,
                            label = paste0("dual = ", Dual_localized)), inherit.aes = FALSE, vjust = -0.35, size = 3) +
  scale_fill_manual(values = c(pal["grey"], pal["orange"], pal["purple"])) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.10))) +
  labs(title = "EffectorP 3 composition of the final secretome",
       subtitle = "Primary classes are exclusive; dual-localized proteins overlap them and are annotated above each bar",
       x = NULL, y = "Final secretome proteins", fill = "Primary classification")
save_plot(p_eff, "figure_09_effectorp_composition", 9, 5.5)

# 10: dbCAN evidence support.
dbcan <- dt[!is.na(cazy_number_of_tools), .N, by = .(sample_id, tools = as.integer(cazy_number_of_tools))]
dbcan <- dbcan[tools %in% 1:3]
dbcan[, evidence := factor(tools, levels = 1:3, labels = c("One tool", "Two tools", "Three tools"))]
write_plot_data(dbcan, "figure_10_dbcan_support.tsv")
p_dbcan <- ggplot(dbcan, aes(sample_id, N, fill = evidence)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = c(pal["grey"], pal["orange"], pal["green"])) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) +
  labs(title = "dbCAN support by number of agreeing tools", subtitle = "Two- and three-tool calls constitute the consensus CAZyme set",
       x = NULL, y = "Proteins with dbCAN evidence", fill = "Evidence")
save_plot(p_dbcan, "figure_10_dbcan_evidence_support", 9, 5.5)

# 00: compact multi-panel overview for a manuscript or presentation.
overview <- (p_flow + p_heat) / (p_compare + p_eff) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Functional annotation, secretome and effector overview",
                  subtitle = "Five PacBio proteomes of Moniliophthora roreri",
                  tag_levels = "A") &
  theme(legend.position = "bottom")
save_plot(overview, "figure_00_functional_annotation_overview", 14, 10)

manifest <- data.table(
  figure = c("00", sprintf("%02d", 1:10)),
  description = c(
    "Four-panel functional annotation overview", "Most represented Pfam domains",
    "Most represented consensus CAZy families", "Most represented GO terms",
    "Most represented KEGG Orthology identifiers", "Protein lengths and proteins below 100 aa",
    "SignalP, secretome and effector comparison", "Functional coverage heatmap",
    "Sequential secretome-to-effector flow", "EffectorP primary composition and dual subset",
    "dbCAN evidence support"
  )
)
write_plot_data(manifest, "figure_manifest.tsv")
write_plot_data(data.table(
  input_master = normalizePath(master_path),
  input_md5 = unname(tools::md5sum(master_path)),
  generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  top_n = top_n,
  short_protein_threshold = short_threshold,
  png_dpi = png_dpi,
  go_obo = if (nzchar(go_obo)) normalizePath(go_obo) else ""
), "run_metadata.tsv")
capture.output(sessionInfo(), file = file.path(output_dir, "R_sessionInfo.txt"))

cat("FUNCTIONAL ANNOTATION FIGURES COMPLETED\n")
cat("Proteins:", format(nrow(dt), big.mark = ","), "\n")
cat("Plot directory:", plot_dir, "\n")
cat("Plot-data directory:", data_dir, "\n")
