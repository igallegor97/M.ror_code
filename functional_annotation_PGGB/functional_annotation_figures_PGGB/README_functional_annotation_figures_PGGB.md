# Functional annotation figures for the PGGB project

This package generates reproducible, publication-quality figures from the final integrated functional annotation table for the five PacBio proteomes. Every figure is exported as PDF, SVG and 300-dpi PNG. The numerical data used in each plot are also written as TSV files, allowing every displayed value to be audited independently.

## Input

The principal input is:

```text
$RESULTS/07_master/final_integrated/functional_annotation_master_final.tsv
```

The script expects one row per protein and the following columns:

- `sample_id` and `protein_id`;
- `protein_length` and `short_protein_flag`;
- `eggnog_description`, `GO_terms` and `KEGG_ko`;
- `pfam_domains`;
- `cazy_families` and `cazy_number_of_tools`;
- `signalp_positive`;
- `secretome_candidate`;
- `effector_candidate`, `effectorp_primary_class` and `effectorp_dual_localized`.

The expected proteomes are `MrorB3`, `MrorC26`, `MrorCO8`, `MrorCO84` and `MrorE7`.

An optional `go-basic.obo` file can be supplied through `GO_OBO`. When it is available, the GO figure includes ontology names. Without it, the plot remains valid but uses GO identifiers only.

## Software requirements

The job requires R and these packages:

```r
data.table
ggplot2
scales
patchwork
svglite
```

If the preflight reports missing packages and the cluster permits user-level CRAN installations, run:

```bash
Rscript 00_install_plotting_R_packages_PGGB.R
```

If compute nodes have no internet access, run this once on the permitted login/software-installation node or ask the HPC administrator to provide the packages. Do not repeatedly install packages inside an SGE job.

The plots use a colour-blind-friendly palette and do not rely on colour alone for critical interpretations. Counts and percentages are kept distinct. Dual-localized EffectorP candidates are not stacked as an independent class because they overlap the apoplastic-primary and cytoplasmic-primary classifications.

## Configuration

Edit `plotting_environment_PGGB.env`. Verify `PROJECT_ROOT`, `RESULTS`, `FINAL_MASTER`, `FIGURE_ROOT` and `LOGS`.

Find the installed R module with:

```bash
module spider R
```

Enter the exact module name in `R_MODULE`. Leave it empty only if `Rscript` is already available.

Optionally set:

```bash
GO_OBO=/absolute/path/to/go-basic.obo
TOP_N=20
SHORT_PROTEIN_THRESHOLD=100
PNG_DPI=300
```

`TOP_N` controls the number of Pfam, CAZy, GO and KEGG terms shown. The complete master table is not reduced; this setting affects display only.

## Preflight

Copy or extract the package files directly into `PROJECT_ROOT`. Then, from that directory:

```bash
bash 00_preflight_functional_annotation_figures_PGGB.sh
```

The preflight verifies the master table, required columns, R executable and R packages. Do not submit the figure job until it reports `FIGURE PREFLIGHT COMPLETED`.

## SGE submission

Load the existing pipeline configuration only if it is safe in the current shell. The plotting job reads its own environment file, so sourcing is not required.

```bash
qsub -terse \
  -pe smp 4 \
  -l h_vmem=24G \
  -l h_rt=12:00:00 \
  -o /Storage/data1/isabella.gallego/MAESTRIA/code/functional_annotation_PGGB/logs \
  -e /Storage/data1/isabella.gallego/MAESTRIA/code/functional_annotation_PGGB/logs \
  09_plot_functional_annotation_PGGB.sge.sh
```

Resource names can differ among SGE installations. If `h_vmem` or `h_rt` is not accepted by the cluster, use the same resource syntax as the successfully completed annotation jobs.

Monitor the job with:

```bash
qstat -u "$USER"
```

Inspect the most recent log with:

```bash
ls -lht logs/pggb_func_plots* | head
```

## Outputs

The default output root is:

```text
$RESULTS/08_figures/
```

`plots/` contains each figure in PDF, SVG and PNG:

1. `figure_00_functional_annotation_overview`: four-panel overview.
2. `figure_01_top_pfam`: most represented Pfam domains.
3. `figure_02_top_cazy`: most represented consensus CAZy families.
4. `figure_03_top_go`: most represented GO terms.
5. `figure_04_top_kegg_ko`: most represented KEGG Orthology identifiers.
6. `figure_05_protein_length_and_short_proteins`: length densities and proteins shorter than 100 aa.
7. `figure_06_secretome_effector_comparison`: SignalP-positive, final secretome and effector coverage.
8. `figure_07_functional_coverage_heatmap`: annotation presence and coverage.
9. `figure_08_secretome_effector_flow`: sequential SignalP–DeepTMHMM–EffectorP filtering.
10. `figure_09_effectorp_composition`: apoplastic-primary, cytoplasmic-primary and non-effector proteins, with the overlapping dual subset reported separately.
11. `figure_10_dbcan_evidence_support`: one-, two- and three-tool dbCAN support.

`plot_data/` contains the summarized TSV data behind the figures and `figure_manifest.tsv`.
It also contains `run_metadata.tsv`, including the MD5 checksum of the exact master table used. `R_sessionInfo.txt` records the R and package versions.

## Statistical and biological interpretation

- Pfam, GO, KEGG and CAZy panels count unique proteins per term. A protein carrying several terms contributes once to each relevant term.
- The top terms are selected using the combined count across the five proteomes, after which the same term set is shown for every proteome. This avoids choosing a different and incomparable top list for each sample.
- CAZy family plots use the consensus field, which contains calls supported by at least two dbCAN tools.
- The coverage heatmap reports both the count and percentage of each complete proteome.
- The secretome-to-effector flow includes only nested sets: all proteins, SignalP-positive proteins, DeepTMHMM-filtered secretome candidates and EffectorP candidates.
- EffectorP primary classes are mutually exclusive. `dual_localized` is an overlapping property and therefore appears as an annotation rather than an additional stacked segment.
- Proteins shorter than 100 aa remain in every relevant analysis and are displayed separately; they are not removed.
- These figures are descriptive. They do not by themselves test enrichment or differential representation. Formal comparisons require an explicit statistical model and correction for multiple testing.

## Reproducibility checks

After completion:

```bash
find "$RESULTS/08_figures/plots" -type f -size +0c | sort
```

There should be 33 non-empty plot files: 11 figures multiplied by three formats.

```bash
find "$RESULTS/08_figures/plots" -type f \( -name '*.pdf' -o -name '*.svg' -o -name '*.png' \) | wc -l
```

The expected result is `33`.
