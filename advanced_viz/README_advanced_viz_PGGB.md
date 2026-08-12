# Advanced functional-annotation visualizations

`advanced_viz_PGGB.R` creates eight complementary visualizations from `functional_annotation_master_final.tsv`:

1. annotation-layer UpSet plot;
2. Pfam, CAZy, GO and KEGG KO dot plots;
3. clustered functional-profile heatmap;
4. PCA of Hellinger-transformed functional profiles;
5. SignalP–DeepTMHMM–EffectorP alluvial flow;
6. protein-length raincloud distributions;
7. EffectorP probability plane;
8. hierarchical CAZy sunburst.

Every main figure is exported as PDF, SVG and 300-dpi PNG. The summarized data behind each figure, an input MD5 checksum and the R session information are retained.

## Installation and configuration

Extract the package files directly into:

```text
/Storage/data1/isabella.gallego/MAESTRIA/code/functional_annotation_PGGB
```

Edit `advanced_viz_environment_PGGB.env`, especially `R_MODULE`. The required packages are:

```text
data.table ggplot2 scales patchwork svglite ComplexUpset ggalluvial ggrepel ggdist
```

If user-level CRAN installation is allowed:

```bash
module load YOUR_R_MODULE
Rscript 00_install_advanced_viz_R_packages_PGGB.R
```

Do not install packages repeatedly inside an SGE job.

## Preflight and submission

```bash
cd /Storage/data1/isabella.gallego/MAESTRIA/code/functional_annotation_PGGB
bash 00_preflight_advanced_viz_PGGB.sh
```

After `ADVANCED VIZ PREFLIGHT COMPLETED`:

```bash
source ./advanced_viz_environment_PGGB.env

qsub -terse \
  -pe smp 4 \
  -l h_vmem=32G \
  -l h_rt=12:00:00 \
  -o "$LOGS" \
  -e "$LOGS" \
  "$PROJECT_ROOT/advanced_viz_PGGB.sge.sh"
```

Outputs are written to:

```text
$RESULTS/09_advanced_viz/plots
$RESULTS/09_advanced_viz/plot_data
```

## Interpretation safeguards

- The UpSet plot describes intersections, not statistical enrichment.
- Dot size is an absolute protein count; dot colour is proteome coverage.
- Heatmap values are row z-scores and must not be interpreted as absolute abundance.
- PCA uses Hellinger-transformed term profiles and is descriptive because there are only five proteomes.
- The alluvial plot begins with SignalP-positive proteins; it does not imply that all proteins were evaluated by DeepTMHMM or EffectorP.
- Raincloud axes stop at the 99th percentile for visibility; all values remain in the exported TSV.
- EffectorP secondary probabilities absent from downloaded FASTA headers are censored at `≤0.5`. They are plotted on the 0.5 boundary and are not treated as exact probabilities.
- The CAZy sunburst uses consensus families and groups families outside the top ten of each class as `Other`.

