#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_FILE="${1:-$SCRIPT_DIR/advanced_viz_environment_PGGB.env}"
[[ -r "$ENV_FILE" ]] || { echo "ERROR: cannot read $ENV_FILE" >&2; exit 1; }
set -a; source "$ENV_FILE"; set +a

[[ -s "$FINAL_MASTER" ]] || { echo "ERROR: missing final master: $FINAL_MASTER" >&2; exit 1; }
mkdir -p "$ADVANCED_VIZ_ROOT" "$LOGS"
if [[ -n "${R_MODULE:-}" ]]; then module load "$R_MODULE"; fi
command -v Rscript >/dev/null 2>&1 || { echo "ERROR: Rscript unavailable; set R_MODULE" >&2; exit 1; }

Rscript -e 'p <- c("data.table","ggplot2","scales","patchwork","svglite","ComplexUpset","ggalluvial","ggrepel","ggdist"); m <- p[!vapply(p,requireNamespace,logical(1),quietly=TRUE)]; if(length(m)) stop("Missing R packages: ",paste(m,collapse=", "))'

header="$(head -n 1 "$FINAL_MASTER")"
for column in sample_id protein_id protein_length GO_terms KEGG_ko pfam_domains cazy_families signalp_positive secretome_candidate effector_candidate effectorp_primary_class effectorp_dual_localized effectorp_apoplastic_probability effectorp_cytoplasmic_probability; do
  printf '%s\n' "$header" | tr '\t' '\n' | grep -Fxq "$column" || { echo "ERROR: missing column: $column" >&2; exit 1; }
done

echo "ADVANCED VIZ PREFLIGHT COMPLETED"
echo "Master : $FINAL_MASTER"
echo "Output : $ADVANCED_VIZ_ROOT"
echo "Rscript: $(command -v Rscript)"

