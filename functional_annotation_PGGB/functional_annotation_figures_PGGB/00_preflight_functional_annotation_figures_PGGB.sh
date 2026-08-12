#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_FILE="${1:-$SCRIPT_DIR/plotting_environment_PGGB.env}"

[[ -r "$ENV_FILE" ]] || { echo "ERROR: cannot read $ENV_FILE" >&2; exit 1; }
set -a
source "$ENV_FILE"
set +a

[[ -s "$FINAL_MASTER" ]] || { echo "ERROR: final master not found or empty: $FINAL_MASTER" >&2; exit 1; }
mkdir -p "$FIGURE_ROOT" "$LOGS"

if [[ -n "${R_MODULE:-}" ]]; then
    module load "$R_MODULE"
fi
command -v Rscript >/dev/null 2>&1 || { echo "ERROR: Rscript is unavailable. Set R_MODULE." >&2; exit 1; }

Rscript -e 'required <- c("data.table","ggplot2","scales","patchwork","svglite"); missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]; if(length(missing)) stop("Missing R packages: ", paste(missing, collapse=", "))'

header="$(head -n 1 "$FINAL_MASTER")"
for column in sample_id protein_id protein_length short_protein_flag GO_terms KEGG_ko pfam_domains cazy_families signalp_positive secretome_candidate effector_candidate; do
    printf '%s\n' "$header" | tr '\t' '\n' | grep -Fxq "$column" || {
        echo "ERROR: missing column in final master: $column" >&2
        exit 1
    }
done

echo "FIGURE PREFLIGHT COMPLETED"
echo "Final master : $FINAL_MASTER"
echo "Figure root : $FIGURE_ROOT"
echo "Rscript     : $(command -v Rscript)"

