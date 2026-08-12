#!/usr/bin/env bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -N pggb_func_plots
#$ -j y

set -euo pipefail

# SGE runs a spool copy of this file, so $0 must not be used to locate the
# project scripts. Submit from PROJECT_ROOT or export PLOTTING_ENV explicitly.
SUBMISSION_DIR="${SGE_O_WORKDIR:-$(pwd)}"
ENV_FILE="${PLOTTING_ENV:-$SUBMISSION_DIR/plotting_environment_PGGB.env}"
[[ -r "$ENV_FILE" ]] || { echo "ERROR: cannot read $ENV_FILE" >&2; exit 1; }

set -a
source "$ENV_FILE"
set +a

SCRIPT_PATH="$PROJECT_ROOT/09_plot_functional_annotation_PGGB.R"
[[ -r "$SCRIPT_PATH" ]] || { echo "ERROR: cannot read $SCRIPT_PATH" >&2; exit 1; }

if [[ -n "${R_MODULE:-}" ]]; then
    module load "$R_MODULE"
fi

export OMP_NUM_THREADS="${THREADS:-1}"
export OPENBLAS_NUM_THREADS="${THREADS:-1}"
export MKL_NUM_THREADS="${THREADS:-1}"

command=(
    Rscript "$SCRIPT_PATH"
    --master "$FINAL_MASTER"
    --output-dir "$FIGURE_ROOT"
    --top-n "$TOP_N"
    --short-threshold "$SHORT_PROTEIN_THRESHOLD"
    --png-dpi "$PNG_DPI"
)
if [[ -n "${GO_OBO:-}" ]]; then
    command+=(--go-obo "$GO_OBO")
fi
"${command[@]}"
