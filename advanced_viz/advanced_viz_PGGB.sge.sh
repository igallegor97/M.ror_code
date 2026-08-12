#!/usr/bin/env bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -N pggb_adv_viz
#$ -j y

set -euo pipefail

SUBMISSION_DIR="${SGE_O_WORKDIR:-$(pwd)}"
ENV_FILE="${ADVANCED_VIZ_ENV:-$SUBMISSION_DIR/advanced_viz_environment_PGGB.env}"
[[ -r "$ENV_FILE" ]] || { echo "ERROR: cannot read $ENV_FILE" >&2; exit 1; }
set -a; source "$ENV_FILE"; set +a

SCRIPT_PATH="$PROJECT_ROOT/advanced_viz_PGGB.R"
[[ -r "$SCRIPT_PATH" ]] || { echo "ERROR: cannot read $SCRIPT_PATH" >&2; exit 1; }
if [[ -n "${R_MODULE:-}" ]]; then module load "$R_MODULE"; fi
export OMP_NUM_THREADS="${THREADS:-1}"
export OPENBLAS_NUM_THREADS="${THREADS:-1}"
export MKL_NUM_THREADS="${THREADS:-1}"

command=(Rscript "$SCRIPT_PATH" --master "$FINAL_MASTER" --output-dir "$ADVANCED_VIZ_ROOT" --top-n "$TOP_N" --png-dpi "$PNG_DPI")
if [[ -n "${GO_OBO:-}" ]]; then command+=(--go-obo "$GO_OBO"); fi
"${command[@]}"

