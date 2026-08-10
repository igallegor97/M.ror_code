#!/usr/bin/env bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -j y
#$ -N pggb_func_eggnog
#$ -pe smp 8

# =============================================================
# eggNOG-mapper 2.1.13 annotation (SGE array)
#
# Resumes only when a valid hits file exists. New or incomplete
# runs start a fresh DIAMOND search with --override.
# =============================================================

set -euo pipefail

: "${PROJECT_ROOT:?PROJECT_ROOT was not exported by the submission driver}"
source "$PROJECT_ROOT/pipeline_config_PGGB.sh"
initialize_array_sample
load_modules_init
module purge
module load "$EGGNOG_MODULE"

OUTPUT_DIR="$RESULTS/01_eggnog/$SAMPLE"
mkdir -p "$OUTPUT_DIR"

if [[ -s "$OUTPUT_DIR/${SAMPLE}.emapper.hits" ]]; then
    RUN_MODE=(--resume)
    RUN_MODE_LABEL="resume from existing hits"
else
    RUN_MODE=(--override)
    RUN_MODE_LABEL="fresh search"
fi

echo "============================================================="
echo "EGGNOG-MAPPER"
echo "Sample       : $SAMPLE"
echo "Input        : $QC_FASTA"
echo "Database     : $EGGNOG_DATA"
echo "Threads      : ${NSLOTS:-$THREADS}"
echo "Run mode     : $RUN_MODE_LABEL"
echo "============================================================="

emapper.py \
    -i "$QC_FASTA" \
    --itype proteins \
    -m diamond \
    --data_dir "$EGGNOG_DATA" \
    --cpu "${NSLOTS:-$THREADS}" \
    --output "$SAMPLE" \
    --output_dir "$OUTPUT_DIR" \
    "${RUN_MODE[@]}"

echo "EGGNOG-MAPPER COMPLETED: $SAMPLE"
