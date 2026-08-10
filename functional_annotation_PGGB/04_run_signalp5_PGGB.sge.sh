#!/usr/bin/env bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -j y
#$ -N pggb_func_signalp
#$ -pe smp 4

# =============================================================
# Eukaryotic signal-peptide prediction with SignalP 5.0b
# =============================================================

set -euo pipefail

: "${PROJECT_ROOT:?PROJECT_ROOT was not exported by the submission driver}"
source "$PROJECT_ROOT/pipeline_config_PGGB.sh"
initialize_array_sample
load_modules_init
module purge
module load "$SIGNALP_MODULE"

OUTPUT_DIR="$RESULTS/04_signalp/$SAMPLE"
TEMP_DIR="$OUTPUT_DIR/tmp"
mkdir -p "$OUTPUT_DIR" "$TEMP_DIR"

echo "============================================================="
echo "SIGNALP 5"
echo "Sample       : $SAMPLE"
echo "Input        : $QC_FASTA"
echo "============================================================="

signalp \
    -fasta "$QC_FASTA" \
    -org euk \
    -format short \
    -mature \
    -gff3 \
    -prefix "$OUTPUT_DIR/$SAMPLE" \
    -tmp "$TEMP_DIR" \
    -batch 1000 \
    -verbose=false

echo "SIGNALP 5 COMPLETED: $SAMPLE"

