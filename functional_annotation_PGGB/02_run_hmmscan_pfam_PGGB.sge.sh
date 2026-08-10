#!/usr/bin/env bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -j y
#$ -N pggb_func_pfam
#$ -pe smp 8

# =============================================================
# Pfam-A v35 domain annotation with HMMER 3.2.1
# =============================================================

set -euo pipefail

: "${PROJECT_ROOT:?PROJECT_ROOT was not exported by the submission driver}"
source "$PROJECT_ROOT/pipeline_config_PGGB.sh"
initialize_array_sample
load_modules_init
module purge
module load "$HMMER_MODULE"

OUTPUT_DIR="$RESULTS/02_pfam/$SAMPLE"
mkdir -p "$OUTPUT_DIR"

echo "============================================================="
echo "PFAM HMMSCAN"
echo "Sample       : $SAMPLE"
echo "Input        : $QC_FASTA"
echo "Database     : $PFAM_HMM"
echo "============================================================="

hmmscan \
    --cpu "${NSLOTS:-$THREADS}" \
    --cut_ga \
    --noali \
    --domtblout "$OUTPUT_DIR/${SAMPLE}.pfam.domtblout" \
    -o "$OUTPUT_DIR/${SAMPLE}.pfam.txt" \
    "$PFAM_HMM" \
    "$QC_FASTA"

echo "PFAM HMMSCAN COMPLETED: $SAMPLE"

