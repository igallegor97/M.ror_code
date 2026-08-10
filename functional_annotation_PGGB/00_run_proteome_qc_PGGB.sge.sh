#!/usr/bin/env bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -j y
#$ -N pggb_func_qc
#$ -pe smp 1

# =============================================================
# PacBio proteome validation and normalization (SGE array)
#
# Usage:
#   qsub -t 1-5 00_run_proteome_qc_PGGB.sge.sh
# =============================================================

set -euo pipefail

: "${PROJECT_ROOT:?PROJECT_ROOT was not exported by the submission driver}"
source "$PROJECT_ROOT/pipeline_config_PGGB.sh"
initialize_array_sample

mkdir -p "$RESULTS/00_qc/fasta" "$RESULTS/00_qc/id_maps" "$RESULTS/00_qc/stats"

echo "============================================================="
echo "PROTEOME QC"
echo "Sample       : $SAMPLE"
echo "Input        : $INPUT_FASTA"
echo "============================================================="

python3 "$PROJECT_ROOT/00_validate_proteomes_PGGB.py" \
    "$INPUT_FASTA" \
    "$QC_FASTA" \
    "$RESULTS/00_qc/id_maps/${SAMPLE}.id_map.tsv" \
    "$RESULTS/00_qc/stats/${SAMPLE}.qc.tsv"

echo "PROTEOME QC COMPLETED: $SAMPLE"

