#!/usr/bin/env bash
#$ -cwd
#$ -V
#$ -j y
#$ -N PB_QC
#$ -pe smp 1
set -euo pipefail
source pipeline_config_PGGB.sh; init_sample
mkdir -p "$RESULTS/00_qc/fasta" "$RESULTS/00_qc/id_maps" "$RESULTS/00_qc/stats"
python3 00_validate_proteomes_PGGB.py "$INPUT_FASTA" "$QC_FASTA" "$RESULTS/00_qc/id_maps/${SAMPLE}.id_map.tsv" "$RESULTS/00_qc/stats/${SAMPLE}.qc.tsv"
