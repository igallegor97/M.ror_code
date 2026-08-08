#!/usr/bin/env bash
#$ -cwd
#$ -V
#$ -j y
#$ -N PB_PFAM
#$ -pe smp 8
set -euo pipefail
source pipeline_config_PGGB.sh; init_sample; load_modules_init; module purge; module load "$HMMER_MODULE"
OUT="$RESULTS/02_pfam/$SAMPLE"; mkdir -p "$OUT"
hmmscan --cpu "${NSLOTS:-$THREADS}" --cut_ga --noali --domtblout "$OUT/${SAMPLE}.pfam.domtblout" -o "$OUT/${SAMPLE}.pfam.txt" "$PFAM_HMM" "$QC_FASTA"
