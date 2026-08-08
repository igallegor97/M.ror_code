#!/usr/bin/env bash
#$ -cwd
#$ -V
#$ -j y
#$ -N PB_SIGNALP
#$ -pe smp 4
set -euo pipefail
source pipeline_config_PGGB.sh; init_sample; load_modules_init; module purge; module load "$SIGNALP_MODULE"
OUT="$RESULTS/04_signalp/$SAMPLE"; TMP="$OUT/tmp"; mkdir -p "$OUT" "$TMP"
signalp -fasta "$QC_FASTA" -org euk -format short -mature -gff3 -prefix "$OUT/$SAMPLE" -tmp "$TMP" -batch 1000 -verbose=false
