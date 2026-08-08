#!/usr/bin/env bash
#$ -cwd
#$ -V
#$ -j y
#$ -N PB_EGGNOG
#$ -pe smp 8
set -euo pipefail
source pipeline_config_PGGB.sh; init_sample; load_modules_init; module purge; module load "$EGGNOG_MODULE"
OUT="$RESULTS/01_eggnog/$SAMPLE"; mkdir -p "$OUT"
emapper.py -i "$QC_FASTA" --itype proteins -m diamond --data_dir "$EGGNOG_DATA" --cpu "${NSLOTS:-$THREADS}" --output "$SAMPLE" --output_dir "$OUT" --resume
