#!/bin/bash
#$ -N dbcan_C26
#$ -q all.q
#$ -cwd
#$ -pe smp 16
#$ -V
#$ -o /Storage/data1/isabella.gallego/MAESTRIA/code/functional_annotation/logs/dbcan_$JOB_ID.out
#$ -e /Storage/data1/isabella.gallego/MAESTRIA/code/functional_annotation/logs/dbcan_$JOB_ID.err
set -euo pipefail
BASE="/Storage/data1/isabella.gallego/MAESTRIA/functional_annotation_C26"
ENV="${BASE}/software/envs/run_dbcan"
DB="${BASE}/databases/dbcan_v5_2_pinned"
INPUT="${BASE}/00_qc/MrorC26_proteins.cleaned.fa"
OUT="${BASE}/08_dbcan"
mkdir -p "$OUT"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV"
run_dbcan CAZyme_annotation --input_raw_data "$INPUT" --mode protein \
  --output_dir "$OUT" --db_dir "$DB" --threads "${NSLOTS:-16}"
