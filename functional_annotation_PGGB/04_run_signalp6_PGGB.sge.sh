#!/bin/bash
#$ -N signalp6_C26
#$ -q all.q
#$ -cwd
#$ -pe smp 8
#$ -V
#$ -o /Storage/data1/isabella.gallego/MAESTRIA/code/functional_annotation/logs/signalp6_$JOB_ID.out
#$ -e /Storage/data1/isabella.gallego/MAESTRIA/code/functional_annotation/logs/signalp6_$JOB_ID.err
set -euo pipefail
BASE="/Storage/data1/isabella.gallego/MAESTRIA/functional_annotation_C26"
ENV="${BASE}/software/envs/signalp6"
INPUT="${BASE}/00_qc/MrorC26_proteins.cleaned.fa"
OUT="${BASE}/04_signalp6"
rm -rf "$OUT"; mkdir -p "$OUT"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV"
signalp6 --fastafile "$INPUT" --organism eukarya --output_dir "$OUT" \
  --format none --mode fast --torch_num_threads "${NSLOTS:-8}" --write_procs 1
