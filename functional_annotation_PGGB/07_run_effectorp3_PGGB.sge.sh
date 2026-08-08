#!/bin/bash
#$ -N effectorp_C26
#$ -q all.q
#$ -cwd
#$ -pe smp 2
#$ -V
#$ -o /Storage/data1/isabella.gallego/MAESTRIA/code/functional_annotation/logs/effectorp_$JOB_ID.out
#$ -e /Storage/data1/isabella.gallego/MAESTRIA/code/functional_annotation/logs/effectorp_$JOB_ID.err
set -euo pipefail
module purge
module load Python/3.12.3 2>/dev/null || true
module load Java/11 2>/dev/null || true
BASE="/Storage/data1/isabella.gallego/MAESTRIA/functional_annotation_C26"
SOFT="${BASE}/software/EffectorP-3.0"
INPUT="${BASE}/06_secretome/MrorC26_secreted_candidates.fa"
OUT="${BASE}/07_effectorp3"
mkdir -p "$OUT"
cd "$SOFT"
python EffectorP.py -i "$INPUT" > "${OUT}/MrorC26_EffectorP3.txt"
