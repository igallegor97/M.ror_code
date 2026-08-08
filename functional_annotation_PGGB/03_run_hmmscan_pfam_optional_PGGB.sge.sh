#!/bin/bash
#$ -N hmmscan_C26
#$ -q all.q
#$ -cwd
#$ -pe smp 8
#$ -V
#$ -o /Storage/data1/isabella.gallego/MAESTRIA/code/functional_annotation/logs/hmmscan_$JOB_ID.out
#$ -e /Storage/data1/isabella.gallego/MAESTRIA/code/functional_annotation/logs/hmmscan_$JOB_ID.err
set -euo pipefail
module purge
module load HMMER/3.4 2>/dev/null || module load hmmer/3.4
BASE="/Storage/data1/isabella.gallego/MAESTRIA/functional_annotation_C26"
INPUT="${BASE}/00_qc/MrorC26_proteins.cleaned.fa"
PFAM="${BASE}/databases/pfam/Pfam-A.hmm"
OUT="${BASE}/03_hmmscan_pfam_optional"
mkdir -p "$OUT"
[[ -f "${PFAM}.h3f" ]] || hmmpress "$PFAM"
hmmscan --cpu "${NSLOTS:-8}" --cut_ga --noali \
  --domtblout "${OUT}/MrorC26_Pfam.domtblout" "$PFAM" "$INPUT" \
  > "${OUT}/MrorC26_Pfam.hmmscan.txt"
