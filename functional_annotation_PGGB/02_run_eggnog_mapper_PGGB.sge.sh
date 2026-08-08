#!/bin/bash
#$ -N eggnog_C26
#$ -q all.q
#$ -cwd
#$ -pe smp 16
#$ -V
#$ -o /Storage/data1/isabella.gallego/MAESTRIA/code/functional_annotation/logs/eggnog_$JOB_ID.out
#$ -e /Storage/data1/isabella.gallego/MAESTRIA/code/functional_annotation/logs/eggnog_$JOB_ID.err
set -euo pipefail
BASE="/Storage/data1/isabella.gallego/MAESTRIA/functional_annotation_C26"
ENV="${BASE}/software/envs/eggnog_mapper_2.1.9"
DB="${BASE}/databases/eggnog_2.1.9"
INPUT="${BASE}/00_qc/MrorC26_proteins.cleaned.fa"
OUT="${BASE}/02_eggnog_mapper"
mkdir -p "$OUT"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV"
emapper.py -i "$INPUT" --itype proteins -m diamond \
  --data_dir "$DB" --tax_scope Fungi --target_orthologs all \
  --go_evidence all --pfam_realign none \
  --cpu "${NSLOTS:-16}" --output MrorC26 --output_dir "$OUT" --override
