#!/bin/bash
#$ -N iprscan6_C26
#$ -q all.q
#$ -cwd
#$ -pe smp 16
#$ -V
#$ -o /Storage/data1/isabella.gallego/MAESTRIA/code/functional_annotation/logs/iprscan6_$JOB_ID.out
#$ -e /Storage/data1/isabella.gallego/MAESTRIA/code/functional_annotation/logs/iprscan6_$JOB_ID.err
set -euo pipefail
module purge
module load singularity-ce/3.11.2
module load Nextflow/25.10.4 2>/dev/null || true
BASE="/Storage/data1/isabella.gallego/MAESTRIA/functional_annotation_C26"
INPUT="${BASE}/00_qc/MrorC26_proteins.cleaned.fa"
OUT="${BASE}/01_interproscan6"
DATA="${BASE}/databases/interproscan6"
WORK="${BASE}/work/interproscan6"
export NXF_HOME="${BASE}/software/nextflow"
mkdir -p "$OUT" "$DATA" "$WORK" "$NXF_HOME"
command -v nextflow >/dev/null || { echo "[ERROR] Nextflow 25.10.4+ required"; exit 1; }
nextflow run ebi-pf-team/interproscan6 \
  -r "${INTERPROSCAN_VERSION:-6.0.1}" \
  -profile singularity \
  -w "$WORK" \
  --input "$INPUT" \
  --datadir "$DATA" \
  --interpro "${INTERPRO_RELEASE:-latest}" \
  --outdir "$OUT" \
  --outprefix "MrorC26_interproscan6" \
  --formats "tsv,json,gff3" \
  --goterms --pathways \
  --cpus "${NSLOTS:-16}" \
  --max-workers 1 \
  -resume
