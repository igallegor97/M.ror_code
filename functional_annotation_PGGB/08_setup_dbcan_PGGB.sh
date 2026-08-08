#!/bin/bash
set -euo pipefail
BASE="/Storage/data1/isabella.gallego/MAESTRIA/functional_annotation_C26"
ENV="${BASE}/software/envs/run_dbcan"
DB="${BASE}/databases/dbcan_v5_2_pinned"
mkdir -p "$(dirname "$ENV")" "$DB"
conda create -y -p "$ENV" -c conda-forge -c bioconda python=3.11 dbcan diamond hmmer
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV"
run_dbcan database --db_dir "$DB" --aws_s3 --no-cgc --timeout 60 --retries 5 --resume
