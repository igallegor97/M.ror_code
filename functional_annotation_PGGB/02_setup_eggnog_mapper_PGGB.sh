#!/bin/bash
set -euo pipefail
BASE="/Storage/data1/isabella.gallego/MAESTRIA/functional_annotation_C26"
ENV="${BASE}/software/envs/eggnog_mapper_2.1.9"
DB="${BASE}/databases/eggnog_2.1.9"
mkdir -p "$(dirname "$ENV")" "$DB"
conda create -y -p "$ENV" -c conda-forge -c bioconda python=3.10 eggnog-mapper=2.1.9 diamond
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV"
download_eggnog_data.py --data_dir "$DB" -y
emapper.py --version
