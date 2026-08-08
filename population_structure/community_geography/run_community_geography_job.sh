#!/bin/bash

#$ -N comm_geo_24G
#$ -q all.q
#$ -cwd
#$ -pe smp 4
#$ -V
#$ -o /Storage/data1/isabella.gallego/MAESTRIA/code/population_structure/community_geography/logs/comm_geo_24G_$JOB_ID.out
#$ -e /Storage/data1/isabella.gallego/MAESTRIA/code/population_structure/community_geography/logs/comm_geo_24G_$JOB_ID.err

set -euo pipefail

CODE_DIR="/Storage/data1/isabella.gallego/MAESTRIA/code/population_structure/community_geography"

module purge
module load Python/3.13.5
module load R/4.5.1

export OMP_NUM_THREADS="${NSLOTS:-4}"
export OPENBLAS_NUM_THREADS="${NSLOTS:-4}"
export MKL_NUM_THREADS="${NSLOTS:-4}"

python -m py_compile \
    "${CODE_DIR}/rank_communities.py" \
    "${CODE_DIR}/validate_community_chromosome_map.py" \
    "${CODE_DIR}/annotate_community_snps.py" \
    "${CODE_DIR}/enrich_top_community_genes.py"

python -c "import pandas, scipy, numpy" 2>/dev/null || {
    echo "[ERROR] Python packages pandas, scipy and numpy are required."
    exit 1
}

Rscript -e '
required <- c("vegan", "ecodist", "ggplot2", "scales", "gridExtra")
missing <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing) > 0) {
  stop(paste("Missing R packages:", paste(missing, collapse=", ")))
}
'

cd "${CODE_DIR}"

bash "${CODE_DIR}/run_community_geography_pipeline.sh"
