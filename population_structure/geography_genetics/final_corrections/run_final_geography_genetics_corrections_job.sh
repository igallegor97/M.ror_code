#!/bin/bash

#$ -N pggb24_geo_final
#$ -q all.q
#$ -cwd
#$ -pe smp 4
#$ -V
#$ -o /Storage/data1/isabella.gallego/MAESTRIA/code/population_structure/geography_genetics/final_corrections/logs/pggb24_geo_final_$JOB_ID.out
#$ -e /Storage/data1/isabella.gallego/MAESTRIA/code/population_structure/geography_genetics/final_corrections/logs/pggb24_geo_final_$JOB_ID.err

# =============================================================================
# Final PGGB24 geography-genetics corrections - SGE job
# =============================================================================

set -euo pipefail

CODE_DIR="/Storage/data1/isabella.gallego/MAESTRIA/code/population_structure/geography_genetics/final_corrections"

PIPELINE="${CODE_DIR}/run_final_geography_genetics_corrections.sh"

LOG_DIR="/Storage/data1/isabella.gallego/MAESTRIA/code/population_structure/geography_genetics/final_corrections/logs"

mkdir -p "${LOG_DIR}"

module purge
module load Python/3.13.5
module load R/4.5.1

export OMP_NUM_THREADS="${NSLOTS:-4}"
export OPENBLAS_NUM_THREADS="${NSLOTS:-4}"
export MKL_NUM_THREADS="${NSLOTS:-4}"

python -m py_compile \
    "${CODE_DIR}/audit_community_profiles.py" \
    "${CODE_DIR}/build_pipeline_master_summary.py"

python -c "import pandas" 2>/dev/null || {
    echo "[ERROR] Python package pandas is required."
    exit 1
}

Rscript -e '
required <- c(
  "geosphere",
  "vegan",
  "ecodist",
  "ggplot2",
  "ggrepel",
  "sf",
  "rnaturalearth",
  "rnaturalearthdata",
  "ggspatial"
)

missing <- required[
  !vapply(
    required,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing) > 0) {
  stop(
    paste(
      "Missing R packages:",
      paste(missing, collapse = ", ")
    )
  )
}
'

echo "======================================================================"
echo "FINAL PGGB24 GEOGRAPHY-GENETICS - SGE JOB"
echo "Start      : $(date '+%Y-%m-%d %H:%M:%S')"
echo "Job ID     : ${JOB_ID:-not_available}"
echo "Node       : $(hostname)"
echo "Threads    : ${NSLOTS:-4}"
echo "======================================================================"

cd "${CODE_DIR}"

bash "${PIPELINE}"
