#!/bin/bash

#$ -N pggb24_geo_gen
#$ -q all.q
#$ -cwd
#$ -pe smp 4
#$ -V
#$ -o /Storage/data1/isabella.gallego/MAESTRIA/code/population_structure/geography_genetics/logs/pggb24_geo_gen_$JOB_ID.out
#$ -e /Storage/data1/isabella.gallego/MAESTRIA/code/population_structure/geography_genetics/logs/pggb24_geo_gen_$JOB_ID.err

# =============================================================
# PGGB24 geography-versus-genetics pipeline — SGE job
# =============================================================

set -euo pipefail

CODE_DIR="/Storage/data1/isabella.gallego/MAESTRIA/code/population_structure/geography_genetics"
PIPELINE="${CODE_DIR}/run_geography_genetics_pipeline.sh"
LOG_DIR="/Storage/data1/isabella.gallego/MAESTRIA/code/population_structure/geography_genetics/logs"

mkdir -p "${LOG_DIR}"

module purge
module load Python/3.13.5
module load R/4.5.1

export OMP_NUM_THREADS="${NSLOTS:-4}"
export OPENBLAS_NUM_THREADS="${NSLOTS:-4}"
export MKL_NUM_THREADS="${NSLOTS:-4}"

export REPLICATES="${REPLICATES:-100}"
export BASE_SEED="${BASE_SEED:-20260803}"
export N_PER_COMMUNITY="${N_PER_COMMUNITY:-0}"

[[ -f "${PIPELINE}" ]] || {
    echo "[ERROR] Pipeline not found: ${PIPELINE}"
    exit 1
}

python -m py_compile "${CODE_DIR}/validate_coordinates.py"
python -c "import pandas" 2>/dev/null || {
    echo "[ERROR] Python package pandas is required."
    exit 1
}

Rscript -e '
required <- c(
  "geosphere","vegan","ecodist","ggplot2","ggrepel",
  "sf","rnaturalearth","rnaturalearthdata","ggspatial"
)
missing <- required[
  !vapply(required,requireNamespace,logical(1),quietly=TRUE)
]
if (length(missing)>0) {
  stop(paste("Missing R packages:",paste(missing,collapse=", ")))
}
'

echo "============================================================="
echo "PGGB24 GEOGRAPHY-GENETICS — SGE JOB"
echo "Start      : $(date '+%Y-%m-%d %H:%M:%S')"
echo "Job ID     : ${JOB_ID:-not_available}"
echo "Node       : $(hostname)"
echo "Threads    : ${NSLOTS:-4}"
echo "Replicates : ${REPLICATES}"
echo "============================================================="

cd "${CODE_DIR}"
bash "${PIPELINE}"
