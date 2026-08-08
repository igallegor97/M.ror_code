#!/bin/bash

#$ -N pggb24_pop_v2
#$ -q all.q
#$ -cwd
#$ -pe smp 8
#$ -V
#$ -o /Storage/data1/isabella.gallego/MAESTRIA/code/population_structure/population_structure_24G_v2/logs/pggb24_pop_v2_$JOB_ID.out
#$ -e /Storage/data1/isabella.gallego/MAESTRIA/code/population_structure/population_structure_24G_v2/logs/pggb24_pop_v2_$JOB_ID.err

# =============================================================
# PGGB 24-genome population-structure pipeline v2 — SGE job
# =============================================================

set -euo pipefail

CODE_DIR="/Storage/data1/isabella.gallego/MAESTRIA/code/population_structure/population_structure_24G_v2"
PIPELINE="${CODE_DIR}/run_24G_population_pipeline.sh"
LOG_DIR="/Storage/data1/isabella.gallego/MAESTRIA/code/PGGB_24_genomes/logs"
THREADS="${NSLOTS:-8}"

mkdir -p "${LOG_DIR}"

module purge
module load Python/3.13.5
module load R/4.5.1

export PYTHON="python"
export RSCRIPT="Rscript"
export MAX_MISSING="${MAX_MISSING:-0.20}"

export OMP_NUM_THREADS="${THREADS}"
export OPENBLAS_NUM_THREADS="${THREADS}"
export MKL_NUM_THREADS="${THREADS}"
export NUMEXPR_NUM_THREADS="${THREADS}"

[[ -f "${PIPELINE}" ]] || { echo "[ERROR] Pipeline not found: ${PIPELINE}"; exit 1; }

python -m py_compile "${CODE_DIR}/build_24G_snp_matrix.py"
python -c "import numpy, pandas" 2>/dev/null || {
    echo "[ERROR] Python packages numpy and pandas are required."
    exit 1
}

Rscript -e '
required <- c("ggplot2","ggrepel")
missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if (length(missing)>0) stop(paste("Missing R packages:",paste(missing,collapse=", ")))
'

echo "============================================================="
echo "PGGB 24G POPULATION STRUCTURE V2 — SGE JOB"
echo "Start       : $(date '+%Y-%m-%d %H:%M:%S')"
echo "Job ID      : ${JOB_ID:-not_available}"
echo "Node        : $(hostname)"
echo "Threads     : ${THREADS}"
echo "MAX_MISSING : ${MAX_MISSING}"
echo "Python      : $(python --version 2>&1)"
echo "R           : $(Rscript --version 2>&1 | head -1)"
echo "============================================================="

cd "${CODE_DIR}"
bash "${PIPELINE}"
