#!/bin/bash

#$ -N balanced_pca
#$ -q all.q
#$ -cwd
#$ -pe smp 4
#$ -V
#$ -o logs/balanced_pca_$JOB_ID.out
#$ -e logs/balanced_pca_$JOB_ID.err

# =============================================================
# Partition-balanced PCA — SGE job
# =============================================================

set -euo pipefail

CODE_DIR="/Storage/data1/isabella.gallego/MAESTRIA/code/population_structure"
PIPELINE_SCRIPT="${CODE_DIR}/run_balanced_pca_pipeline.sh"
LOG_DIR="${CODE_DIR}/logs"
THREADS="${NSLOTS:-4}"

export TOTAL_SNPS="${TOTAL_SNPS:-1000}"
export RANDOM_SEED="${RANDOM_SEED:-42}"

mkdir -p "${LOG_DIR}"

start_time="$(date '+%Y-%m-%d %H:%M:%S')"

echo "============================================================="
echo "PARTITION-BALANCED PCA — SGE JOB"
echo "Start       : ${start_time}"
echo "Job ID      : ${JOB_ID:-not_available}"
echo "Node        : $(hostname)"
echo "CODE_DIR    : ${CODE_DIR}"
echo "TOTAL_SNPS  : ${TOTAL_SNPS}"
echo "RANDOM_SEED : ${RANDOM_SEED}"
echo "THREADS     : ${THREADS}"
echo "============================================================="
echo ""

module purge
module load Python/3.13.5
module load R/4.5.1

export PYTHON="python"
export RSCRIPT="Rscript"
export OMP_NUM_THREADS="${THREADS}"
export OPENBLAS_NUM_THREADS="${THREADS}"
export MKL_NUM_THREADS="${THREADS}"
export NUMEXPR_NUM_THREADS="${THREADS}"

echo "Python : $(python --version 2>&1)"
echo "R      : $(Rscript --version 2>&1 | head -1)"
echo ""

if [[ ! -f "${PIPELINE_SCRIPT}" ]]; then
    echo "[ERROR] Pipeline script not found: ${PIPELINE_SCRIPT}"
    exit 1
fi

for required_file in \
    "${CODE_DIR}/sample_balanced_snps.py" \
    "${CODE_DIR}/analyze_population_structure.R" \
    "${CODE_DIR}/compare_pggb_cactus_pca.R" \
    "${CODE_DIR}/sample_metadata.tsv"
do
    if [[ ! -f "${required_file}" ]]; then
        echo "[ERROR] Required file not found: ${required_file}"
        exit 1
    fi
done

python -m py_compile "${CODE_DIR}/sample_balanced_snps.py"
python -c "import numpy, pandas"

Rscript -e '
required <- c("ggplot2", "ggrepel", "scales")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
    stop(paste("Missing R packages:", paste(missing, collapse = ", ")))
}
'

cd "${CODE_DIR}"

bash "${PIPELINE_SCRIPT}"

exit_code=$?
end_time="$(date '+%Y-%m-%d %H:%M:%S')"

echo ""
echo "============================================================="
echo "BALANCED PCA JOB SUMMARY"
echo "Start       : ${start_time}"
echo "End         : ${end_time}"
echo "Exit code   : ${exit_code}"
echo "Job ID      : ${JOB_ID:-not_available}"
echo "============================================================="

exit "${exit_code}"
