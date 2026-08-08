#!/bin/bash

#$ -N balanced_pca_rep
#$ -q all.q
#$ -cwd
#$ -pe smp 2
#$ -V
#$ -t 1-100
#$ -o logs/robustness_rep_$TASK_ID_$JOB_ID.out
#$ -e logs/robustness_rep_$TASK_ID_$JOB_ID.err

# =============================================================
# Balanced PCA robustness replicates — SGE array job
#
# Each array task runs one random seed.
#
# Default:
#   100 replicates
#   1,000 SNPs per method and replicate
#
# Usage:
#   qsub run_balanced_pca_robustness_array.sh
# =============================================================

set -euo pipefail

# =========================
# CONFIGURATION
# =========================

CODE_DIR="/Storage/data1/isabella.gallego/MAESTRIA/code/population_structure"

POPULATION_BASE="/Storage/data1/isabella.gallego/MAESTRIA/population_structure"

SCRIPT="${CODE_DIR}/run_balanced_pca_seed.py"

METADATA="${CODE_DIR}/sample_metadata.tsv"

TOTAL_SNPS="${TOTAL_SNPS:-1000}"

# The task ID itself is used as the reproducible random seed.
SEED="${SGE_TASK_ID}"

OUTPUT_DIR="${POPULATION_BASE}/robustness/seed_${SEED}"

PGGB_MATRIX="${POPULATION_BASE}/PGGB/PGGB_snp_matrix_imputed.tsv"
PGGB_FEATURES="${POPULATION_BASE}/PGGB/PGGB_snp_feature_metadata.tsv"

CACTUS_MATRIX="${POPULATION_BASE}/Cactus/Cactus_snp_matrix_imputed.tsv"
CACTUS_FEATURES="${POPULATION_BASE}/Cactus/Cactus_snp_feature_metadata.tsv"

THREADS="${NSLOTS:-2}"

# =========================
# ENVIRONMENT
# =========================

module purge
module load Python/3.13.5

export OMP_NUM_THREADS="${THREADS}"
export OPENBLAS_NUM_THREADS="${THREADS}"
export MKL_NUM_THREADS="${THREADS}"
export NUMEXPR_NUM_THREADS="${THREADS}"

# =========================
# PRE-RUN CHECKS
# =========================

mkdir -p "${CODE_DIR}/logs"

for required_file in \
    "${SCRIPT}" \
    "${METADATA}" \
    "${PGGB_MATRIX}" \
    "${PGGB_FEATURES}" \
    "${CACTUS_MATRIX}" \
    "${CACTUS_FEATURES}"
do
    if [[ ! -f "${required_file}" ]]; then
        echo "[ERROR] Required file not found:"
        echo "        ${required_file}"
        exit 1
    fi
done

python -m py_compile "${SCRIPT}"

python -c "import numpy, pandas" 2>/dev/null || {
    echo "[ERROR] Python packages numpy and pandas are required."
    exit 1
}

# =========================
# RUN REPLICATE
# =========================

echo "============================================================="
echo "BALANCED PCA ROBUSTNESS ARRAY TASK"
echo "Job ID      : ${JOB_ID:-not_available}"
echo "Task ID     : ${SGE_TASK_ID}"
echo "Seed        : ${SEED}"
echo "Total SNPs  : ${TOTAL_SNPS}"
echo "Node        : $(hostname)"
echo "Output      : ${OUTPUT_DIR}"
echo "============================================================="

python "${SCRIPT}" \
    --seed "${SEED}" \
    --total-snps "${TOTAL_SNPS}" \
    --metadata "${METADATA}" \
    --pggb-matrix "${PGGB_MATRIX}" \
    --pggb-features "${PGGB_FEATURES}" \
    --cactus-matrix "${CACTUS_MATRIX}" \
    --cactus-features "${CACTUS_FEATURES}" \
    --output-dir "${OUTPUT_DIR}"
