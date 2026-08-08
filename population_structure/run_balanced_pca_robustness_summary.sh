#!/bin/bash

#$ -N balanced_pca_summary
#$ -q all.q
#$ -cwd
#$ -pe smp 2
#$ -V
#$ -o logs/robustness_summary_$JOB_ID.out
#$ -e logs/robustness_summary_$JOB_ID.err

# =============================================================
# Balanced PCA robustness summary — SGE job
#
# Run this job after the array job has completed successfully.
#
# Usage:
#   qsub run_balanced_pca_robustness_summary.sh
# =============================================================

set -euo pipefail

# =========================
# CONFIGURATION
# =========================

CODE_DIR="/Storage/data1/isabella.gallego/MAESTRIA/code/population_structure"

POPULATION_BASE="/Storage/data1/isabella.gallego/MAESTRIA/population_structure"

ROBUSTNESS_DIR="${POPULATION_BASE}/robustness"

OUTPUT_DIR="${ROBUSTNESS_DIR}/summary"

SUMMARY_SCRIPT="${CODE_DIR}/summarize_balanced_pca_robustness.R"

FULL_PGGB="${POPULATION_BASE}/PGGB/PGGB_pca_coordinates.tsv"

FULL_CACTUS="${POPULATION_BASE}/Cactus/Cactus_pca_coordinates.tsv"

# =========================
# ENVIRONMENT
# =========================

module purge
module load R/4.5.1

mkdir -p "${CODE_DIR}/logs"
mkdir -p "${OUTPUT_DIR}"

# =========================
# PRE-RUN CHECKS
# =========================

for required_file in \
    "${SUMMARY_SCRIPT}" \
    "${FULL_PGGB}" \
    "${FULL_CACTUS}"
do
    if [[ ! -f "${required_file}" ]]; then
        echo "[ERROR] Required file not found:"
        echo "        ${required_file}"
        exit 1
    fi
done

Rscript -e '
required <- c("ggplot2", "scales")
missing <- required[
    !vapply(required, requireNamespace, logical(1), quietly = TRUE)
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

PGGB_COUNT="$(
    find "${ROBUSTNESS_DIR}" \
        -path "*/PGGB/pca_coordinates.tsv" \
        -type f \
        | wc -l
)"

CACTUS_COUNT="$(
    find "${ROBUSTNESS_DIR}" \
        -path "*/Cactus/pca_coordinates.tsv" \
        -type f \
        | wc -l
)"

echo "PGGB replicates found   : ${PGGB_COUNT}"
echo "Cactus replicates found : ${CACTUS_COUNT}"

if [[ "${PGGB_COUNT}" -eq 0 || "${CACTUS_COUNT}" -eq 0 ]]; then
    echo "[ERROR] No complete robustness replicates were found."
    exit 1
fi

if [[ "${PGGB_COUNT}" -ne "${CACTUS_COUNT}" ]]; then
    echo "[ERROR] PGGB and Cactus replicate counts differ."
    exit 1
fi

# =========================
# RUN SUMMARY
# =========================

Rscript "${SUMMARY_SCRIPT}" \
    "${ROBUSTNESS_DIR}" \
    "${FULL_PGGB}" \
    "${FULL_CACTUS}" \
    "${OUTPUT_DIR}"
