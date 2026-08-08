#!/bin/bash

#$ -N population_structure
#$ -q all.q
#$ -cwd
#$ -pe smp 8
#$ -V
#$ -o logs/population_structure_$JOB_ID.out
#$ -e logs/population_structure_$JOB_ID.err

# =============================================================
# Population-structure pipeline — SGE job
#
# Runs the complete PGGB and minigraph-cactus population-analysis
# pipeline on the HPC cluster.
#
# Analyses:
#   1. Build the PGGB all-SNP matrix.
#   2. Build the Cactus all-SNP matrix.
#   3. Run PCA for both graph methods.
#   4. Calculate Euclidean and Hamming distances.
#   5. Generate heatmaps and UPGMA dendrograms.
#   6. Generate the PGGB-versus-Cactus PCA comparison.
#
# PCA colors:
#   Geographic region from sample_metadata.tsv
#
# Usage:
#   qsub run_population_structure_job.sh
# =============================================================

set -euo pipefail

# =========================
# CONFIGURATION
# =========================

CODE_DIR="/Storage/data1/isabella.gallego/MAESTRIA/code/population_structure"

PIPELINE_SCRIPT="${CODE_DIR}/run_population_structure_pipeline.sh"

METADATA_FILE="${CODE_DIR}/sample_metadata.tsv"

LOG_DIR="${CODE_DIR}/logs"

# Number of CPU cores requested from SGE.
THREADS="${NSLOTS:-8}"

# =========================
# SETUP
# =========================

mkdir -p "${LOG_DIR}"

start_time="$(date '+%Y-%m-%d %H:%M:%S')"

echo "============================================================="
echo "POPULATION-STRUCTURE PIPELINE — SGE JOB"
echo "Start          : ${start_time}"
echo "Job ID         : ${JOB_ID:-not_available}"
echo "Task ID        : ${SGE_TASK_ID:-not_available}"
echo "Node           : $(hostname)"
echo "Working dir    : $(pwd)"
echo "CODE_DIR       : ${CODE_DIR}"
echo "PIPELINE       : ${PIPELINE_SCRIPT}"
echo "METADATA       : ${METADATA_FILE}"
echo "THREADS        : ${THREADS}"
echo "============================================================="
echo ""

# =========================
# LOAD MODULES
# =========================

module purge
module load Python/3.13.5
module load R/4.5.1

echo "Python : $(python --version 2>&1)"
echo "R      : $(Rscript --version 2>&1 | head -1)"
echo ""

# =========================
# PRE-RUN CHECKS
# =========================

if [[ ! -d "${CODE_DIR}" ]]; then
    echo "[ERROR] CODE_DIR does not exist:"
    echo "        ${CODE_DIR}"
    exit 1
fi

if [[ ! -f "${PIPELINE_SCRIPT}" ]]; then
    echo "[ERROR] Pipeline script not found:"
    echo "        ${PIPELINE_SCRIPT}"
    exit 1
fi

if [[ ! -f "${METADATA_FILE}" ]]; then
    echo "[ERROR] Metadata file not found:"
    echo "        ${METADATA_FILE}"
    exit 1
fi

for required_file in \
    "${CODE_DIR}/build_snp_matrix.py" \
    "${CODE_DIR}/analyze_population_structure.R" \
    "${CODE_DIR}/compare_pggb_cactus_pca.R"
do
    if [[ ! -f "${required_file}" ]]; then
        echo "[ERROR] Required pipeline file not found:"
        echo "        ${required_file}"
        exit 1
    fi
done

command -v python >/dev/null 2>&1 || {
    echo "[ERROR] python is not available in PATH."
    exit 1
}

command -v Rscript >/dev/null 2>&1 || {
    echo "[ERROR] Rscript is not available in PATH."
    exit 1
}

python -c "import numpy, pandas" 2>/dev/null || {
    echo "[ERROR] Python packages numpy and pandas are required."
    echo "        Install them in your user environment before submitting the job."
    exit 1
}

# Check the R packages without installing them on the compute node.
Rscript -e '
required <- c("ggplot2", "ggrepel", "scales")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
    stop(
        paste(
            "Missing R packages:",
            paste(missing, collapse = ", ")
        )
    )
}
' || {
    echo "[ERROR] One or more required R packages are unavailable."
    echo "        Required: ggplot2, ggrepel, scales"
    exit 1
}

# =========================
# RUN PIPELINE
# =========================

cd "${CODE_DIR}"

export PYTHON="python"
export RSCRIPT="Rscript"
export OMP_NUM_THREADS="${THREADS}"
export OPENBLAS_NUM_THREADS="${THREADS}"
export MKL_NUM_THREADS="${THREADS}"
export NUMEXPR_NUM_THREADS="${THREADS}"

echo "Starting population-structure pipeline..."
echo ""

bash "${PIPELINE_SCRIPT}"

exit_code=$?

# =========================
# FINAL REPORT
# =========================

end_time="$(date '+%Y-%m-%d %H:%M:%S')"

echo ""
echo "============================================================="
echo "POPULATION-STRUCTURE JOB SUMMARY"
echo "Start       : ${start_time}"
echo "End         : ${end_time}"
echo "Exit code   : ${exit_code}"
echo "Job ID      : ${JOB_ID:-not_available}"
echo "Node        : $(hostname)"
echo "============================================================="

exit "${exit_code}"
