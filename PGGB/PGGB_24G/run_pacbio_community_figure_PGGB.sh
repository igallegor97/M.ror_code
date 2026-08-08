#!/bin/bash

#$ -N pggb_24G_comm_figure
#$ -q all.q
#$ -cwd
#$ -pe smp 2
#$ -V
#$ -o /Storage/data1/isabella.gallego/MAESTRIA/code/PGGB_24_genomes/logs/24G_comm_figure_$JOB_ID.out
#$ -e /Storage/data1/isabella.gallego/MAESTRIA/code/PGGB_24_genomes/logs/24G_comm_figure_$JOB_ID.err

# =============================================================
# PacBio chromosome-to-community figure — SGE job
#
# Creates a two-panel publication figure:
#
#   A. Group1-Group11 -> PGGB communities
#   B. Ungrouped PacBio sequences -> PGGB communities
#
# Usage:
#   qsub run_pacbio_community_figure.sh
# =============================================================

set -euo pipefail

# =========================
# CONFIGURATION
# =========================

CODE_DIR="/Storage/data1/isabella.gallego/MAESTRIA/code/PGGB_24_genomes"

R_SCRIPT="${CODE_DIR}/plot_pacbio_community_chromosome_map.R"

export MAP_TSV="/Storage/data1/isabella.gallego/MAESTRIA/PGGB_results_24G/community_sequence_tables/pacbio_chromosome_to_communities.tsv"

export COMMUNITY_SUMMARY="/Storage/data1/isabella.gallego/MAESTRIA/PGGB_results_24G/community_sequence_tables/community_summary.tsv"

export OUTPUT_DIR="/Storage/data1/isabella.gallego/MAESTRIA/PGGB_results_24G/community_sequence_figures"

# =========================
# ENVIRONMENT
# =========================

module purge
module load R/4.5.1

# =========================
# PRE-RUN CHECKS
# =========================

mkdir -p "${OUTPUT_DIR}"

if [[ ! -f "${R_SCRIPT}" ]]; then
    echo "[ERROR] R script not found:"
    echo "        ${R_SCRIPT}"
    exit 1
fi

if [[ ! -f "${MAP_TSV}" ]]; then
    echo "[ERROR] Mapping table not found:"
    echo "        ${MAP_TSV}"
    exit 1
fi

if [[ ! -f "${COMMUNITY_SUMMARY}" ]]; then
    echo "[ERROR] Community summary not found:"
    echo "        ${COMMUNITY_SUMMARY}"
    exit 1
fi

Rscript -e '
if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Required R package is missing: ggplot2")
}
'

# =========================
# RUN INFORMATION
# =========================

echo "============================================================="
echo "PACBIO CHROMOSOME-TO-COMMUNITY FIGURE"
echo "Start             : $(date '+%Y-%m-%d %H:%M:%S')"
echo "Job ID            : ${JOB_ID:-not_available}"
echo "Node              : $(hostname)"
echo "MAP_TSV           : ${MAP_TSV}"
echo "COMMUNITY_SUMMARY : ${COMMUNITY_SUMMARY}"
echo "OUTPUT_DIR        : ${OUTPUT_DIR}"
echo "R                 : $(Rscript --version 2>&1 | head -1)"
echo "============================================================="
echo ""

# =========================
# CREATE FIGURE
# =========================

Rscript "${R_SCRIPT}"

exit_code=$?

# =========================
# FINAL REPORT
# =========================

echo ""
echo "============================================================="
echo "FIGURE JOB COMPLETED"
echo "End       : $(date '+%Y-%m-%d %H:%M:%S')"
echo "Exit code : ${exit_code}"
echo "Results   : ${OUTPUT_DIR}"
echo "============================================================="

exit "${exit_code}"
