#!/bin/bash

#$ -N pggb_24G_inventory
#$ -q all.q
#$ -cwd
#$ -pe smp 2
#$ -V
#$ -o /Storage/data1/isabella.gallego/MAESTRIA/code/PGGB_24_genomes/logs/24G_inventory_$JOB_ID.out
#$ -e /Storage/data1/isabella.gallego/MAESTRIA/code/PGGB_24_genomes/logs/24G_inventory_$JOB_ID.err

# =============================================================
# PGGB 24-genome community sequence inventory — SGE job
#
# Reads the partition-before-pggb community FASTA files and creates:
#
#   1. community_sequence_inventory.tsv
#   2. community_summary.tsv
#   3. pacbio_chromosome_to_communities.tsv
#
# Additional outputs:
#
#   - candidate_communities_all_5_pacbio.tsv
#   - inventory_validation_summary.tsv
#
# Usage:
#   qsub run_community_sequence_inventory.sh
# =============================================================

set -euo pipefail

# =========================
# CONFIGURATION
# =========================

CODE_DIR="/Storage/data1/isabella.gallego/MAESTRIA/code/PGGB_24_genomes"

PYTHON_SCRIPT="${CODE_DIR}/build_community_sequence_tables.py"

export PARTITION_DIR="/Storage/data1/isabella.gallego/MAESTRIA/PGGB_results_24G/24G_partition_results"

export OUTPUT_DIR="/Storage/data1/isabella.gallego/MAESTRIA/PGGB_results_24G/community_sequence_tables"

export PACBIO_SAMPLES="B3,C26,CO8,CO84,E7"

export COMMUNITY_GLOB="*.community.*.fa"

# =========================
# ENVIRONMENT
# =========================

module purge
module load Python/3.13.5

# =========================
# PRE-RUN CHECKS
# =========================

if [[ ! -f "${PYTHON_SCRIPT}" ]]; then
    echo "[ERROR] Python script not found:"
    echo "        ${PYTHON_SCRIPT}"
    exit 1
fi

if [[ ! -d "${PARTITION_DIR}" ]]; then
    echo "[ERROR] PARTITION_DIR does not exist:"
    echo "        ${PARTITION_DIR}"
    exit 1
fi

python -m py_compile "${PYTHON_SCRIPT}"

mkdir -p "${OUTPUT_DIR}"

# =========================
# RUN INFORMATION
# =========================

echo "============================================================="
echo "PGGB 24-GENOME COMMUNITY INVENTORY"
echo "Start          : $(date '+%Y-%m-%d %H:%M:%S')"
echo "Job ID         : ${JOB_ID:-not_available}"
echo "Node           : $(hostname)"
echo "PARTITION_DIR  : ${PARTITION_DIR}"
echo "OUTPUT_DIR     : ${OUTPUT_DIR}"
echo "PACBIO_SAMPLES : ${PACBIO_SAMPLES}"
echo "Python         : $(python --version 2>&1)"
echo "============================================================="
echo ""

# =========================
# RUN INVENTORY
# =========================

python "${PYTHON_SCRIPT}"

exit_code=$?

# =========================
# FINAL REPORT
# =========================

echo ""
echo "============================================================="
echo "PGGB COMMUNITY INVENTORY JOB COMPLETED"
echo "End        : $(date '+%Y-%m-%d %H:%M:%S')"
echo "Exit code  : ${exit_code}"
echo "Results    : ${OUTPUT_DIR}"
echo "============================================================="

exit "${exit_code}"
