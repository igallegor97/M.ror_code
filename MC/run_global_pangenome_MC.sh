#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -pe smp 10
#$ -t 1

# =============================================================
# Cactus Pangenome Pipeline
#
# Runs the Cactus pangenome workflow using a Singularity container.
#
# Input:
#   Genome list file (seqFile format)
#
# Output:
#   Complete pangenome graph including:
#     - GBZ
#     - GFA
#     - VCF
#     - ODGI graphs
#     - Chromosome-level VG graphs
#     - Giraffe indexes
#     - Visualization files
#
# Usage:
#   qsub run_cactus_pangenome.sh
#
# Container:
#   quay.io/comparative-genomics-toolkit/cactus:latest
# =============================================================

set -euo pipefail

# =========================
# CONFIGURATION
# =========================

GENOME_LIST="data/Mror_Genomes.txt"

OUTPUT_DIR="results/cactus_out"

WORK_DIR="tempdir"

LOG_FILE="PGv1.log"

OUTPUT_NAME="Mror_Pangenome_v1"

REFERENCE_GENOME="Mror_1466_REFE"

THREADS=10

MINIGRAPH_MEMORY="100Gi"

CONTAINER_IMAGE="docker://quay.io/comparative-genomics-toolkit/cactus:latest"

# =========================
# LOAD MODULES
# =========================

module load singularity-ce/3.11.2

# =========================
# RUN INFORMATION
# =========================

echo "============================================================="
echo "CACTUS PANGENOME PIPELINE"
echo "Start            : $(date '+%Y-%m-%d %H:%M:%S')"
echo "Job ID           : ${JOB_ID:-not_available}"
echo "Task ID          : ${SGE_TASK_ID:-not_available}"
echo "Node             : $(hostname)"
echo "GENOME_LIST      : $GENOME_LIST"
echo "OUTPUT_DIR       : $OUTPUT_DIR"
echo "WORK_DIR         : $WORK_DIR"
echo "OUTPUT_NAME      : $OUTPUT_NAME"
echo "REFERENCE_GENOME : $REFERENCE_GENOME"
echo "THREADS          : $THREADS"
echo "MEMORY           : $MINIGRAPH_MEMORY"
echo "CONTAINER        : $CONTAINER_IMAGE"
echo "============================================================="
echo ""

# =========================
# RUN CACTUS
# =========================

singularity exec \
    -H "$(pwd)" \
    "$CONTAINER_IMAGE" \
    cactus-pangenome \
        ./tt \
        "$GENOME_LIST" \
        --outDir "$OUTPUT_DIR" \
        --logFile "$LOG_FILE" \
        --outName "$OUTPUT_NAME" \
        --reference "$REFERENCE_GENOME" \
        --giraffe \
        --viz \
        --odgi \
        --chrom-vg \
        --chrom-og \
        --gbz \
        --gfa \
        --vcf \
        --workDir "$WORK_DIR" \
        --consCores "$THREADS" \
        --mgMemory "$MINIGRAPH_MEMORY"

# =========================
# FINAL REPORT
# =========================

echo ""
echo "============================================================="
echo "CACTUS PANGENOME PIPELINE COMPLETED"
echo "End          : $(date '+%Y-%m-%d %H:%M:%S')"
echo "Results      : $OUTPUT_DIR"
echo "Log file     : $LOG_FILE"
echo "============================================================="
