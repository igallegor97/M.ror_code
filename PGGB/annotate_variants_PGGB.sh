#!/bin/bash

#$ -N pggb_region_annotation         # Job name
#$ -q all.q                          # Queue
#$ -cwd                              # Run from the current working directory
#$ -pe smp 8                         # 8 CPUs per task (bedtools is I/O-bound)
#$ -V                                # Export environment variables
#$ -t 1-10                           # Array job: one task per community (community.0 to community.9)

# =============================================================
# PGGB — Variant annotation by genomic region
# =============================================================

set -euo pipefail

mkdir -p logs

# -------------------------------------------------------------
# 1. Environment
# -------------------------------------------------------------
module purge
module load Python/3.13.5
module load bedtools/2.28.0
module load Samtools/1.22       # Used to generate .fai files if they do not exist

command -v bedtools >/dev/null 2>&1 || {
    echo "[ERROR] bedtools is not available in PATH"
    exit 1
}

command -v python >/dev/null 2>&1 || {
    echo "[ERROR] python is not available in PATH"
    exit 1
}

python -c "import pandas" 2>/dev/null || {
    echo "[ERROR] pandas is not installed"
    exit 1
}

echo "bedtools : $(bedtools --version | head -1)"
echo "python   : $(python --version)"
echo "samtools : $(samtools --version | head -1)"

# -------------------------------------------------------------
# 2. Community list — one-based indexing for the array job
# -------------------------------------------------------------
COMMUNITY_LIST=(
    community.0
    community.1
    community.2
    community.3
    community.4
    community.5
    community.6
    community.7
    community.8
    community.9
)

# SGE_TASK_ID starts at 1, whereas Bash arrays start at 0
COMMUNITY="${COMMUNITY_LIST[$((SGE_TASK_ID - 1))]}"

# -------------------------------------------------------------
# 3. Environment variables for the Python script
# -------------------------------------------------------------
export GFF3_DIR="/Storage/data1/isabella.gallego/MAESTRIA/data/pacbio"
export PGGB_BASE="/Storage/data1/isabella.gallego/MAESTRIA/PGGB_results_pacbio/pggb_partitioned_results"
export GRAPH_HASH="bf3285f"
export PER_REF_BASE="/Storage/data1/isabella.gallego/MAESTRIA/PGGB_results_pacbio/results_VCF/per_reference_vcfs"
export OUTPUT_DIR="/Storage/data1/isabella.gallego/MAESTRIA/PGGB_results_pacbio/results_VCF/region_annotation"
export COMMUNITY="$COMMUNITY"

PYTHON_SCRIPT="/Storage/data1/isabella.gallego/MAESTRIA/code/annotate_variants_by_region.py"

# -------------------------------------------------------------
# 4. Diagnostic information
# -------------------------------------------------------------
echo "======================================================"
echo "Job ID       : $JOB_ID  (task $SGE_TASK_ID / ${#COMMUNITY_LIST[@]})"
echo "Node         : $(hostname)"
echo "Start time   : $(date '+%Y-%m-%d %H:%M:%S')"
echo "Community    : $COMMUNITY"
echo "GFF3_DIR     : $GFF3_DIR"
echo "PGGB_BASE    : $PGGB_BASE"
echo "PER_REF_BASE : $PER_REF_BASE"
echo "OUTPUT_DIR   : $OUTPUT_DIR"
echo "======================================================"

# -------------------------------------------------------------
# 5. Pre-run checks
# -------------------------------------------------------------
[[ ! -f "$PYTHON_SCRIPT" ]] && {
    echo "[ERROR] Script not found: $PYTHON_SCRIPT"
    exit 1
}

[[ ! -d "$GFF3_DIR" ]] && {
    echo "[ERROR] GFF3_DIR does not exist: $GFF3_DIR"
    exit 1
}

[[ ! -d "$PGGB_BASE" ]] && {
    echo "[ERROR] PGGB_BASE does not exist: $PGGB_BASE"
    exit 1
}

mkdir -p "$OUTPUT_DIR"

# -------------------------------------------------------------
# 6. Run variant annotation
# -------------------------------------------------------------
echo ""
echo "Starting variant annotation by genomic region..."
echo ""

python "$PYTHON_SCRIPT"

EXIT_CODE=$?

# -------------------------------------------------------------
# 7. Job completion
# -------------------------------------------------------------
echo ""
echo "======================================================"
echo "Community     : $COMMUNITY"
echo "Job finished  : $(date '+%Y-%m-%d %H:%M:%S')"
echo "Exit code     : $EXIT_CODE"
echo "Results saved : $OUTPUT_DIR/$COMMUNITY"
echo "======================================================"

exit $EXIT_CODE
