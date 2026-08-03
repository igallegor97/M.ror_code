#!/bin/bash

#$ -N pggb_classify_variants
#$ -q all.q
#$ -cwd
#$ -pe smp 16
#$ -V

# =============================================================
# PGGB — Variant classification by community
# =============================================================

set -euo pipefail   # Exit immediately on error, undefined variable, or broken pipe

# -------------------------------------------------------------
# 0. Create log directory if it does not exist
# -------------------------------------------------------------
mkdir -p logs

# -------------------------------------------------------------
# 1. Load cluster modules (adjust according to your HPC system)
#    Run `module avail python` to see available Python versions
# -------------------------------------------------------------
module purge
module load Python/3.13.5          # or anaconda3, miniforge3, etc.

# -------------------------------------------------------------
# 2. Verify that pandas is available
# -------------------------------------------------------------
python -c "import pandas" 2>/dev/null || {
    echo "[ERROR] pandas is not installed in the active Python environment."
    echo "        Install it with: pip install pandas --user"
    exit 1
}

# -------------------------------------------------------------
# 3. Environment variables — edit here to change paths
# -------------------------------------------------------------
export BASE_DIR="/Storage/data1/isabella.gallego/MAESTRIA/PGGB_results_pacbio/pggb_partitioned_results"
export VCF_NAME="variants.vcf.gz"
export OUTPUT_DIR="/Storage/data1/isabella.gallego/MAESTRIA/PGGB_results_pacbio/results_VCF"
export COMM_PATTERN="*.community.*"   # Glob pattern for community folders

# Path to the Python script (absolute path)
PYTHON_SCRIPT="/Storage/data1/isabella.gallego/MAESTRIA/code/classify_variants_batch.py"

# -------------------------------------------------------------
# 4. Diagnostic information for the log
# -------------------------------------------------------------
echo "======================================================"
echo "Job ID       : $JOB_ID"
echo "Node         : $(hostname)"
echo "Start time   : $(date '+%Y-%m-%d %H:%M:%S')"
echo "Python       : $(which python)"
echo "BASE_DIR     : $BASE_DIR"
echo "OUTPUT_DIR   : $OUTPUT_DIR"
echo "VCF_NAME     : $VCF_NAME"
echo "COMM_PATTERN : $COMM_PATTERN"
echo "Script       : $PYTHON_SCRIPT"
echo "======================================================"

# -------------------------------------------------------------
# 5. Pre-run checks
# -------------------------------------------------------------
if [[ ! -f "$PYTHON_SCRIPT" ]]; then
    echo "[ERROR] Script not found: $PYTHON_SCRIPT"
    exit 1
fi

if [[ ! -d "$BASE_DIR" ]]; then
    echo "[ERROR] BASE_DIR does not exist: $BASE_DIR"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# -------------------------------------------------------------
# 6. Run variant classification
# -------------------------------------------------------------
echo ""
echo "Starting variant classification..."
echo ""

python "$PYTHON_SCRIPT"

EXIT_CODE=$?

# -------------------------------------------------------------
# 7. Job completion
# -------------------------------------------------------------
echo ""
echo "======================================================"
echo "Job finished  : $(date '+%Y-%m-%d %H:%M:%S')"
echo "Exit code     : $EXIT_CODE"
echo "Results saved : $OUTPUT_DIR"
echo "======================================================"

exit $EXIT_CODE
