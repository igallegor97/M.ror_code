#!/bin/bash

# =============================================================
# PGGB — Per-reference VCF generation for reference-bias analysis
#
# Description:
#   Runs generate_per_reference_vcfs.py for one PGGB community.
#   The pipeline generates one VCF for each genome used as the
#   reference path in the graph.
#
# Usage:
#   bash run_per_reference_vcfs_local.sh
#
# Requirements:
#   vg
#   Python 3
#   generate_per_reference_vcfs.py
# =============================================================

set -euo pipefail

# =========================
# CONFIGURATION
# =========================

COMMUNITY="community.8"
GRAPH_HASH="bf3285f"

BASE_DIR="/home/isabella_gallego/OneDrive/Documentos/Maestria/PGGB/comunidades"

COMMUNITY_DIR="${BASE_DIR}/all_pacbio_pansn.fasta.${GRAPH_HASH}.${COMMUNITY}"

FASTA_DIR=(
    "/home/isabella_gallego/OneDrive/Documentos/Maestria/"
    "cactus/data/results_pacbio/raw_data"
)

OUTPUT_DIR="${BASE_DIR}/per_reference_vcfs/${COMMUNITY}"

PYTHON_SCRIPT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/generate_per_reference_vcfs.py"

# Maximum number of threads to use.
# Set THREADS in the environment to override this value.
DEFAULT_THREADS="$(nproc)"
THREADS="${THREADS:-$DEFAULT_THREADS}"

# File patterns used to identify graph files
GFA_PATTERN="*.gfa"
XG_PATTERN="*.xg"

# =========================
# OPTIONAL CONDA ENVIRONMENT
# =========================

# Uncomment and adjust one of the following lines when vg is installed
# inside a Conda environment.

# source "$HOME/miniconda3/etc/profile.d/conda.sh"
# conda activate pggb_env

# source "$HOME/anaconda3/etc/profile.d/conda.sh"
# conda activate pggb_env

# =========================
# FUNCTIONS
# =========================

find_single_file() {
    local search_directory="$1"
    local file_pattern="$2"
    local file_label="$3"

    local matching_files=()

    shopt -s nullglob
    matching_files=(
        "$search_directory"/$file_pattern
    )
    shopt -u nullglob

    if [[ ${#matching_files[@]} -eq 0 ]]; then
        echo "[ERROR] No ${file_label} file matching '${file_pattern}' was found in:"
        echo "        $search_directory"
        return 1
    fi

    if [[ ${#matching_files[@]} -gt 1 ]]; then
        echo "[WARNING] Multiple ${file_label} files were found:"
        printf "          %s\n" "${matching_files[@]}"
        echo "          Using: ${matching_files[0]}"
    fi

    printf "%s\n" "${matching_files[0]}"
}


validate_positive_integer() {
    local value="$1"
    local variable_name="$2"

    if [[ ! "$value" =~ ^[1-9][0-9]*$ ]]; then
        echo "[ERROR] ${variable_name} must be a positive integer: $value"
        exit 1
    fi
}

# =========================
# PRE-RUN CHECKS
# =========================

start_time="$(date '+%Y-%m-%d %H:%M:%S')"

if [[ ! -d "$COMMUNITY_DIR" ]]; then
    echo "[ERROR] Community directory does not exist:"
    echo "        $COMMUNITY_DIR"
    exit 1
fi

if [[ ! -d "$FASTA_DIR" ]]; then
    echo "[ERROR] FASTA_DIR does not exist or is not a directory:"
    echo "        $FASTA_DIR"
    exit 1
fi

if [[ ! -f "$PYTHON_SCRIPT" ]]; then
    echo "[ERROR] Python script not found:"
    echo "        $PYTHON_SCRIPT"
    exit 1
fi

command -v vg >/dev/null 2>&1 || {
    echo "[ERROR] vg is not available in PATH."
    exit 1
}

command -v python >/dev/null 2>&1 || {
    echo "[ERROR] python is not available in PATH."
    exit 1
}

validate_positive_integer "$THREADS" "THREADS"

# Locate graph files only after validating the community directory.
GFA_FILE="$(find_single_file "$COMMUNITY_DIR" "$GFA_PATTERN" "GFA")"
XG_FILE="$(find_single_file "$COMMUNITY_DIR" "$XG_PATTERN" "XG")"

export GFA_FILE
export XG_FILE
export FASTA_DIR
export OUTPUT_DIR
export THREADS

mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/logs"

# =========================
# SOFTWARE INFORMATION
# =========================

VG_VERSION="$(
    vg version 2>/dev/null \
        | head -n 1 \
        || true
)"

PYTHON_VERSION="$(
    python --version 2>&1 \
        || true
)"

# =========================
# RUN INFORMATION
# =========================

echo "============================================================="
echo "PER-REFERENCE VCF GENERATION"
echo "Start         : $start_time"
echo "Community     : $COMMUNITY"
echo "Graph hash    : $GRAPH_HASH"
echo "Community dir : $COMMUNITY_DIR"
echo "GFA_FILE      : $GFA_FILE"
echo "XG_FILE       : $XG_FILE"
echo "FASTA_DIR     : $FASTA_DIR"
echo "OUTPUT_DIR    : $OUTPUT_DIR"
echo "THREADS       : $THREADS"
echo "Python script : $PYTHON_SCRIPT"
echo "vg             : ${VG_VERSION:-unknown}"
echo "Python         : ${PYTHON_VERSION:-unknown}"
echo "============================================================="
echo ""

# =========================
# RUN PIPELINE
# =========================

log_file="$OUTPUT_DIR/logs/run_$(date '+%Y%m%d_%H%M%S').log"

echo "Starting per-reference VCF generation..."
echo "Log file: $log_file"
echo ""

set +e

python "$PYTHON_SCRIPT" 2>&1 \
    | tee "$log_file"

pipeline_status=(
    "${PIPESTATUS[@]}"
)

set -e

python_exit_code="${pipeline_status[0]}"
tee_exit_code="${pipeline_status[1]}"

if [[ "$python_exit_code" -ne 0 ]]; then
    exit_code="$python_exit_code"
elif [[ "$tee_exit_code" -ne 0 ]]; then
    exit_code="$tee_exit_code"
else
    exit_code=0
fi

# =========================
# FINAL REPORT
# =========================

end_time="$(date '+%Y-%m-%d %H:%M:%S')"

echo ""
echo "============================================================="
echo "PER-REFERENCE VCF GENERATION SUMMARY"
echo "Community        : $COMMUNITY"
echo "Finished         : $end_time"
echo "Python exit code : $python_exit_code"
echo "tee exit code    : $tee_exit_code"
echo "Final exit code  : $exit_code"
echo "VCFs saved to    : $OUTPUT_DIR"
echo "Log saved to     : $log_file"
echo "============================================================="

if [[ "$exit_code" -eq 0 ]]; then
    echo ""
    echo "Next step:"
    echo ""
    echo "  Run the per-reference variant-classification script on:"
    echo "  $OUTPUT_DIR"
else
    echo ""
    echo "[ERROR] Per-reference VCF generation did not complete successfully."
    echo "        Check the log file:"
    echo "        $log_file"
fi

exit "$exit_code"
