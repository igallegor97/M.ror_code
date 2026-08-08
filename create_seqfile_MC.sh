#!/bin/bash

#$ -N create_seqfile                 
#$ -q all.q                          
#$ -cwd                              
#$ -V                                

# =============================================================
# Sequence file generator
#
# Creates a tab-separated seqFile.txt containing one row per FASTA file.
#
# Output columns:
#   sample_name
#   fasta_path
#
# Expected input files:
#   *.fasta
#
# Example output:
#   B3    /path/to/B3.fasta
#   C26   /path/to/C26.fasta
#
# Usage:
#   qsub create_seqfile.sh
# =============================================================

set -euo pipefail

# =========================
# CONFIGURATION
# =========================

INPUT_DIR="/Storage/data1/isabella.gallego/MAESTRIA/data/cleaned_fasta"

OUTPUT_FILE="/Storage/data1/isabella.gallego/MAESTRIA/data/seqFile.txt"

FILE_PATTERN="*.fasta"

# =========================
# SETUP
# =========================

start_time=$(date '+%Y-%m-%d %H:%M:%S')

echo "============================================================="
echo "SEQUENCE FILE GENERATION"
echo "Start        : $start_time"
echo "Job ID       : ${JOB_ID:-not_available}"
echo "Node         : $(hostname)"
echo "INPUT_DIR    : $INPUT_DIR"
echo "OUTPUT_FILE  : $OUTPUT_FILE"
echo "FILE_PATTERN : $FILE_PATTERN"
echo "============================================================="
echo ""

# =========================
# PRE-RUN CHECKS
# =========================

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "[ERROR] INPUT_DIR does not exist or is not a directory:"
    echo "        $INPUT_DIR"
    exit 1
fi

if ! command -v realpath >/dev/null 2>&1; then
    echo "[ERROR] realpath is not available in PATH."
    exit 1
fi

output_parent=$(dirname "$OUTPUT_FILE")

if [[ ! -d "$output_parent" ]]; then
    echo "[INFO] Creating output directory:"
    echo "       $output_parent"

    mkdir -p "$output_parent"
fi

cd "$INPUT_DIR"

# Prevent an unmatched glob pattern from being interpreted literally.
shopt -s nullglob

fasta_files=($FILE_PATTERN)

if [[ ${#fasta_files[@]} -eq 0 ]]; then
    echo "[ERROR] No FASTA files matching '$FILE_PATTERN' were found in:"
    echo "        $INPUT_DIR"
    exit 1
fi

echo "FASTA files found: ${#fasta_files[@]}"
echo ""

# =========================
# INITIALIZE OUTPUT FILE
# =========================

# Create or overwrite the output file.
: > "$OUTPUT_FILE"

# =========================
# PROCESS FASTA FILES
# =========================

processed_files=0

for fasta_file in "${fasta_files[@]}"; do
    sample_name=$(basename "$fasta_file" .fasta)
    fasta_path=$(realpath "$fasta_file")

    printf "%s\t%s\n" \
        "$sample_name" \
        "$fasta_path" \
        >> "$OUTPUT_FILE"

    echo "  [OK] $sample_name → $fasta_path"

    ((processed_files += 1))
done

# =========================
# FINAL REPORT
# =========================

end_time=$(date '+%Y-%m-%d %H:%M:%S')

echo ""
echo "============================================================="
echo "SEQUENCE FILE SUMMARY"
echo "FASTA files found : ${#fasta_files[@]}"
echo "Entries written   : $processed_files"
echo "Output file       : $OUTPUT_FILE"
echo "End               : $end_time"
echo "============================================================="
