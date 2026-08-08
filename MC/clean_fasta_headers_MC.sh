#!/bin/bash

#$ -q all.q
#$ -V
#$ -cwd

# =============================================================
# FASTA header simplification
#
# Rewrites FASTA headers using a standardized sequential naming
# scheme while preserving the nucleotide or protein sequences.
#
# Input:
#   Directory containing genome FASTA files.
#
# Output:
#   One cleaned FASTA file per input genome.
#
# Header format:
#   >SAMPLE_0001
#   >SAMPLE_0002
#   >SAMPLE_0003
#   ...
#
# Example:
#
#   Original:
#       >chr1 some long description...
#
#   Converted:
#       >Mror_0001
#
# Usage:
#   qsub clean_fasta_headers.sh
# =============================================================

set -euo pipefail

# =========================
# CONFIGURATION
# =========================

INPUT_DIR="/Storage/data1/isabella.gallego/MAESTRIA/data/genomes"

OUTPUT_DIR="/Storage/data1/isabella.gallego/MAESTRIA/data/cleaned_fasta"

# =========================
# SETUP
# =========================

mkdir -p "$OUTPUT_DIR"

echo "============================================================="
echo "FASTA HEADER SIMPLIFICATION"
echo "Start        : $(date '+%Y-%m-%d %H:%M:%S')"
echo "INPUT_DIR    : $INPUT_DIR"
echo "OUTPUT_DIR   : $OUTPUT_DIR"
echo "============================================================="
echo ""

# =========================
# PROCESS FASTA FILES
# =========================

for file in "$INPUT_DIR"/*.fasta; do
    [[ -e "$file" ]] || continue

    echo "Processing: $file"

    # Extract the filename prefix (e.g. Mror from Mror.fasta)
    sample=$(basename "$file" | cut -d '.' -f 1)

    output="$OUTPUT_DIR/${sample}.cleaned.fasta"

    # Rewrite FASTA headers using sequential identifiers
    awk -v prefix="$sample" '
        BEGIN { i = 1 }

        /^>/ {
            printf(">%s_%04d\n", prefix, i++)
            next
        }

        {
            print
        }
    ' "$file" > "$output"

done

# =========================
# FINAL REPORT
# =========================

echo ""
echo "============================================================="
echo "FASTA HEADER SIMPLIFICATION COMPLETED"
echo "End          : $(date '+%Y-%m-%d %H:%M:%S')"
echo "Results saved to:"
echo "  $OUTPUT_DIR"
echo "============================================================="
