#!/bin/bash

#$ -N pggb_community_sequence_summary   # Job name
#$ -q all.q                             # Queue
#$ -cwd                                 # Run from the current working directory
#$ -pe smp 2                            # Request 2 CPU cores
#$ -V                                   # Export environment variables

# =============================================================
# PGGB — Community sequence summary
#
# Lists all sequence headers found in each community FASTA file
# and combines them into a single tab-separated summary table.
#
# Expected input filename format:
#   *.community.*.fa
#
# Output columns:
#   community
#   sequence_name
#
# Usage:
#   qsub community_sequence_summary.sh
# =============================================================

set -euo pipefail

# =========================
# CONFIGURATION
# =========================

PARTITION_DIR="/Storage/data1/isabella.gallego/MAESTRIA/PGGB_results_pacbio/partition_results"

OUTPUT_FILE="${PARTITION_DIR}/community_sequence_summary.tsv"

FILE_PATTERN="*.community.*.fa"

# =========================
# SETUP
# =========================

start_time=$(date '+%Y-%m-%d %H:%M:%S')

echo "============================================================="
echo "PGGB COMMUNITY SEQUENCE SUMMARY"
echo "Start          : $start_time"
echo "Job ID         : ${JOB_ID:-not_available}"
echo "Node           : $(hostname)"
echo "PARTITION_DIR  : $PARTITION_DIR"
echo "OUTPUT_FILE    : $OUTPUT_FILE"
echo "FILE_PATTERN   : $FILE_PATTERN"
echo "============================================================="
echo ""

# =========================
# PRE-RUN CHECKS
# =========================

if [[ ! -d "$PARTITION_DIR" ]]; then
    echo "[ERROR] PARTITION_DIR does not exist: $PARTITION_DIR"
    exit 1
fi

cd "$PARTITION_DIR"

# Enable nullglob so that an unmatched pattern expands to an empty array
# instead of remaining as the literal pattern.
shopt -s nullglob

community_files=($FILE_PATTERN)

if [[ ${#community_files[@]} -eq 0 ]]; then
    echo "[ERROR] No FASTA files matching '$FILE_PATTERN' were found in:"
    echo "        $PARTITION_DIR"
    exit 1
fi

echo "Community FASTA files found: ${#community_files[@]}"
echo ""

# =========================
# INITIALIZE OUTPUT TABLE
# =========================

printf "community\tsequence_name\n" > "$OUTPUT_FILE"

# =========================
# PROCESS COMMUNITY FILES
# =========================

processed_files=0
total_sequences=0
empty_files=()

for fasta_file in "${community_files[@]}"; do
    community_name=$(basename "$fasta_file" .fa)

    echo "[$community_name]"
    echo "  FASTA: $fasta_file"

    sequence_count=0

    while IFS= read -r header; do
        # Remove the leading ">" character from each FASTA header
        sequence_name="${header#>}"

        printf "%s\t%s\n" \
            "$community_name" \
            "$sequence_name" \
            >> "$OUTPUT_FILE"

        ((sequence_count += 1))
        ((total_sequences += 1))
    done < <(grep '^>' "$fasta_file" || true)

    if [[ $sequence_count -eq 0 ]]; then
        echo "  [WARNING] No FASTA headers were found."
        empty_files+=("$fasta_file")
    else
        echo "  Sequences found: $sequence_count"
    fi

    ((processed_files += 1))

    echo ""
done

# =========================
# FINAL REPORT
# =========================

end_time=$(date '+%Y-%m-%d %H:%M:%S')

echo "============================================================="
echo "COMMUNITY SEQUENCE SUMMARY"
echo "Files found       : ${#community_files[@]}"
echo "Files processed   : $processed_files"
echo "Total sequences   : $total_sequences"
echo "Files without data: ${#empty_files[@]}"
echo "Results saved to  : $OUTPUT_FILE"
echo "End               : $end_time"
echo "============================================================="

if [[ ${#empty_files[@]} -gt 0 ]]; then
    echo ""
    echo "Files without FASTA headers:"

    for empty_file in "${empty_files[@]}"; do
        echo "  $empty_file"
    done
fi
