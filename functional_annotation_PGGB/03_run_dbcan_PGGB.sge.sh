#!/usr/bin/env bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -j y
#$ -N pggb_func_dbcan
#$ -pe smp 8

# =============================================================
# CAZyme annotation with run_dbCAN 4.0
# =============================================================

set -euo pipefail

: "${PROJECT_ROOT:?PROJECT_ROOT was not exported by the submission driver}"
source "$PROJECT_ROOT/pipeline_config_PGGB.sh"
initialize_array_sample

RUN_DBCAN="$DBCAN_ENV/bin/run_dbcan"
OUTPUT_DIR="$RESULTS/03_dbcan/$SAMPLE"
mkdir -p "$OUTPUT_DIR"

help_text="$("$RUN_DBCAN" CAZyme_annotation --help 2>&1)"
has_option() { grep -q -- "$1" <<< "$help_text"; }

arguments=(CAZyme_annotation)

if has_option --input_raw_data; then
    arguments+=(--input_raw_data "$QC_FASTA")
elif has_option --input; then
    arguments+=(--input "$QC_FASTA")
else
    echo "ERROR: unsupported run_dbCAN input option." >&2
    exit 2
fi

if has_option --mode; then
    arguments+=(--mode protein)
elif has_option --input_type; then
    arguments+=(--input_type protein)
fi

if has_option --output_dir; then
    arguments+=(--output_dir "$OUTPUT_DIR")
elif has_option --out_dir; then
    arguments+=(--out_dir "$OUTPUT_DIR")
else
    echo "ERROR: unsupported run_dbCAN output option." >&2
    exit 2
fi

if has_option --db_dir; then
    arguments+=(--db_dir "$DBCAN_DB")
elif has_option --dbcan_db_dir; then
    arguments+=(--dbcan_db_dir "$DBCAN_DB")
else
    echo "ERROR: unsupported run_dbCAN database option." >&2
    exit 2
fi

if has_option --threads; then
    arguments+=(--threads "${NSLOTS:-$THREADS}")
elif has_option --cpu; then
    arguments+=(--cpu "${NSLOTS:-$THREADS}")
fi

echo "============================================================="
echo "RUN_DBCAN"
echo "Sample       : $SAMPLE"
echo "Input        : $QC_FASTA"
echo "Database     : $DBCAN_DB"
echo "============================================================="

"$RUN_DBCAN" "${arguments[@]}"

echo "RUN_DBCAN COMPLETED: $SAMPLE"

