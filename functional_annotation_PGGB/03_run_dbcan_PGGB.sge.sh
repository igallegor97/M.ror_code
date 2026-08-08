#!/usr/bin/env bash
#$ -cwd
#$ -V
#$ -j y
#$ -N PB_DBCAN
#$ -pe smp 8
set -euo pipefail
source pipeline_config_PGGB.sh; init_sample
RUN="$DBCAN_ENV/bin/run_dbcan"; OUT="$RESULTS/03_dbcan/$SAMPLE"; mkdir -p "$OUT"
HELP="$($RUN CAZyme_annotation --help 2>&1)"
has(){ grep -q -- "$1" <<< "$HELP"; }
ARGS=(CAZyme_annotation)
if has --input_raw_data; then ARGS+=(--input_raw_data "$QC_FASTA"); elif has --input; then ARGS+=(--input "$QC_FASTA"); else echo "ERROR: opción de entrada desconocida; guarde: $RUN CAZyme_annotation --help" >&2; exit 2; fi
if has --mode; then ARGS+=(--mode protein); elif has --input_type; then ARGS+=(--input_type protein); fi
if has --output_dir; then ARGS+=(--output_dir "$OUT"); elif has --out_dir; then ARGS+=(--out_dir "$OUT"); else echo "ERROR: opción de salida desconocida" >&2; exit 2; fi
if has --db_dir; then ARGS+=(--db_dir "$DBCAN_DB"); elif has --dbcan_db_dir; then ARGS+=(--dbcan_db_dir "$DBCAN_DB"); else echo "ERROR: opción de base desconocida" >&2; exit 2; fi
if has --threads; then ARGS+=(--threads "${NSLOTS:-$THREADS}"); elif has --cpu; then ARGS+=(--cpu "${NSLOTS:-$THREADS}"); fi
"$RUN" "${ARGS[@]}"
