#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/pipeline_config_PGGB.sh"
mkdir -p "$LOGS" "$RESULTS"
N="$(count_samples)"; [[ "$N" -eq 5 ]] || { echo "ERROR: run 00_preflight_functional_annotation_PGGB.sh first" >&2; exit 1; }
Q=(); [[ -n "$QUEUE" ]] && Q=(-q "$QUEUE")
submit(){ qsub -terse -o "$LOGS" -e "$LOGS" "${Q[@]}" "$@"; }
JQC="$(submit -t 1-$N 00_run_proteome_qc_PGGB.sge.sh)"
JE="$(submit -hold_jid "$JQC" -t 1-$N 01_run_eggnog_mapper_PGGB.sge.sh)"
JP="$(submit -hold_jid "$JQC" -t 1-$N 02_run_hmmscan_pfam_PGGB.sge.sh)"
JD="$(submit -hold_jid "$JQC" -t 1-$N 03_run_dbcan_PGGB.sge.sh)"
JS="$(submit -hold_jid "$JQC" -t 1-$N 04_run_signalp5_PGGB.sge.sh)"
JI="$(submit -hold_jid "$JE,$JP,$JD,$JS" 05_integrate_functional_annotations_PGGB.sge.sh)"
printf 'QC=%s\neggNOG=%s\nPfam=%s\ndbCAN=%s\nSignalP=%s\nIntegration=%s\n' "$JQC" "$JE" "$JP" "$JD" "$JS" "$JI"
