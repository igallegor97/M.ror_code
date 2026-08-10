#!/usr/bin/env bash

# =============================================================
# Functional annotation pipeline submission driver
#
# Submits QC, eggNOG, Pfam, dbCAN, SignalP and integration
# jobs with SGE array dependencies.
#
# Usage:
#   bash run_functional_annotation_pipeline_PGGB.sh
# =============================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$SCRIPT_DIR/pipeline_config_PGGB.sh"

mkdir -p "$LOGS" "$RESULTS"

sample_count="$(count_samples)"
[[ "$sample_count" -eq 5 ]] || {
    echo "ERROR: run 00_preflight_functional_annotation_PGGB.sh first." >&2
    exit 1
}

queue_options=()
[[ -n "${QUEUE:-}" ]] && queue_options=(-q "$QUEUE")

submit_job() {
    qsub -terse -o "$LOGS" -e "$LOGS" "${queue_options[@]}" "$@"
}

# qsub -terse returns array IDs such as 25581.1-5:1.
# This SGE installation requires the base ID (25581) for -hold_jid.
base_job_id() {
    printf '%s\n' "${1%%.*}"
}

echo "============================================================="
echo "FUNCTIONAL ANNOTATION PIPELINE"
echo "Start        : $(date '+%Y-%m-%d %H:%M:%S')"
echo "Samples      : $sample_count"
echo "Results      : $RESULTS"
echo "Logs         : $LOGS"
echo "============================================================="

qc_raw="$(submit_job -t "1-$sample_count" "$PROJECT_ROOT/00_run_proteome_qc_PGGB.sge.sh")"
qc_job="$(base_job_id "$qc_raw")"

eggnog_raw="$(submit_job -hold_jid "$qc_job" -t "1-$sample_count" "$PROJECT_ROOT/01_run_eggnog_mapper_PGGB.sge.sh")"
eggnog_job="$(base_job_id "$eggnog_raw")"

pfam_raw="$(submit_job -hold_jid "$qc_job" -t "1-$sample_count" "$PROJECT_ROOT/02_run_hmmscan_pfam_PGGB.sge.sh")"
pfam_job="$(base_job_id "$pfam_raw")"

dbcan_raw="$(submit_job -hold_jid "$qc_job" -t "1-$sample_count" "$PROJECT_ROOT/03_run_dbcan_PGGB.sge.sh")"
dbcan_job="$(base_job_id "$dbcan_raw")"

signalp_raw="$(submit_job -hold_jid "$qc_job" -t "1-$sample_count" "$PROJECT_ROOT/04_run_signalp5_PGGB.sge.sh")"
signalp_job="$(base_job_id "$signalp_raw")"

integration_raw="$(submit_job -hold_jid "$eggnog_job,$pfam_job,$dbcan_job,$signalp_job" "$PROJECT_ROOT/05_integrate_functional_annotations_PGGB.sge.sh")"
integration_job="$(base_job_id "$integration_raw")"

echo ""
echo "============================================================="
echo "ALL JOBS SUBMITTED"
echo "QC           : $qc_job"
echo "eggNOG       : $eggnog_job"
echo "Pfam         : $pfam_job"
echo "dbCAN        : $dbcan_job"
echo "SignalP      : $signalp_job"
echo "Integration  : $integration_job"
echo "Monitor with : qstat"
echo "============================================================="

