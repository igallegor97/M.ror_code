#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "$ROOT/pipeline_environment_PGGB.env"

load_modules_init() {
  if ! type module >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    source /etc/profile.d/modules.sh 2>/dev/null || true
  fi
  type module >/dev/null 2>&1 || { echo "ERROR: comando module no disponible" >&2; exit 1; }
}

sample_line() {
  local n="${SGE_TASK_ID:?SGE_TASK_ID no definido}"
  awk -F '\t' -v n="$n" 'NR==n+1 {print $1 "\t" $2}' "$SAMPLES_TSV"
}

init_sample() {
  local line
  line="$(sample_line)"
  IFS=$'\t' read -r SAMPLE INPUT_FASTA <<< "$line"
  [[ -n "${SAMPLE:-}" && -n "${INPUT_FASTA:-}" ]] || { echo "ERROR: tarea $SGE_TASK_ID sin muestra" >&2; exit 1; }
  QC_FASTA="$RESULTS/00_qc/fasta/${SAMPLE}.faa"
}

count_samples() { awk -F '\t' 'NR>1 && $1!="" {n++} END{print n+0}' "$SAMPLES_TSV"; }
