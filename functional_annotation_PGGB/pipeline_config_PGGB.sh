#!/usr/bin/env bash

# Shared functions for all functional-annotation jobs.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Export every environment variable so qsub -V transmits it to compute nodes.
set -a
source "$ROOT/pipeline_environment_PGGB.env"
set +a

load_modules_init() {
    if ! type module >/dev/null 2>&1; then
        source /etc/profile.d/modules.sh 2>/dev/null || true
    fi
    type module >/dev/null 2>&1 || {
        echo "ERROR: the environment module command is unavailable." >&2
        exit 1
    }
}

count_samples() {
    awk -F '\t' 'NR > 1 && $1 != "" { n++ } END { print n + 0 }' "$SAMPLES_TSV"
}

initialize_array_sample() {
    local line
    local task_id="${SGE_TASK_ID:?SGE_TASK_ID is not defined}"

    line="$(awk -F '\t' -v n="$task_id" 'NR == n + 1 { print $1 "\t" $2 }' "$SAMPLES_TSV")"
    IFS=$'\t' read -r SAMPLE INPUT_FASTA <<< "$line"

    [[ -n "${SAMPLE:-}" && -n "${INPUT_FASTA:-}" ]] || {
        echo "ERROR: no sample corresponds to SGE task $task_id." >&2
        exit 1
    }

    QC_FASTA="$RESULTS/00_qc/fasta/${SAMPLE}.faa"
}

