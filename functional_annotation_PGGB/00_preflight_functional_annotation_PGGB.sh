#!/usr/bin/env bash

# =============================================================
# Functional annotation preflight validation
#
# Validates the five-proteome manifest, input files, HPC modules
# and shared databases without submitting analysis jobs.
#
# Usage:
#   bash 00_preflight_functional_annotation_PGGB.sh
# =============================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$SCRIPT_DIR/pipeline_config_PGGB.sh"

echo "============================================================="
echo "FUNCTIONAL ANNOTATION PREFLIGHT"
echo "Start        : $(date '+%Y-%m-%d %H:%M:%S')"
echo "PROJECT_ROOT : $PROJECT_ROOT"
echo "RESULTS      : $RESULTS"
echo "MANIFEST     : $SAMPLES_TSV"
echo "============================================================="

command -v python3 >/dev/null || { echo "ERROR: python3 is unavailable." >&2; exit 1; }
load_modules_init

[[ -r "$SAMPLES_TSV" ]] || { echo "ERROR: cannot read $SAMPLES_TSV" >&2; exit 1; }

sample_count="$(count_samples)"
[[ "$sample_count" -eq 5 ]] || {
    echo "ERROR: expected exactly five samples; found $sample_count." >&2
    exit 1
}

awk -F '\t' '
    NR == 1 {
        if ($1 != "sample_id" || $2 != "proteome_fasta") exit 1
        next
    }
    {
        if ($1 !~ /^[A-Za-z0-9_.-]+$/ || seen[$1]++) exit 1
    }
' "$SAMPLES_TSV" || { echo "ERROR: invalid manifest header or sample IDs." >&2; exit 1; }

while IFS=$'\t' read -r sample fasta; do
    [[ "$sample" == "sample_id" ]] && continue
    [[ -r "$fasta" ]] || { echo "ERROR: cannot read $fasta ($sample)." >&2; exit 1; }
    echo "Input OK     : $sample -> $fasta"
done < "$SAMPLES_TSV"

module purge
module load "$EGGNOG_MODULE"
command -v emapper.py
emapper.py --version
[[ -r "$EGGNOG_DATA/eggnog.db" && -r "$EGGNOG_DATA/eggnog_proteins.dmnd" ]] || {
    echo "ERROR: the shared eggNOG database is incomplete." >&2
    exit 1
}

module purge
module load "$HMMER_MODULE"
command -v hmmscan
[[ -r "$PFAM_HMM" ]] || { echo "ERROR: cannot read $PFAM_HMM" >&2; exit 1; }
for suffix in h3f h3i h3m h3p; do
    [[ -r "$PFAM_HMM.$suffix" ]] || { echo "ERROR: missing $PFAM_HMM.$suffix" >&2; exit 1; }
done

[[ -x "$DBCAN_ENV/bin/run_dbcan" ]] || { echo "ERROR: run_dbcan is unavailable." >&2; exit 1; }
[[ -r "$DBCAN_DB/dbCAN.hmm" && -r "$DBCAN_DB/dbCAN-sub.hmm" ]] || {
    echo "ERROR: the shared dbCAN database is incomplete." >&2
    exit 1
}

module purge
module load "$SIGNALP_MODULE"
command -v signalp
signalp -version

echo "============================================================="
echo "FUNCTIONAL ANNOTATION PREFLIGHT COMPLETED"
echo "End          : $(date '+%Y-%m-%d %H:%M:%S')"
echo "Five samples, modules and databases were verified."
echo "============================================================="

