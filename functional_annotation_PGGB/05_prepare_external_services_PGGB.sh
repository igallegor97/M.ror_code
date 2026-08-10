#!/usr/bin/env bash

# =============================================================
# External-service FASTA preparation
#
# Prepares all-protein, SignalP-secreted and selective candidate
# FASTA files for DeepTMHMM, EffectorP 3 and InterProScan.
# =============================================================

set -euo pipefail

if [[ -n "${PROJECT_ROOT:-}" && -r "$PROJECT_ROOT/pipeline_config_PGGB.sh" ]]; then
    source "$PROJECT_ROOT/pipeline_config_PGGB.sh"
else
    SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
    source "$SCRIPT_DIR/pipeline_config_PGGB.sh"
fi

OUTPUT_DIR="$RESULTS/06_external"
mkdir -p "$OUTPUT_DIR"

awk -F '\t' 'NR > 1 { print $1 }' "$SAMPLES_TSV" |
while read -r sample; do
    cp "$RESULTS/00_qc/fasta/${sample}.faa" "$OUTPUT_DIR/${sample}.all_proteins.faa"

    mature_fasta="$(find "$RESULTS/04_signalp/$sample" -maxdepth 1 -type f \
        \( -name '*mature*.fasta' -o -name '*mature*.fa' \) | head -n 1 || true)"

    if [[ -n "$mature_fasta" ]]; then
        cp "$mature_fasta" "$OUTPUT_DIR/${sample}.signalp_secreted.faa"
    else
        touch "$OUTPUT_DIR/${sample}.signalp_secreted.faa.MISSING"
    fi
done

clean_ids="$OUTPUT_DIR/.priority_gene_ids.clean"
awk '!/^#/ && NF { print $1 }' "$PROJECT_ROOT/priority_gene_ids_PGGB.txt" > "$clean_ids"

if [[ -s "$clean_ids" ]]; then
    awk '
        NR == FNR { requested[$1] = 1; next }
        /^>/ {
            identifier = substr($1, 2)
            keep = identifier in requested
        }
        keep
    ' "$clean_ids" "$RESULTS"/00_qc/fasta/*.faa > "$OUTPUT_DIR/priority_candidates.faa"
fi

echo "External-service files saved to: $OUTPUT_DIR"

