#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/pipeline_config_PGGB.sh"
OUT="$RESULTS/06_external"; mkdir -p "$OUT"
awk -F '\t' 'NR>1{print $1}' "$SAMPLES_TSV" | while read -r s; do
  cp "$RESULTS/00_qc/fasta/${s}.faa" "$OUT/${s}.all_proteins.faa"
  mature="$(find "$RESULTS/04_signalp/$s" -maxdepth 1 -type f \( -name '*mature*.fasta' -o -name '*mature*.fa' \) | head -n1 || true)"
  if [[ -n "$mature" ]]; then cp "$mature" "$OUT/${s}.signalp_secreted.faa"; else : > "$OUT/${s}.signalp_secreted.faa.MISSING"; fi
done
IDS="$PROJECT_ROOT/priority_gene_ids_PGGB.txt"
awk '!/^#/ && NF{print $1}' "$IDS" > "$OUT/.priority_ids.clean"
if [[ -s "$OUT/.priority_ids.clean" ]]; then
  awk 'NR==FNR{want[$1]=1;next} /^>/{id=substr($1,2); keep=(id in want)} keep' "$OUT/.priority_ids.clean" "$RESULTS"/00_qc/fasta/*.faa > "$OUT/priority_candidates.faa"
fi
