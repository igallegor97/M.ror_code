#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/pipeline_config_PGGB.sh"
load_modules_init
command -v python3 >/dev/null || { echo "ERROR: python3 no disponible"; exit 1; }
[[ -r "$SAMPLES_TSV" ]] || { echo "ERROR: no se lee $SAMPLES_TSV"; exit 1; }
N="$(count_samples)"; [[ "$N" -eq 5 ]] || { echo "ERROR: se esperaban 5 muestras; hay $N"; exit 1; }
awk -F '\t' 'NR==1 {if($1!="sample_id"||$2!="proteome_fasta") exit 1} NR>1 {if($1!~/^[A-Za-z0-9_.-]+$/||seen[$1]++) exit 1}' "$SAMPLES_TSV" || { echo "ERROR: manifiesto inválido"; exit 1; }
while IFS=$'\t' read -r s f; do [[ "$s" == sample_id ]] && continue; [[ -r "$f" ]] || { echo "ERROR: no se lee $f ($s)"; exit 1; }; done < "$SAMPLES_TSV"

module purge; module load "$EGGNOG_MODULE"; command -v emapper.py; emapper.py --version
[[ -r "$EGGNOG_DATA/eggnog.db" && -r "$EGGNOG_DATA/eggnog_proteins.dmnd" ]] || { echo "ERROR: base eggNOG incompleta"; exit 1; }
module purge; module load "$HMMER_MODULE"; command -v hmmscan; hmmscan -h | head -n 2
[[ -r "$PFAM_HMM" ]] || { echo "ERROR: falta $PFAM_HMM"; exit 1; }
for x in h3f h3i h3m h3p; do [[ -r "$PFAM_HMM.$x" ]] || { echo "ERROR: falta índice $PFAM_HMM.$x"; exit 1; }; done
[[ -x "$DBCAN_ENV/bin/run_dbcan" ]] || { echo "ERROR: no ejecutable run_dbcan"; exit 1; }
"$DBCAN_ENV/bin/run_dbcan" version || true
[[ -r "$DBCAN_DB/dbCAN.hmm" && -r "$DBCAN_DB/dbCAN-sub.hmm" ]] || { echo "ERROR: base dbCAN incompleta"; exit 1; }
module purge; module load "$SIGNALP_MODULE"; command -v signalp; signalp -version
echo "OK: cinco muestras, módulos y bases verificados."
