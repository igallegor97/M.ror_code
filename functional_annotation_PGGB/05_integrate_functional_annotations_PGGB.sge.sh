#!/usr/bin/env bash
#$ -cwd
#$ -V
#$ -j y
#$ -N PB_INTEGRATE
#$ -pe smp 1
set -euo pipefail
source pipeline_config_PGGB.sh
python3 05_normalize_functional_annotations_PGGB.py --root "$RESULTS" --samples "$SAMPLES_TSV"
bash 05_prepare_external_services_PGGB.sh
