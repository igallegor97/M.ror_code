#!/usr/bin/env bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -j y
#$ -N pggb_func_integrate
#$ -pe smp 1

# =============================================================
# Functional-annotation normalization and integration
# =============================================================

set -euo pipefail

: "${PROJECT_ROOT:?PROJECT_ROOT was not exported by the submission driver}"
source "$PROJECT_ROOT/pipeline_config_PGGB.sh"

echo "============================================================="
echo "FUNCTIONAL ANNOTATION INTEGRATION"
echo "Start        : $(date '+%Y-%m-%d %H:%M:%S')"
echo "Results      : $RESULTS"
echo "============================================================="

python3 "$PROJECT_ROOT/05_normalize_functional_annotations_PGGB.py" \
    --root "$RESULTS" \
    --samples "$SAMPLES_TSV"

bash "$PROJECT_ROOT/05_prepare_external_services_PGGB.sh"

echo "============================================================="
echo "FUNCTIONAL ANNOTATION INTEGRATION COMPLETED"
echo "End          : $(date '+%Y-%m-%d %H:%M:%S')"
echo "Master table : $RESULTS/07_master/functional_annotation_master.tsv"
echo "============================================================="

