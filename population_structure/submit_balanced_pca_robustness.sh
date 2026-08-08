#!/bin/bash

# =============================================================
# submit_balanced_pca_robustness.sh
#
# Submits the 100-replicate SGE array and then submits the summary
# job with a dependency on the array job.
#
# Usage:
#   bash submit_balanced_pca_robustness.sh
# =============================================================

set -euo pipefail

CODE_DIR="/Storage/data1/isabella.gallego/MAESTRIA/code/population_structure"

cd "${CODE_DIR}"

mkdir -p logs

ARRAY_JOB_ID="$(
    qsub \
        -terse \
        run_balanced_pca_robustness_array.sh
)"

echo "Array job submitted: ${ARRAY_JOB_ID}"

SUMMARY_JOB_ID="$(
    qsub \
        -terse \
        -hold_jid "${ARRAY_JOB_ID}" \
        run_balanced_pca_robustness_summary.sh
)"

echo "Summary job submitted: ${SUMMARY_JOB_ID}"
echo "The summary job will start after the array job finishes."
