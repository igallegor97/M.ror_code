#!/bin/bash

# =============================================================
# run_24G_population_pipeline.sh
#
# Runs:
#   1. all_11 profile
#   2. conservative_9 profile
# =============================================================

set -euo pipefail

CODE_DIR="/Storage/data1/isabella.gallego/MAESTRIA/code/population_structure/population_structure_24G_v2"
PGGB_BASE="/Storage/data1/isabella.gallego/MAESTRIA/PGGB_results_24G/pggb_partitioned_results_24G"
REFERENCE_SUMMARY="${PGGB_BASE}/vcf_reference_chromosome_communities.tsv"
OUTPUT_BASE="/Storage/data1/isabella.gallego/MAESTRIA/population_structure_24G_v2"

NAME_MAP="${CODE_DIR}/sample_name_map_24G.tsv"
METADATA="${CODE_DIR}/metadata_24G.tsv"
MANIFEST="${CODE_DIR}/community_manifest.tsv"
BUILD_SCRIPT="${CODE_DIR}/build_24G_snp_matrix.py"
ANALYSIS_SCRIPT="${CODE_DIR}/analyze_24G_population_structure.R"

MAX_MISSING="${MAX_MISSING:-0.20}"
PYTHON="${PYTHON:-python}"
RSCRIPT="${RSCRIPT:-Rscript}"

for f in "${REFERENCE_SUMMARY}" "${NAME_MAP}" "${METADATA}" "${MANIFEST}" \
         "${BUILD_SCRIPT}" "${ANALYSIS_SCRIPT}"; do
    [[ -f "$f" ]] || { echo "[ERROR] Missing required file: $f"; exit 1; }
done

mkdir -p "${OUTPUT_BASE}"

run_profile() {
    local profile="$1"
    local label="$2"
    local dir="${OUTPUT_BASE}/${profile}"
    mkdir -p "${dir}"

    echo "============================================================="
    echo "PROFILE: ${profile}"
    echo "============================================================="

    "${PYTHON}" "${BUILD_SCRIPT}" \
        --profile "${profile}" \
        --base-dir "${PGGB_BASE}" \
        --reference-summary "${REFERENCE_SUMMARY}" \
        --manifest "${MANIFEST}" \
        --name-map "${NAME_MAP}" \
        --metadata "${METADATA}" \
        --output-dir "${dir}" \
        --max-missing "${MAX_MISSING}"

    "${RSCRIPT}" "${ANALYSIS_SCRIPT}" \
        "${dir}/PGGB24_${profile}_snp_matrix_imputed.tsv" \
        "${METADATA}" \
        "${label}" \
        "${dir}/analysis"
}

echo "============================================================="
echo "PGGB 24-GENOME POPULATION-STRUCTURE PIPELINE V2"
echo "Start       : $(date '+%Y-%m-%d %H:%M:%S')"
echo "MAX_MISSING : ${MAX_MISSING}"
echo "Output      : ${OUTPUT_BASE}"
echo "============================================================="

run_profile "all_11" "PGGB24_all11"
run_profile "conservative_9" "PGGB24_conservative9"

echo "============================================================="
echo "PIPELINE COMPLETED"
echo "End     : $(date '+%Y-%m-%d %H:%M:%S')"
echo "Results : ${OUTPUT_BASE}"
echo "============================================================="
