#!/bin/bash

# =============================================================
# run_geography_genetics_pipeline.sh
#
# Complete PGGB24 geography-versus-genetics workflow.
#
# Before running:
#   1. Copy sample_coordinates_24G_template.tsv to
#      sample_coordinates_24G.tsv
#   2. Fill latitude, longitude, precision, and source.
#   3. Review every coordinate manually.
# =============================================================

set -euo pipefail

CODE_DIR="/Storage/data1/isabella.gallego/MAESTRIA/code/population_structure/geography_genetics"
POP_DIR="/Storage/data1/isabella.gallego/MAESTRIA/population_structure_24G_v2"
OUTPUT_DIR="/Storage/data1/isabella.gallego/MAESTRIA/geography_genetics_24G"

METADATA="/Storage/data1/isabella.gallego/MAESTRIA/code/population_structure/population_structure_24G_v2/metadata_24G.tsv"
COORDINATES="${CODE_DIR}/sample_coordinates_24G.tsv"

ALL11_DISTANCE="${POP_DIR}/all_11/analysis/all_24_genomes/PGGB24_all11_all_24_genomes_distance_matrix.tsv"
CONS9_DISTANCE="${POP_DIR}/conservative_9/analysis/all_24_genomes/PGGB24_conservative9_all_24_genomes_distance_matrix.tsv"

ALL11_MATRIX="${POP_DIR}/all_11/PGGB24_all_11_snp_matrix_imputed.tsv"
ALL11_FEATURES="${POP_DIR}/all_11/PGGB24_all_11_snp_feature_metadata.tsv"
CONS9_MATRIX="${POP_DIR}/conservative_9/PGGB24_conservative_9_snp_matrix_imputed.tsv"
CONS9_FEATURES="${POP_DIR}/conservative_9/PGGB24_conservative_9_snp_feature_metadata.tsv"

VALIDATE_SCRIPT="${CODE_DIR}/validate_coordinates.py"
STATS_SCRIPT="${CODE_DIR}/build_and_analyze_spatial_genetics.R"
PLOT_SCRIPT="${CODE_DIR}/make_spatial_genetic_figures.R"
ROBUSTNESS_SCRIPT="${CODE_DIR}/run_balanced_ibd_robustness.R"

REPLICATES="${REPLICATES:-100}"
BASE_SEED="${BASE_SEED:-20260803}"
N_PER_COMMUNITY="${N_PER_COMMUNITY:-0}"

VALIDATION_DIR="${OUTPUT_DIR}/01_coordinate_validation"
STATS_DIR="${OUTPUT_DIR}/02_spatial_statistics"
FIGURE_DIR="${OUTPUT_DIR}/03_figures"
ROBUSTNESS_DIR="${OUTPUT_DIR}/04_balanced_robustness"

mkdir -p "${VALIDATION_DIR}" "${STATS_DIR}" "${FIGURE_DIR}" "${ROBUSTNESS_DIR}"

for file in \
    "${METADATA}" "${COORDINATES}" \
    "${ALL11_DISTANCE}" "${CONS9_DISTANCE}" \
    "${ALL11_MATRIX}" "${ALL11_FEATURES}" \
    "${CONS9_MATRIX}" "${CONS9_FEATURES}" \
    "${VALIDATE_SCRIPT}" "${STATS_SCRIPT}" \
    "${PLOT_SCRIPT}" "${ROBUSTNESS_SCRIPT}"
do
    [[ -f "${file}" ]] || {
        echo "[ERROR] Required file not found:"
        echo "        ${file}"
        exit 1
    }
done

echo "============================================================="
echo "PGGB24 GEOGRAPHY-VERSUS-GENETICS PIPELINE"
echo "Start          : $(date '+%Y-%m-%d %H:%M:%S')"
echo "Coordinates    : ${COORDINATES}"
echo "Replicates     : ${REPLICATES}"
echo "Base seed      : ${BASE_SEED}"
echo "SNPs/community : ${N_PER_COMMUNITY} (0 = automatic minimum)"
echo "Output         : ${OUTPUT_DIR}"
echo "============================================================="

# -------------------------------------------------------------
# 1. Validate coordinate table
# -------------------------------------------------------------

python "${VALIDATE_SCRIPT}" \
    --coordinates "${COORDINATES}" \
    --metadata "${METADATA}" \
    --output-dir "${VALIDATION_DIR}"

VALID_COORDS="${VALIDATION_DIR}/sample_coordinates_24G_validated.tsv"

# -------------------------------------------------------------
# 2. Build geographic matrices and run spatial statistics
# -------------------------------------------------------------

Rscript "${STATS_SCRIPT}" \
    "${VALID_COORDS}" \
    "${METADATA}" \
    "${ALL11_DISTANCE}" \
    "${CONS9_DISTANCE}" \
    "${STATS_DIR}"

# -------------------------------------------------------------
# 3. Create maps and figures
# -------------------------------------------------------------

Rscript "${PLOT_SCRIPT}" \
    "${VALID_COORDS}" \
    "${METADATA}" \
    "${STATS_DIR}/geographic_genetic_pair_table.tsv" \
    "${FIGURE_DIR}"

# -------------------------------------------------------------
# 4. Balanced community-level robustness
# -------------------------------------------------------------

Rscript "${ROBUSTNESS_SCRIPT}" \
    "${ALL11_MATRIX}" \
    "${ALL11_FEATURES}" \
    "${STATS_DIR}/geographic_distance_matrix_km.tsv" \
    "all11" \
    "${ROBUSTNESS_DIR}/all11" \
    "${REPLICATES}" \
    "${BASE_SEED}" \
    "${N_PER_COMMUNITY}"

Rscript "${ROBUSTNESS_SCRIPT}" \
    "${CONS9_MATRIX}" \
    "${CONS9_FEATURES}" \
    "${STATS_DIR}/geographic_distance_matrix_km.tsv" \
    "conservative9" \
    "${ROBUSTNESS_DIR}/conservative9" \
    "${REPLICATES}" \
    "$((BASE_SEED + REPLICATES))" \
    "${N_PER_COMMUNITY}"

echo ""
echo "============================================================="
echo "PIPELINE COMPLETED"
echo "End     : $(date '+%Y-%m-%d %H:%M:%S')"
echo "Results : ${OUTPUT_DIR}"
echo "============================================================="
