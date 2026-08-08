#!/bin/bash

# =============================================================================
# run_final_geography_genetics_corrections.sh
#
# Applies all final corrections to the completed PGGB24
# geography-versus-genetics analysis.
#
# This workflow reuses existing distance matrices and robustness outputs.
# It does not rebuild the PGGB SNP matrices.
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# 1. Paths
# -----------------------------------------------------------------------------

CODE_DIR="/Storage/data1/isabella.gallego/MAESTRIA/code/population_structure/geography_genetics/final_corrections"

GEOGRAPHY_BASE="/Storage/data1/isabella.gallego/MAESTRIA/geography_genetics_24G"

POPULATION_BASE="/Storage/data1/isabella.gallego/MAESTRIA/population_structure_24G_v2"

POPULATION_CODE="/Storage/data1/isabella.gallego/MAESTRIA/code/population_structure/population_structure_24G_v2"

OUTPUT_DIR="${GEOGRAPHY_BASE}/final_corrected"

COORDINATES="${GEOGRAPHY_BASE}/01_coordinate_validation/sample_coordinates_24G_validated.tsv"

METADATA="${POPULATION_CODE}/metadata_24G.tsv"

MANIFEST="${POPULATION_CODE}/community_manifest.tsv"

ALL11_DISTANCE="${POPULATION_BASE}/all_11/analysis/all_24_genomes/PGGB24_all11_all_24_genomes_distance_matrix.tsv"

CONSERVATIVE9_DISTANCE="${POPULATION_BASE}/conservative_9/analysis/all_24_genomes/PGGB24_conservative9_all_24_genomes_distance_matrix.tsv"

ALL11_VCF_QC="${POPULATION_BASE}/all_11/PGGB24_all_11_vcf_qc.tsv"

CONSERVATIVE9_VCF_QC="${POPULATION_BASE}/conservative_9/PGGB24_conservative_9_vcf_qc.tsv"

ALL11_FEATURES="${POPULATION_BASE}/all_11/PGGB24_all_11_snp_feature_metadata.tsv"

CONSERVATIVE9_FEATURES="${POPULATION_BASE}/conservative_9/PGGB24_conservative_9_snp_feature_metadata.tsv"

ALL11_MATRIX_SUMMARY="${POPULATION_BASE}/all_11/PGGB24_all_11_matrix_summary.tsv"

CONSERVATIVE9_MATRIX_SUMMARY="${POPULATION_BASE}/conservative_9/PGGB24_conservative_9_matrix_summary.tsv"

ALL11_ROBUSTNESS="${GEOGRAPHY_BASE}/04_balanced_robustness/all11/all11_balanced_ibd_summary.tsv"

CONSERVATIVE9_ROBUSTNESS="${GEOGRAPHY_BASE}/04_balanced_robustness/conservative9/conservative9_balanced_ibd_summary.tsv"

COORDINATE_SUMMARY="${GEOGRAPHY_BASE}/01_coordinate_validation/coordinate_validation_summary.tsv"

STATISTICS_SCRIPT="${CODE_DIR}/recalculate_spatial_statistics_final.R"

FIGURE_SCRIPT="${CODE_DIR}/make_spatial_genetic_figures_final.R"

COMMUNITY_AUDIT_SCRIPT="${CODE_DIR}/audit_community_profiles.py"

MASTER_SUMMARY_SCRIPT="${CODE_DIR}/build_pipeline_master_summary.py"

STATISTICS_DIR="${OUTPUT_DIR}/01_corrected_statistics"

FIGURE_DIR="${OUTPUT_DIR}/02_corrected_figures"

AUDIT_DIR="${OUTPUT_DIR}/03_community_audit"

SUMMARY_DIR="${OUTPUT_DIR}/04_master_summary"

mkdir -p \
    "${STATISTICS_DIR}" \
    "${FIGURE_DIR}" \
    "${AUDIT_DIR}" \
    "${SUMMARY_DIR}"

# -----------------------------------------------------------------------------
# 2. Input checks
# -----------------------------------------------------------------------------

required_files=(
    "${COORDINATES}"
    "${METADATA}"
    "${MANIFEST}"
    "${ALL11_DISTANCE}"
    "${CONSERVATIVE9_DISTANCE}"
    "${ALL11_VCF_QC}"
    "${CONSERVATIVE9_VCF_QC}"
    "${ALL11_FEATURES}"
    "${CONSERVATIVE9_FEATURES}"
    "${ALL11_MATRIX_SUMMARY}"
    "${CONSERVATIVE9_MATRIX_SUMMARY}"
    "${ALL11_ROBUSTNESS}"
    "${CONSERVATIVE9_ROBUSTNESS}"
    "${COORDINATE_SUMMARY}"
    "${STATISTICS_SCRIPT}"
    "${FIGURE_SCRIPT}"
    "${COMMUNITY_AUDIT_SCRIPT}"
    "${MASTER_SUMMARY_SCRIPT}"
)

for file_path in "${required_files[@]}"; do
    if [[ ! -f "${file_path}" ]]; then
        echo "[ERROR] Required file not found:"
        echo "        ${file_path}"
        exit 1
    fi
done

echo "======================================================================"
echo "FINAL PGGB24 GEOGRAPHY-GENETICS CORRECTIONS"
echo "Start      : $(date '+%Y-%m-%d %H:%M:%S')"
echo "Output     : ${OUTPUT_DIR}"
echo "======================================================================"

# -----------------------------------------------------------------------------
# 3. Recalculate corrected statistics
# -----------------------------------------------------------------------------

echo ""
echo "[1/4] Recalculating corrected spatial statistics..."

Rscript "${STATISTICS_SCRIPT}" \
    "${COORDINATES}" \
    "${METADATA}" \
    "${ALL11_DISTANCE}" \
    "${CONSERVATIVE9_DISTANCE}" \
    "${STATISTICS_DIR}"

# -----------------------------------------------------------------------------
# 4. Generate corrected figures
# -----------------------------------------------------------------------------

echo ""
echo "[2/4] Generating corrected figures..."

Rscript "${FIGURE_SCRIPT}" \
    "${COORDINATES}" \
    "${METADATA}" \
    "${STATISTICS_DIR}/geographic_genetic_pair_table_final.tsv" \
    "${STATISTICS_DIR}/mantel_results_final.tsv" \
    "${STATISTICS_DIR}/mrm_coefficients_final.tsv" \
    "${FIGURE_DIR}"

# -----------------------------------------------------------------------------
# 5. Audit community labels and retained contributions
# -----------------------------------------------------------------------------

echo ""
echo "[3/4] Auditing community profiles..."

python "${COMMUNITY_AUDIT_SCRIPT}" \
    --manifest "${MANIFEST}" \
    --all11-vcf-qc "${ALL11_VCF_QC}" \
    --all11-features "${ALL11_FEATURES}" \
    --conservative9-vcf-qc "${CONSERVATIVE9_VCF_QC}" \
    --conservative9-features "${CONSERVATIVE9_FEATURES}" \
    --output-dir "${AUDIT_DIR}"

# -----------------------------------------------------------------------------
# 6. Build one master summary
# -----------------------------------------------------------------------------

echo ""
echo "[4/4] Building master summary..."

python "${MASTER_SUMMARY_SCRIPT}" \
    --coordinate-summary "${COORDINATE_SUMMARY}" \
    --pair-table "${STATISTICS_DIR}/geographic_genetic_pair_table_final.tsv" \
    --mantel "${STATISTICS_DIR}/mantel_results_final.tsv" \
    --mrm-coefficients "${STATISTICS_DIR}/mrm_coefficients_final.tsv" \
    --mrm-models "${STATISTICS_DIR}/mrm_model_summary_final.tsv" \
    --regression "${STATISTICS_DIR}/descriptive_regression_results_final.tsv" \
    --all11-robustness "${ALL11_ROBUSTNESS}" \
    --conservative9-robustness "${CONSERVATIVE9_ROBUSTNESS}" \
    --all11-matrix-summary "${ALL11_MATRIX_SUMMARY}" \
    --conservative9-matrix-summary "${CONSERVATIVE9_MATRIX_SUMMARY}" \
    --community-audit "${AUDIT_DIR}/community_profile_audit.tsv" \
    --colocation-summary "${STATISTICS_DIR}/colocation_representative_mantel_summary.tsv" \
    --output-dir "${SUMMARY_DIR}"

echo ""
echo "======================================================================"
echo "CORRECTIONS COMPLETED"
echo "End        : $(date '+%Y-%m-%d %H:%M:%S')"
echo "Results    : ${OUTPUT_DIR}"
echo "======================================================================"
