#!/bin/bash

# =============================================================
# run_population_structure_pipeline.sh
#
# Unified population-structure analysis for PGGB and
# chromosome-level minigraph-cactus VCFs.
#
# Verified VCF samples:
#   PGGB:   C26, CO8, CO84, E7
#   Cactus: Mror_SAMPLE_GroupN
#
# Verified reference coordinates:
#   PGGB:   B3#0#GroupN
#   Cactus: Mror_B3_GroupN
#
# B3 is reconstructed as 0 for every retained SNP ALT allele.
# =============================================================

set -euo pipefail

# =========================
# CONFIGURATION
# =========================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
METADATA_FILE="${SCRIPT_DIR}/sample_metadata.tsv"

PGGB_BASE="/Storage/data1/isabella.gallego/MAESTRIA/PGGB_results_pacbio/pggb_partitioned_results"
CACTUS_BASE="/Storage/data1/isabella.gallego/MAESTRIA/cactus_chroms"
OUTPUT_BASE="/Storage/data1/isabella.gallego/MAESTRIA/population_structure"

PYTHON="${PYTHON:-python}"
RSCRIPT="${RSCRIPT:-Rscript}"
REFERENCE_SAMPLE="B3"
MAX_MISSING="0.20"

PGGB_OUTPUT="${OUTPUT_BASE}/PGGB"
CACTUS_OUTPUT="${OUTPUT_BASE}/Cactus"
COMPARISON_OUTPUT="${OUTPUT_BASE}/comparison"
MANIFEST_DIR="${OUTPUT_BASE}/manifests"

mkdir -p "${PGGB_OUTPUT}" "${CACTUS_OUTPUT}" "${COMPARISON_OUTPUT}" "${MANIFEST_DIR}"

PGGB_VCF_LIST="${MANIFEST_DIR}/pggb_vcfs.txt"
CACTUS_VCF_LIST="${MANIFEST_DIR}/cactus_vcfs.txt"

# =========================
# PRE-RUN CHECKS
# =========================

for required_file in \
    "${METADATA_FILE}" \
    "${SCRIPT_DIR}/build_snp_matrix.py" \
    "${SCRIPT_DIR}/analyze_population_structure.R" \
    "${SCRIPT_DIR}/compare_pggb_cactus_pca.R"
do
    [[ -f "${required_file}" ]] || {
        echo "[ERROR] Required file not found: ${required_file}"
        exit 1
    }
done

[[ -d "${PGGB_BASE}" ]] || { echo "[ERROR] PGGB_BASE does not exist: ${PGGB_BASE}"; exit 1; }
[[ -d "${CACTUS_BASE}" ]] || { echo "[ERROR] CACTUS_BASE does not exist: ${CACTUS_BASE}"; exit 1; }

command -v "${PYTHON}" >/dev/null 2>&1 || { echo "[ERROR] Python not found: ${PYTHON}"; exit 1; }
command -v "${RSCRIPT}" >/dev/null 2>&1 || { echo "[ERROR] Rscript not found: ${RSCRIPT}"; exit 1; }

# =========================
# DISCOVER VCF FILES
# =========================

echo "[1/5] Discovering VCF files..."

find "${PGGB_BASE}" -maxdepth 2 -type f \
    \( -name "variants.vcf" -o -name "variants.vcf.gz" \) \
    | sort -V > "${PGGB_VCF_LIST}"

find "${CACTUS_BASE}" -maxdepth 2 -type f \
    -name "cactus_group*_v1.vcf.gz" \
    | sort -V > "${CACTUS_VCF_LIST}"

PGGB_VCF_COUNT="$(wc -l < "${PGGB_VCF_LIST}")"
CACTUS_VCF_COUNT="$(wc -l < "${CACTUS_VCF_LIST}")"

[[ "${PGGB_VCF_COUNT}" -gt 0 ]] || { echo "[ERROR] No PGGB VCFs found."; exit 1; }
[[ "${CACTUS_VCF_COUNT}" -gt 0 ]] || { echo "[ERROR] No Cactus VCFs found."; exit 1; }

echo "  PGGB VCFs   : ${PGGB_VCF_COUNT}"
echo "  Cactus VCFs : ${CACTUS_VCF_COUNT}"

# =========================
# BUILD SNP MATRICES
# =========================

echo ""
echo "[2/5] Building PGGB SNP matrix..."

"${PYTHON}" "${SCRIPT_DIR}/build_snp_matrix.py" \
    --source "PGGB" \
    --vcf-list "${PGGB_VCF_LIST}" \
    --metadata "${METADATA_FILE}" \
    --sample-name-mode "pggb" \
    --reference-sample "${REFERENCE_SAMPLE}" \
    --output-prefix "${PGGB_OUTPUT}/PGGB" \
    --max-missing "${MAX_MISSING}"

echo ""
echo "[3/5] Building Cactus SNP matrix..."

"${PYTHON}" "${SCRIPT_DIR}/build_snp_matrix.py" \
    --source "Cactus" \
    --vcf-list "${CACTUS_VCF_LIST}" \
    --metadata "${METADATA_FILE}" \
    --sample-name-mode "cactus" \
    --reference-sample "${REFERENCE_SAMPLE}" \
    --output-prefix "${CACTUS_OUTPUT}/Cactus" \
    --max-missing "${MAX_MISSING}"

# =========================
# POPULATION ANALYSES
# =========================

echo ""
echo "[4/5] Running PCA, distances, and clustering..."

"${RSCRIPT}" "${SCRIPT_DIR}/analyze_population_structure.R" \
    "${PGGB_OUTPUT}/PGGB_snp_matrix_imputed.tsv" \
    "${METADATA_FILE}" \
    "PGGB" \
    "${PGGB_OUTPUT}"

"${RSCRIPT}" "${SCRIPT_DIR}/analyze_population_structure.R" \
    "${CACTUS_OUTPUT}/Cactus_snp_matrix_imputed.tsv" \
    "${METADATA_FILE}" \
    "Cactus" \
    "${CACTUS_OUTPUT}"

# =========================
# METHOD COMPARISON
# =========================

echo ""
echo "[5/5] Generating PGGB-versus-Cactus PCA panel..."

"${RSCRIPT}" "${SCRIPT_DIR}/compare_pggb_cactus_pca.R" \
    "${PGGB_OUTPUT}/PGGB_pca_coordinates.tsv" \
    "${CACTUS_OUTPUT}/Cactus_pca_coordinates.tsv" \
    "${COMPARISON_OUTPUT}/PGGB_vs_Cactus_PCA.pdf"

echo ""
echo "============================================================="
echo "POPULATION-STRUCTURE PIPELINE COMPLETED"
echo "Reference reconstructed: ${REFERENCE_SAMPLE}"
echo "PGGB results          : ${PGGB_OUTPUT}"
echo "Cactus results        : ${CACTUS_OUTPUT}"
echo "Comparison results    : ${COMPARISON_OUTPUT}"
echo "============================================================="
