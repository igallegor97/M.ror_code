#!/bin/bash

# =============================================================================
# run_community_geography_pipeline.sh
#
# Complete community-level genomic geography pipeline.
# =============================================================================

set -euo pipefail

CODE_DIR="/Storage/data1/isabella.gallego/MAESTRIA/code/population_structure/community_geography"

POP_BASE="/Storage/data1/isabella.gallego/MAESTRIA/population_structure_24G_v2"
GEO_BASE="/Storage/data1/isabella.gallego/MAESTRIA/geography_genetics_24G"
POP_CODE="/Storage/data1/isabella.gallego/MAESTRIA/code/population_structure/population_structure_24G_v2"

OUTPUT_DIR="/Storage/data1/isabella.gallego/MAESTRIA/community_geography_24G"

PROFILE="${PROFILE:-all11}"
MIN_SNPS="${MIN_SNPS:-100}"
PERMUTATIONS="${PERMUTATIONS:-9999}"

MATRIX="${POP_BASE}/all_11/PGGB24_all_11_snp_matrix_imputed.tsv"
FEATURES="${POP_BASE}/all_11/PGGB24_all_11_snp_feature_metadata.tsv"

if [[ "${PROFILE}" == "conservative9" ]]; then
    MATRIX="${POP_BASE}/conservative_9/PGGB24_conservative_9_snp_matrix_imputed.tsv"
    FEATURES="${POP_BASE}/conservative_9/PGGB24_conservative_9_snp_feature_metadata.tsv"
fi

GEO_MATRIX="${GEO_BASE}/final_corrected/01_corrected_statistics/geographic_distance_matrix_km_final.tsv"
COORDINATES="${GEO_BASE}/01_coordinate_validation/sample_coordinates_24G_validated.tsv"
METADATA="${POP_CODE}/metadata_24G.tsv"

COMMUNITY_MAP="${CODE_DIR}/community_chromosome_map.tsv"
GFF3="${GFF3:-/Storage/data1/isabella.gallego/MAESTRIA/data/genomes/Mror_C26.gff3}"
FUNCTIONAL_ANNOTATIONS="${FUNCTIONAL_ANNOTATIONS:-${CODE_DIR}/gene_function_annotations.tsv}"

PROFILE_DIR="${OUTPUT_DIR}/${PROFILE}"
ANALYSIS_DIR="${PROFILE_DIR}/01_per_community_analysis"
RANKING_DIR="${PROFILE_DIR}/02_ranking"
MAP_DIR="${PROFILE_DIR}/03_mapping"
FIGURE_DIR="${PROFILE_DIR}/04_figures"
BIOLOGY_DIR="${PROFILE_DIR}/05_biological_characterization"

mkdir -p \
    "${ANALYSIS_DIR}" \
    "${RANKING_DIR}" \
    "${MAP_DIR}" \
    "${FIGURE_DIR}" \
    "${BIOLOGY_DIR}"

required_files=(
    "${MATRIX}"
    "${FEATURES}"
    "${GEO_MATRIX}"
    "${COORDINATES}"
    "${METADATA}"
    "${COMMUNITY_MAP}"
    "${GFF3}"
)

for path in "${required_files[@]}"; do
    [[ -f "${path}" ]] || {
        echo "[ERROR] Required file not found: ${path}"
        exit 1
    }
done

echo "======================================================================"
echo "COMMUNITY-LEVEL GENOMIC GEOGRAPHY PIPELINE"
echo "Profile      : ${PROFILE}"
echo "Minimum SNPs : ${MIN_SNPS}"
echo "Permutations : ${PERMUTATIONS}"
echo "Output       : ${PROFILE_DIR}"
echo "======================================================================"

Rscript "${CODE_DIR}/per_community_spatial_analysis.R" \
    "${MATRIX}" \
    "${FEATURES}" \
    "${GEO_MATRIX}" \
    "${COORDINATES}" \
    "${METADATA}" \
    "${PROFILE}" \
    "${ANALYSIS_DIR}" \
    "${MIN_SNPS}" \
    "${PERMUTATIONS}"

python "${CODE_DIR}/rank_communities.py" \
    --mantel "${ANALYSIS_DIR}/${PROFILE}_community_mantel.tsv" \
    --mrm-coefficients "${ANALYSIS_DIR}/${PROFILE}_community_mrm_coefficients.tsv" \
    --mrm-models "${ANALYSIS_DIR}/${PROFILE}_community_mrm_models.tsv" \
    --regression "${ANALYSIS_DIR}/${PROFILE}_community_regression.tsv" \
    --community-summary "${ANALYSIS_DIR}/${PROFILE}_community_analysis_summary.tsv" \
    --community-map "${COMMUNITY_MAP}" \
    --profile "${PROFILE}" \
    --subset "high_precision" \
    --output-dir "${RANKING_DIR}"

python "${CODE_DIR}/validate_community_chromosome_map.py" \
    --map "${COMMUNITY_MAP}" \
    --ranking "${RANKING_DIR}/${PROFILE}_community_geographic_ranking.tsv" \
    --output-dir "${MAP_DIR}"

Rscript "${CODE_DIR}/plot_community_geographic_contribution.R" \
    "${RANKING_DIR}/${PROFILE}_community_geographic_ranking.tsv" \
    "${MAP_DIR}/community_chromosome_map_validated.tsv" \
    "${FIGURE_DIR}" \
    "${PROFILE}"

python "${CODE_DIR}/annotate_community_snps.py" \
    --features "${FEATURES}" \
    --gff3 "${GFF3}" \
    --reference-prefix "Mror_C26" \
    --output-dir "${BIOLOGY_DIR}"

if [[ -f "${FUNCTIONAL_ANNOTATIONS}" ]]; then
    python "${CODE_DIR}/enrich_top_community_genes.py" \
        --gene-summary "${BIOLOGY_DIR}/community_gene_summary.tsv" \
        --annotations "${FUNCTIONAL_ANNOTATIONS}" \
        --top-communities "${RANKING_DIR}/${PROFILE}_top_geographic_communities.tsv" \
        --output-dir "${BIOLOGY_DIR}/functional_enrichment"
else
    echo "[WARNING] Functional annotation table not found:"
    echo "          ${FUNCTIONAL_ANNOTATIONS}"
    echo "          SNP-to-gene annotation was completed, but GO/PFAM/CAZyme"
    echo "          enrichment was skipped."
fi

echo "======================================================================"
echo "PIPELINE COMPLETED"
echo "Results: ${PROFILE_DIR}"
echo "======================================================================"
