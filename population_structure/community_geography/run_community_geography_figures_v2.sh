#!/bin/bash

# =============================================================================
# run_community_geography_figures_v2.sh
#
# Regenerates publication-ready community/chromosome geography figures for
# all11 and conservative9 without recalculating any statistics.
# =============================================================================

set -euo pipefail

CODE_DIR="/Storage/data1/isabella.gallego/MAESTRIA/code/population_structure/community_geography"

RESULTS_BASE="/Storage/data1/isabella.gallego/MAESTRIA/community_geography_24G"

MAP_FILE="${CODE_DIR}/community_chromosome_map.tsv"

PLOT_SCRIPT="${CODE_DIR}/plot_community_geographic_contribution_v2.R"

PROFILES="${PROFILES:-all11 conservative9}"

for required_file in \
    "${MAP_FILE}" \
    "${PLOT_SCRIPT}"
do
    if [[ ! -f "${required_file}" ]]; then
        echo "[ERROR] Required file not found:"
        echo "        ${required_file}"
        exit 1
    fi
done

echo "======================================================================"
echo "COMMUNITY GEOGRAPHY FIGURE PIPELINE V2"
echo "Start    : $(date '+%Y-%m-%d %H:%M:%S')"
echo "Profiles : ${PROFILES}"
echo "======================================================================"

for profile in ${PROFILES}; do

    ranking_file="${RESULTS_BASE}/${profile}/02_ranking/${profile}_community_geographic_ranking.tsv"

    output_dir="${RESULTS_BASE}/${profile}/06_publication_figures_v2"

    if [[ ! -f "${ranking_file}" ]]; then
        echo "[ERROR] Ranking file not found for ${profile}:"
        echo "        ${ranking_file}"
        exit 1
    fi

    mkdir -p "${output_dir}"

    echo ""
    echo "Processing profile: ${profile}"

    Rscript "${PLOT_SCRIPT}" \
        "${ranking_file}" \
        "${MAP_FILE}" \
        "${output_dir}" \
        "${profile}"

    echo "Completed: ${output_dir}"
done

echo ""
echo "======================================================================"
echo "FIGURE PIPELINE COMPLETED"
echo "End: $(date '+%Y-%m-%d %H:%M:%S')"
echo "======================================================================"
