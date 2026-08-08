#!/bin/bash

#$ -N VCF_generation_24G_chr
#$ -q all.q
#$ -cwd
#$ -pe smp 8
#$ -V
#$ -o /Storage/data1/isabella.gallego/MAESTRIA/code/PGGB_24_genomes/logs/vcf_generation_chr_$JOB_ID.out
#$ -e /Storage/data1/isabella.gallego/MAESTRIA/code/PGGB_24_genomes/logs/vcf_generation_chr_$JOB_ID.err

# =============================================================
# PGGB 24 genomes — VCF generation for chromosome communities
#
# Processes only the communities associated with the
# chromosome-scale syntenic groups.
#
# Group10 is represented across community.0 and community.1.
#
# For each selected community, the script:
#   1. Finds the most recent smooth.final.gfa file.
#   2. Converts the GFA graph to VG format.
#   3. Builds an XG index.
#   4. Lists the paths present in the graph.
#   5. Selects a valid reference prefix.
#   6. Runs vg deconstruct.
#   7. Compresses the VCF using BGZF.
#   8. Validates and indexes the VCF.
#
# Usage:
#   qsub VCF_generation_24G_chromosome_communities.sh
# =============================================================

set -euo pipefail

# =========================
# ENVIRONMENT
# =========================

module purge
module load singularity-ce/3.11.2
module load Samtools/1.22
module load bcftools/1.22

# =========================
# CONFIGURATION
# =========================

BASE_DIR="/Storage/data1/isabella.gallego/MAESTRIA/PGGB_results_24G/pggb_partitioned_results_24G"

SIF="/Storage/data1/isabella.gallego/MAESTRIA/pggb_latest.sif"

THREADS="${NSLOTS:-8}"

LOG_DIR="/Storage/data1/isabella.gallego/MAESTRIA/code/PGGB_24_genomes/logs"

REFERENCE_SUMMARY="${BASE_DIR}/vcf_reference_chromosome_communities.tsv"

# Communities associated with chromosome-scale syntenic groups.
COMMUNITIES=(
    community.0
    community.1
    community.2
    community.3
    community.4
    community.5
    community.6
    community.7
    community.8
    community.9
    community.12
)

# Preferred reference order.
#
# B3 is selected whenever it is present.
# community.0 does not contain B3, so C26 should be selected.
PREFERRED_REFERENCES=(
    "Mror_B3"
    "Mror_C26"
    "Mror_CO8"
    "Mror_CO84"
    "Mror_E7"
)

# =========================
# SETUP
# =========================

mkdir -p "${LOG_DIR}"

if [[ ! -d "${BASE_DIR}" ]]; then
    echo "[ERROR] BASE_DIR does not exist:"
    echo "        ${BASE_DIR}"
    exit 1
fi

if [[ ! -f "${SIF}" ]]; then
    echo "[ERROR] PGGB Singularity image not found:"
    echo "        ${SIF}"
    exit 1
fi

command -v singularity >/dev/null 2>&1 || {
    echo "[ERROR] singularity is not available in PATH."
    exit 1
}

command -v bgzip >/dev/null 2>&1 || {
    echo "[ERROR] bgzip is not available in PATH."
    exit 1
}

command -v tabix >/dev/null 2>&1 || {
    echo "[ERROR] tabix is not available in PATH."
    exit 1
}

command -v bcftools >/dev/null 2>&1 || {
    echo "[ERROR] bcftools is not available in PATH."
    exit 1
}

cd "${BASE_DIR}"

printf \
    "community\treference_prefix\treference_path\tgfa\tvariants\tsamples\tstatus\n" \
    > "${REFERENCE_SUMMARY}"

start_time="$(date '+%Y-%m-%d %H:%M:%S')"

echo "============================================================="
echo "PGGB 24-GENOME VCF GENERATION"
echo "Start             : ${start_time}"
echo "Job ID            : ${JOB_ID:-not_available}"
echo "Node              : $(hostname)"
echo "Threads           : ${THREADS}"
echo "BASE_DIR          : ${BASE_DIR}"
echo "REFERENCE_SUMMARY : ${REFERENCE_SUMMARY}"
echo "============================================================="

# =========================
# LOOP THROUGH COMMUNITIES
# =========================

for COMMUNITY in "${COMMUNITIES[@]}"; do

    COMMUNITY_DIR="all_genomes_pansn.fasta.bf3285f.${COMMUNITY}"

    FULL_COMMUNITY_DIR="${BASE_DIR}/${COMMUNITY_DIR}"

    echo ""
    echo "============================================================="
    echo "Processing: ${COMMUNITY}"
    echo "Directory : ${FULL_COMMUNITY_DIR}"
    echo "============================================================="

    # ---------------------------------------------------------
    # Validate community directory
    # ---------------------------------------------------------

    if [[ ! -d "${FULL_COMMUNITY_DIR}" ]]; then

        echo "[WARNING] Community directory not found."

        printf \
            "%s\tNA\tNA\tNA\t0\t0\tNO_DIRECTORY\n" \
            "${COMMUNITY}" \
            >> "${REFERENCE_SUMMARY}"

        continue
    fi

    cd "${FULL_COMMUNITY_DIR}"

    # ---------------------------------------------------------
    # 1. Select the most recent final GFA
    # ---------------------------------------------------------

    GFA="$(
        find . \
            -maxdepth 1 \
            -type f \
            -name "*.smooth.final.gfa" \
            -printf '%T@\t%p\n' \
        | sort -nr \
        | head -n 1 \
        | cut -f2-
    )"

    if [[ -z "${GFA}" || ! -f "${GFA}" ]]; then

        echo "[WARNING] No smooth.final.gfa file was found."

        printf \
            "%s\tNA\tNA\tNA\t0\t0\tNO_GFA\n" \
            "${COMMUNITY}" \
            >> "${REFERENCE_SUMMARY}"

        cd "${BASE_DIR}"
        continue
    fi

    echo "GFA: ${GFA}"

    # ---------------------------------------------------------
    # 2. Convert GFA to VG
    # ---------------------------------------------------------

    rm -f \
        graph.vg \
        graph.xg

    echo "Converting GFA to VG..."

    singularity exec \
        --bind /Storage \
        "${SIF}" \
        vg convert \
            -g "${GFA}" \
        > graph.vg

    if [[ ! -s graph.vg ]]; then

        echo "[ERROR] graph.vg is empty."

        printf \
            "%s\tNA\tNA\t%s\t0\t0\tEMPTY_VG\n" \
            "${COMMUNITY}" \
            "${GFA}" \
            >> "${REFERENCE_SUMMARY}"

        cd "${BASE_DIR}"
        continue
    fi

    # ---------------------------------------------------------
    # 3. Build XG index
    # ---------------------------------------------------------

    echo "Building XG index..."

    singularity exec \
        --bind /Storage \
        "${SIF}" \
        vg index \
            -x graph.xg \
            graph.vg

    if [[ ! -s graph.xg ]]; then

        echo "[ERROR] graph.xg is empty."

        printf \
            "%s\tNA\tNA\t%s\t0\t0\tEMPTY_XG\n" \
            "${COMMUNITY}" \
            "${GFA}" \
            >> "${REFERENCE_SUMMARY}"

        cd "${BASE_DIR}"
        continue
    fi

    # ---------------------------------------------------------
    # 4. Extract graph paths
    # ---------------------------------------------------------

    PATH_LIST="${COMMUNITY}_paths.txt"

    echo "Extracting graph paths..."

    singularity exec \
        --bind /Storage \
        "${SIF}" \
        vg paths \
            -x graph.xg \
            -L \
        > "${PATH_LIST}"

    if [[ ! -s "${PATH_LIST}" ]]; then

        echo "[ERROR] No graph paths were found."

        printf \
            "%s\tNA\tNA\t%s\t0\t0\tNO_PATHS\n" \
            "${COMMUNITY}" \
            "${GFA}" \
            >> "${REFERENCE_SUMMARY}"

        cd "${BASE_DIR}"
        continue
    fi

    N_PATHS="$(
        wc -l < "${PATH_LIST}"
    )"

    echo "Paths found: ${N_PATHS}"

    echo "Sample prefixes:"

    cut -d'#' -f1 "${PATH_LIST}" \
        | sort -u \
        | sed 's/^/  /'

    # ---------------------------------------------------------
    # 5. Select a reference present in the graph
    # ---------------------------------------------------------

    REF_PREFIX=""

    for candidate in "${PREFERRED_REFERENCES[@]}"; do

        if grep -q "^${candidate}#" "${PATH_LIST}"; then
            REF_PREFIX="${candidate}"
            break
        fi

    done

    # If none of the preferred PacBio genomes is present,
    # select the first sample prefix found in the graph.
    if [[ -z "${REF_PREFIX}" ]]; then

        REF_PREFIX="$(
            head -n 1 "${PATH_LIST}" \
            | cut -d'#' -f1
        )"

    fi

    REF_PATH="$(
        grep "^${REF_PREFIX}#" "${PATH_LIST}" \
        | head -n 1
    )"

    if [[ -z "${REF_PREFIX}" || -z "${REF_PATH}" ]]; then

        echo "[ERROR] A valid reference could not be selected."

        printf \
            "%s\tNA\tNA\t%s\t0\t0\tNO_REFERENCE\n" \
            "${COMMUNITY}" \
            "${GFA}" \
            >> "${REFERENCE_SUMMARY}"

        cd "${BASE_DIR}"
        continue
    fi

    echo "Reference prefix: ${REF_PREFIX}"
    echo "First reference path: ${REF_PATH}"

    # ---------------------------------------------------------
    # 6. Generate BGZF-compressed VCF
    # ---------------------------------------------------------

    rm -f \
        variants.vcf \
        variants.vcf.gz \
        variants.vcf.gz.tbi \
        variants.vcf.gz.csi

    DECONSTRUCT_ERROR_LOG="${COMMUNITY}_vg_deconstruct.err"

    rm -f "${DECONSTRUCT_ERROR_LOG}"

    echo "Running vg deconstruct..."

    set +e

    singularity exec \
        --bind /Storage \
        "${SIF}" \
        vg deconstruct \
            -P "${REF_PREFIX}" \
            -t "${THREADS}" \
            graph.xg \
        2> "${DECONSTRUCT_ERROR_LOG}" \
    | bgzip -c \
        > variants.vcf.gz

    # PIPESTATUS must be copied immediately because any later
    # shell command or assignment replaces its contents.
    pipeline_status=(
        "${PIPESTATUS[@]}"
    )

    deconstruct_status="${pipeline_status[0]:-1}"
    bgzip_status="${pipeline_status[1]:-1}"

    set -e

    echo "vg deconstruct exit code : ${deconstruct_status}"
    echo "bgzip exit code          : ${bgzip_status}"

    if [[ "${deconstruct_status}" -ne 0 ]]; then

        echo "[ERROR] vg deconstruct failed for ${COMMUNITY}."

        if [[ -s "${DECONSTRUCT_ERROR_LOG}" ]]; then
            echo "vg deconstruct error output:"
            cat "${DECONSTRUCT_ERROR_LOG}"
        fi

        rm -f \
            variants.vcf.gz \
            variants.vcf.gz.tbi \
            variants.vcf.gz.csi

        printf \
            "%s\t%s\t%s\t%s\t0\t0\tDECONSTRUCT_FAILED\n" \
            "${COMMUNITY}" \
            "${REF_PREFIX}" \
            "${REF_PATH}" \
            "${GFA}" \
            >> "${REFERENCE_SUMMARY}"

        cd "${BASE_DIR}"
        continue
    fi

    if [[ "${bgzip_status}" -ne 0 ]]; then

        echo "[ERROR] bgzip failed for ${COMMUNITY}."

        rm -f \
            variants.vcf.gz \
            variants.vcf.gz.tbi \
            variants.vcf.gz.csi

        printf \
            "%s\t%s\t%s\t%s\t0\t0\tBGZIP_FAILED\n" \
            "${COMMUNITY}" \
            "${REF_PREFIX}" \
            "${REF_PATH}" \
            "${GFA}" \
            >> "${REFERENCE_SUMMARY}"

        cd "${BASE_DIR}"
        continue
    fi

    if [[ ! -s variants.vcf.gz ]]; then

        echo "[ERROR] variants.vcf.gz is empty."

        rm -f \
            variants.vcf.gz \
            variants.vcf.gz.tbi \
            variants.vcf.gz.csi

        printf \
            "%s\t%s\t%s\t%s\t0\t0\tEMPTY_VCF\n" \
            "${COMMUNITY}" \
            "${REF_PREFIX}" \
            "${REF_PATH}" \
            "${GFA}" \
            >> "${REFERENCE_SUMMARY}"

        cd "${BASE_DIR}"
        continue
    fi

    # ---------------------------------------------------------
    # 7. Validate the VCF header
    # ---------------------------------------------------------

    if ! bcftools view \
        -h \
        variants.vcf.gz \
        | grep -q '^#CHROM'; then

        echo "[ERROR] Generated file is not a valid VCF."

        if [[ -s "${DECONSTRUCT_ERROR_LOG}" ]]; then
            echo "vg deconstruct stderr:"
            cat "${DECONSTRUCT_ERROR_LOG}"
        fi

        rm -f \
            variants.vcf.gz \
            variants.vcf.gz.tbi \
            variants.vcf.gz.csi

        printf \
            "%s\t%s\t%s\t%s\t0\t0\tINVALID_VCF\n" \
            "${COMMUNITY}" \
            "${REF_PREFIX}" \
            "${REF_PATH}" \
            "${GFA}" \
            >> "${REFERENCE_SUMMARY}"

        cd "${BASE_DIR}"
        continue
    fi

    # ---------------------------------------------------------
    # 8. Index VCF
    # ---------------------------------------------------------

    echo "Indexing VCF..."

    tabix \
        -f \
        -p vcf \
        variants.vcf.gz

    if [[ ! -f variants.vcf.gz.tbi ]]; then

        echo "[ERROR] VCF index was not created."

        printf \
            "%s\t%s\t%s\t%s\t0\t0\tINDEX_FAILED\n" \
            "${COMMUNITY}" \
            "${REF_PREFIX}" \
            "${REF_PATH}" \
            "${GFA}" \
            >> "${REFERENCE_SUMMARY}"

        cd "${BASE_DIR}"
        continue
    fi

    # ---------------------------------------------------------
    # 9. Count variants and samples
    # ---------------------------------------------------------

    N_VARIANTS="$(
        bcftools view \
            -H \
            variants.vcf.gz \
        | wc -l
    )"

    N_SAMPLES="$(
        bcftools query \
            -l \
            variants.vcf.gz \
        | wc -l
    )"

    echo ""
    echo "VCF generated successfully"
    echo "Community : ${COMMUNITY}"
    echo "Reference : ${REF_PREFIX}"
    echo "Variants  : ${N_VARIANTS}"
    echo "Samples   : ${N_SAMPLES}"
    echo "Output    : ${FULL_COMMUNITY_DIR}/variants.vcf.gz"

    printf \
        "%s\t%s\t%s\t%s\t%s\t%s\tSUCCESS\n" \
        "${COMMUNITY}" \
        "${REF_PREFIX}" \
        "${REF_PATH}" \
        "${GFA}" \
        "${N_VARIANTS}" \
        "${N_SAMPLES}" \
        >> "${REFERENCE_SUMMARY}"

    cd "${BASE_DIR}"

done

# =========================
# FINAL REPORT
# =========================

end_time="$(date '+%Y-%m-%d %H:%M:%S')"

SUCCESS_COUNT="$(
    awk \
        -F '\t' \
        'NR > 1 && $7 == "SUCCESS" {count++} END {print count + 0}' \
        "${REFERENCE_SUMMARY}"
)"

FAILURE_COUNT="$(
    awk \
        -F '\t' \
        'NR > 1 && $7 != "SUCCESS" {count++} END {print count + 0}' \
        "${REFERENCE_SUMMARY}"
)"

echo ""
echo "============================================================="
echo "SELECTED COMMUNITIES COMPLETED"
echo "Start             : ${start_time}"
echo "End               : ${end_time}"
echo "Successful VCFs   : ${SUCCESS_COUNT}"
echo "Failed VCFs       : ${FAILURE_COUNT}"
echo "Reference summary : ${REFERENCE_SUMMARY}"
echo "============================================================="

if [[ "${FAILURE_COUNT}" -gt 0 ]]; then
    exit 1
fi

exit 0
