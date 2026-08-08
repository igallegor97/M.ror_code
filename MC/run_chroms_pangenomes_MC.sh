```bash
#!/bin/bash

# =============================================================================
# run_pangenomes.sh
#
# Orchestrator script that submits two chained SGE jobs for each chromosome
# group using -hold_jid.
#
# Job 1: prep_groupN
#   Downloads the input FASTA, splits it by genome, and simplifies the headers.
#
# Job 2: cactus_groupN
#   Runs cactus-pangenome after the corresponding preparation job finishes.
#
# Usage:
#   bash run_pangenomes.sh [seqfile] [output_base]
#
# Arguments:
#   seqfile
#       Tab-separated file containing:
#         group_name    fasta_url
#
#   output_base
#       Base directory where FASTA files, sequence files, temporary files,
#       and Cactus results will be stored.
# =============================================================================

set -euo pipefail

# =========================
# COMMAND-LINE ARGUMENTS
# =========================

SEQFILE="${1:-/Storage/data1/isabella.gallego/MAESTRIA/cactus_chroms/seqfile_chroms.txt}"

OUTPUT_BASE="${2:-/Storage/data1/isabella.gallego/MAESTRIA/cactus_chroms}"

# =========================
# FIXED PATHS
# =========================

CODE_DIR="/Storage/data1/isabella.gallego/MAESTRIA/code/minigraph_cactus/chroms_cactus"

FASTA_DIR="${OUTPUT_BASE}/fastas_clean"

SEQFILES_DIR="${OUTPUT_BASE}/seqfiles"

LOGS_DIR="${CODE_DIR}/logs"

WORK_DIR="${OUTPUT_BASE}/tempdir"

mkdir -p \
    "${FASTA_DIR}" \
    "${SEQFILES_DIR}" \
    "${LOGS_DIR}" \
    "${WORK_DIR}"

# =========================
# FUNCTIONS
# =========================

# -----------------------------------------------------------------------------
# Submit the preparation job
#
# This job:
#   1. Downloads the group FASTA file.
#   2. Splits the combined FASTA into one file per genome.
#   3. Removes the initial sample prefix from each FASTA header.
#   4. Creates a Cactus seqfile for the group.
# -----------------------------------------------------------------------------

submit_preparation_job() {
    local group="$1"
    local url="$2"

    local raw_fasta="${FASTA_DIR}/${group}_raw.fasta"
    local group_output_dir="${FASTA_DIR}/${group}"
    local group_seqfile="${SEQFILES_DIR}/seqfile_${group}.txt"
    local job_script="${CODE_DIR}/prep_${group}.sh"

    cat > "${job_script}" << EOF
#!/bin/bash

#$ -N prep_${group}
#$ -q all.q
#$ -cwd
#$ -pe smp 2
#$ -V
#$ -o ${LOGS_DIR}/prep_${group}_\$JOB_ID.out
#$ -e ${LOGS_DIR}/prep_${group}_\$JOB_ID.err

# =============================================================
# Preparation job for ${group}
#
# Downloads the combined FASTA file, separates genomes, simplifies
# FASTA headers, and creates the corresponding Cactus seqfile.
# =============================================================

set -euo pipefail

# =========================
# CONFIGURATION
# =========================

RAW_FASTA="${raw_fasta}"

GROUP_OUTPUT_DIR="${group_output_dir}"

GROUP_SEQFILE="${group_seqfile}"

DOWNLOAD_URL="${url}"

# =========================
# SETUP
# =========================

mkdir -p "\${GROUP_OUTPUT_DIR}"

echo "============================================================="
echo "PREPARATION JOB — ${group}"
echo "Start            : \$(date '+%Y-%m-%d %H:%M:%S')"
echo "Job ID           : \${JOB_ID:-not_available}"
echo "Node             : \$(hostname)"
echo "RAW_FASTA        : \${RAW_FASTA}"
echo "GROUP_OUTPUT_DIR : \${GROUP_OUTPUT_DIR}"
echo "GROUP_SEQFILE    : \${GROUP_SEQFILE}"
echo "DOWNLOAD_URL     : \${DOWNLOAD_URL}"
echo "============================================================="
echo ""

# =========================
# DOWNLOAD FASTA
# =========================

if [[ ! -f "\${RAW_FASTA}" ]]; then
    echo "[prep_${group}] Downloading FASTA file..."

    wget \
        -q \
        -O "\${RAW_FASTA}" \
        "\${DOWNLOAD_URL}"

    echo "[prep_${group}] Download completed."
else
    echo "[prep_${group}] FASTA file already exists. Skipping download."
fi

# =========================
# SPLIT GENOMES AND CLEAN HEADERS
# =========================

# Example header conversion:
#
#   >b3_Mror_B3_Group1
#
# becomes:
#
#   >Mror_B3_Group1

echo "[prep_${group}] Splitting genomes and simplifying headers..."

# Create or clear the group seqfile.
: > "\${GROUP_SEQFILE}"

# Write the AWK script using printf to prevent Bash variable expansion.
printf '%s\n' \
    '/^>/ {' \
    '    raw = substr(DOLLAR0, 2)' \
    '    sub(/^[^_]*_/, "", raw)' \
    '    clean = raw' \
    '    current_file = outdir "/" clean ".fasta"' \
    '    print ">" clean > current_file' \
    '    print clean "\t" current_file >> seqfile' \
    '}' \
    '!/^>/ {' \
    '    if (current_file != "")' \
    '        print >> current_file' \
    '}' \
    | sed 's/DOLLAR0/\$0/g' \
    > "/tmp/split_fasta_${group}.awk"

awk \
    -v outdir="\${GROUP_OUTPUT_DIR}" \
    -v seqfile="\${GROUP_SEQFILE}" \
    -f "/tmp/split_fasta_${group}.awk" \
    "\${RAW_FASTA}"

echo "[prep_${group}] Genomes separated:"
sed 's/^/  /' "\${GROUP_SEQFILE}"

echo ""
echo "============================================================="
echo "PREPARATION COMPLETED — ${group}"
echo "End              : \$(date '+%Y-%m-%d %H:%M:%S')"
echo "Sequence file    : \${GROUP_SEQFILE}"
echo "Genome FASTAs    : \${GROUP_OUTPUT_DIR}"
echo "============================================================="
EOF

    echo "  → Submitting preparation job for ${group}..."

    # qsub -terse returns the submitted job ID.
    qsub -terse "${job_script}"
}


# -----------------------------------------------------------------------------
# Submit the Cactus pangenome job
#
# This job waits for the corresponding preparation job using -hold_jid.
# -----------------------------------------------------------------------------

submit_cactus_job() {
    local group="$1"
    local preparation_job_id="$2"

    local group_seqfile="${SEQFILES_DIR}/seqfile_${group}.txt"
    local group_output_dir="${OUTPUT_BASE}/${group}"
    local group_work_dir="${WORK_DIR}/${group}"
    local job_script="${CODE_DIR}/job_${group}.sh"

    mkdir -p \
        "${group_output_dir}" \
        "${group_work_dir}"

    cat > "${job_script}" << EOF
#!/bin/bash

#$ -N cactus_${group}
#$ -q all.q
#$ -cwd
#$ -pe smp 10
#$ -V
#$ -o ${LOGS_DIR}/cactus_${group}_\$JOB_ID.out
#$ -e ${LOGS_DIR}/cactus_${group}_\$JOB_ID.err

# =============================================================
# Cactus pangenome job for ${group}
#
# Runs cactus-pangenome using the sequence file generated by the
# corresponding preparation job.
# =============================================================

set -euo pipefail

# =========================
# CONFIGURATION
# =========================

GROUP_SEQFILE="${group_seqfile}"

OUTPUT_DIR="${group_output_dir}"

WORK_DIR="${group_work_dir}"

JOBSTORE="\${WORK_DIR}/jobstore"

LOG_FILE="${LOGS_DIR}/cactus_${group}.log"

OUTPUT_NAME="cactus_${group}_v1"

REFERENCE_GENOME="Mror_B3_${group^}"

THREADS=10

MINIGRAPH_MEMORY="100Gi"

CONTAINER_IMAGE="docker://quay.io/comparative-genomics-toolkit/cactus:latest"

# =========================
# ENVIRONMENT
# =========================

module load singularity-ce/3.11.2

mkdir -p "\${WORK_DIR}"

# Remove the job store from previous attempts.
rm -rf "\${JOBSTORE}"

# =========================
# RUN INFORMATION
# =========================

echo "============================================================="
echo "CACTUS PANGENOME — ${group}"
echo "Start            : \$(date '+%Y-%m-%d %H:%M:%S')"
echo "Job ID           : \${JOB_ID:-not_available}"
echo "Node             : \$(hostname)"
echo "GROUP_SEQFILE    : \${GROUP_SEQFILE}"
echo "OUTPUT_DIR       : \${OUTPUT_DIR}"
echo "WORK_DIR         : \${WORK_DIR}"
echo "JOBSTORE         : \${JOBSTORE}"
echo "OUTPUT_NAME      : \${OUTPUT_NAME}"
echo "REFERENCE_GENOME : \${REFERENCE_GENOME}"
echo "THREADS          : \${THREADS}"
echo "MEMORY           : \${MINIGRAPH_MEMORY}"
echo "CONTAINER        : \${CONTAINER_IMAGE}"
echo "============================================================="
echo ""

# =========================
# RUN CACTUS
# =========================

singularity exec \
    -H "\$(pwd)" \
    -B /Storage:/Storage \
    "\${CONTAINER_IMAGE}" \
    cactus-pangenome \
        "\${JOBSTORE}" \
        "\${GROUP_SEQFILE}" \
        --outDir "\${OUTPUT_DIR}" \
        --logFile "\${LOG_FILE}" \
        --outName "\${OUTPUT_NAME}" \
        --reference "\${REFERENCE_GENOME}" \
        --giraffe \
        --viz \
        --odgi \
        --chrom-vg \
        --chrom-og \
        --gbz \
        --gfa \
        --vcf \
        --workDir "\${WORK_DIR}" \
        --consCores "\${THREADS}" \
        --mgMemory "\${MINIGRAPH_MEMORY}"

# =========================
# FINAL REPORT
# =========================

echo ""
echo "============================================================="
echo "CACTUS PANGENOME COMPLETED — ${group}"
echo "End          : \$(date '+%Y-%m-%d %H:%M:%S')"
echo "Results      : \${OUTPUT_DIR}"
echo "Log file     : \${LOG_FILE}"
echo "============================================================="
EOF

    echo \
        "  → Submitting Cactus job for ${group} " \
        "(waiting for job ${preparation_job_id})..."

    qsub \
        -hold_jid "${preparation_job_id}" \
        "${job_script}"
}

# =========================
# MAIN WORKFLOW
# =========================

echo "============================================================="
echo "CHROMOSOME-LEVEL PANGENOME JOB ORCHESTRATOR"
echo "Start       : $(date '+%Y-%m-%d %H:%M:%S')"
echo "SEQFILE     : ${SEQFILE}"
echo "OUTPUT_BASE : ${OUTPUT_BASE}"
echo "CODE_DIR    : ${CODE_DIR}"
echo "LOGS_DIR    : ${LOGS_DIR}"
echo "============================================================="

while IFS=$'\t ' read -r group url; do
    # Skip empty lines and comments.
    [[ -z "${group}" || "${group}" =~ ^# ]] && continue

    echo ""
    echo "[${group}]"

    # Job 1: preparation
    preparation_job_id=$(
        submit_preparation_job \
            "${group}" \
            "${url}"
    )

    echo "  ✓ Preparation job ID: ${preparation_job_id}"

    # Job 2: Cactus, held until the preparation job finishes
    submit_cactus_job \
        "${group}" \
        "${preparation_job_id}"

done < "${SEQFILE}"

# =========================
# FINAL REPORT
# =========================

echo ""
echo "============================================================="
echo "ALL JOBS SUBMITTED"
echo "End       : $(date '+%Y-%m-%d %H:%M:%S')"
echo "Queue     : Check job status with: qstat"
echo "Logs      : ${LOGS_DIR}"
echo "============================================================="
```
