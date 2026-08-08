#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -pe smp 16
#$ -V

module load singularity-ce/3.11.2

# ==========================================
# Paths
# ==========================================

WORKDIR="/Storage/data1/isabella.gallego/MAESTRIA/PGGB_results_pacbio"

PARTITIONS="${WORKDIR}/partition_results"

SIF="/Storage/data1/isabella.gallego/MAESTRIA/pggb_latest.sif"

OUTBASE="${WORKDIR}/pggb_partitioned_results"

mkdir -p "$OUTBASE"

cd "$WORKDIR"

echo "========================================="
echo "Starting PGGB for all communities"
echo "========================================="

# ==========================================
# Loop over all community FASTAs
# ==========================================

for fasta in ${PARTITIONS}/*.fa; do

    # Obtain community name
    community=$(basename "$fasta" .fa)

    echo ""
    echo "-----------------------------------------"
    echo "Processing: $community"
    echo "Input: $fasta"
    echo "-----------------------------------------"

    # Output directory for each community
    OUTDIR="${OUTBASE}/${community}"

    mkdir -p "$OUTDIR"

    # ==========================================
    # Run PGGB
    # ==========================================

    singularity exec \
        --bind /Storage \
        "$SIF" pggb \
            -i "$fasta" \
            -o "$OUTDIR" \
            -t 16 \
            -p 95 \
            -n 5 \
            -s 2000 \
            -k 29 \
            -G 5000
            -V 'B3:1000'

    echo ""
    echo "Finished: $community"

done

echo ""
echo "========================================="
echo "All PGGB runs completed"
echo "========================================="
