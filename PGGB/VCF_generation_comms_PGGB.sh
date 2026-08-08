#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -pe smp 8
#$ -V

module load singularity-ce/3.11.2
module load Samtools/1.22

# ==========================================
# Paths
# ==========================================

BASE_DIR="/Storage/data1/isabella.gallego/MAESTRIA/PGGB_results_pacbio/pggb_partitioned_results"

SIF="/Storage/data1/isabella.gallego/MAESTRIA/pggb_latest.sif"

THREADS=8

# ==========================================
# Reference strain
# ==========================================

REF="B3"

# ==========================================
# Loop through communities
# ==========================================

cd "$BASE_DIR"

for COMMUNITY_DIR in all_pacbio_pansn*.community.*; do

    echo ""
    echo "======================================="
    echo "Processing $COMMUNITY_DIR"
    echo "======================================="

    cd "$COMMUNITY_DIR"

    # ==========================================
    # Find GFA
    # ==========================================

    GFA=$(find . -name "*.final.gfa" | head -n 1)

    if [[ ! -f "$GFA" ]]; then
        echo "[WARNING] No final GFA found"
        cd "$BASE_DIR"
        continue
    fi

    echo "GFA found:"
    echo "$GFA"

    # ==========================================
    # Convert GFA -> VG
    # ==========================================

    echo ""
    echo "Converting GFA to VG..."

    singularity exec \
        --bind /Storage \
        "$SIF" vg convert \
            -g "$GFA" > graph.vg

    # ==========================================
    # Index graph
    # ==========================================

    echo ""
    echo "Indexing graph..."

    singularity exec \
        --bind /Storage \
        "$SIF" vg index \
            -x graph.xg \
            graph.vg

    # ==========================================
    # Generate VCF
    # ==========================================

    echo ""
    echo "Generating VCF..."

    singularity exec \
        --bind /Storage \
        "$SIF" vg deconstruct \
            -P "$REF" \
            -t "$THREADS" \
            -e \
            graph.xg > variants.vcf

    # ==========================================
    # Compress + index
    # ==========================================

    bgzip -f variants.vcf
    tabix -f -p vcf variants.vcf.gz

    echo ""
    echo "VCF generated:"
    echo "${COMMUNITY_DIR}/variants.vcf.gz"

    cd "$BASE_DIR"

done

echo ""
echo "======================================="
echo "ALL COMMUNITIES FINISHED"
echo "======================================="
