#!/bin/bash

#$ -N sequences_per_community
#$ -q all.q
#$ -cwd
#$ -pe smp 2
#$ -V
#$ -o /Storage/data1/isabella.gallego/MAESTRIA/code/PGGB_24_genomes/logs/sequences_per_comm_$JOB_ID.out
#$ -e /Storage/data1/isabella.gallego/MAESTRIA/code/PGGB_24_genomes/logs/sequences_per_comm_$JOB_ID.err


# ==========================================
# Paths
# ==========================================

PARTITION_DIR="/Storage/data1/isabella.gallego/MAESTRIA/PGGB_results_24G/24G_partition_results"

OUTFILE="/Storage/data1/isabella.gallego/MAESTRIA/PGGB_results_24G/community_sequence_summary_24G.tsv"

cd "$PARTITION_DIR"

# ==========================================
# Header
# ==========================================

echo -e "community\tsequence_name" > "$OUTFILE"

echo "======================================="
echo "Listing sequences in each community"
echo "======================================="

# ==========================================
# Loop over community fasta files
# ==========================================
for FA in *.community.*.fa; do

    echo "Processing $FA ..."

    COMMUNITY=$(basename "$FA" | sed 's/\.fa//')

    grep "^>" "$FA" | sed 's/^>//' | while read HEADER; do

        echo -e "${COMMUNITY}\t${HEADER}" >> "$OUTFILE"

    done

done

echo ""
echo "======================================="
echo "DONE"
echo "======================================="

echo ""
echo "Results saved in:"
echo "$OUTFILE"
