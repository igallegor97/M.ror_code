#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -pe smp 16
#$ -V

module load singularity-ce/3.11.2
module load Samtools/1.22

# Paths
WORKDIR="/Storage/data1/isabella.gallego/MAESTRIA/PGGB_results_pacbio"
SIF="/Storage/data1/isabella.gallego/MAESTRIA/pggb_latest.sif"

cd "$WORKDIR"

INPUT_FASTA="/Storage/data1/isabella.gallego/MAESTRIA/PGGB_results_pacbio/all_pacbio.fasta"
OUTDIR="${WORKDIR}/partition_results"

mkdir -p "$OUTDIR"

echo "Running partition-before-pggb..."

singularity exec \
    --bind /Storage \
    "$SIF" partition-before-pggb \
        -i "$INPUT_FASTA" \
        -o "$OUTDIR" \
        -n 5 \
        -t 16 \
        -p 90 \
        -s 5k \
        -V 'B3:1000'

echo "Partition finished."