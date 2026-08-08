```bash
#!/bin/bash

#$ -N group_chromosomes              # Job name
#$ -q all.q                          # Queue
#$ -cwd                              # Run from the current working directory
#$ -V                                # Export environment variables

# =============================================================
# Chromosome grouping by syntenic group
#
# Searches for chromosome FASTA files from all strains and groups
# sequences belonging to the same chromosome group into a single
# combined FASTA file.
#
# Expected input filename format:
#   mror_<strain>_group<number>.fasta
#
# Example:
#   mror_B3_group1.fasta
#   mror_C26_group1.fasta
#
# Output:
#   One combined FASTA file per chromosome group.
#
# Example:
#   grouped_by_chromosome/group1.fasta
#   grouped_by_chromosome/group2.fasta
#
# FASTA headers are modified by adding the strain name as a prefix.
#
# Usage:
#   qsub group_chromosomes.sh
# =============================================================

set -euo pipefail

# =========================
# CONFIGURATION
# =========================

# Base directory containing the FASTA files for each strain
STRAIN_DIR="/Storage/data1/isabella.gallego/MAESTRIA/data/pacbio"

# Directory where grouped FASTA files will be created
OUTPUT_DIR="grouped_by_chromosome"

# =========================
# SETUP
# =========================

mkdir -p "$OUTPUT_DIR"

echo "============================================================="
echo "CHROMOSOME GROUPING"
echo "Start        : $(date '+%Y-%m-%d %H:%M:%S')"
echo "Job ID       : ${JOB_ID:-not_available}"
echo "Node         : $(hostname)"
echo "STRAIN_DIR   : $STRAIN_DIR"
echo "OUTPUT_DIR   : $OUTPUT_DIR"
echo "============================================================="
echo ""

# =========================
# GROUP CHROMOSOME FASTA FILES
# =========================

find "$STRAIN_DIR" \
    -type f \
    -name "mror_*_group*.fasta" \
| while read -r fasta_file; do

    # Extract the chromosome group name, for example "group1"
    group_name=$(
        basename "$fasta_file" \
        | grep -o "group[0-9]\+"
    )

    # Extract the strain name from the filename
    strain_name=$(
        basename "$fasta_file" \
        | cut -d'_' -f2
    )

    output_file="${OUTPUT_DIR}/${group_name}.fasta"

    echo "Processing: $fasta_file"
    echo "  Strain : $strain_name"
    echo "  Group  : $group_name"
    echo "  Output : $output_file"

    # Append the sequence to the corresponding chromosome-group FASTA.
    # The strain name is added as a prefix to each FASTA header.
    awk \
        -v strain="$strain_name" \
        '
        /^>/ {
            print ">" strain "_" substr($0, 2)
        }

        !/^>/ {
            print
        }
        ' \
        "$fasta_file" \
        >> "$output_file"

done

# =========================
# FINAL REPORT
# =========================

echo ""
echo "============================================================="
echo "CHROMOSOME GROUPING COMPLETED"
echo "End          : $(date '+%Y-%m-%d %H:%M:%S')"
echo "Results saved to:"
echo "  $OUTPUT_DIR"
echo "============================================================="
```
