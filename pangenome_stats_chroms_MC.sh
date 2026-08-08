```bash
#!/bin/bash

#$ -N pangenome_stats
#$ -q all.q
#$ -cwd
#$ -pe smp 4
#$ -V
#$ -o /Storage/data1/isabella.gallego/MAESTRIA/code/minigraph_cactus/chroms_cactus/logs/pangenome_stats_$JOB_ID.out
#$ -e /Storage/data1/isabella.gallego/MAESTRIA/code/minigraph_cactus/chroms_cactus/logs/pangenome_stats_$JOB_ID.err

# =============================================================
# pangenome_stats.sh
#
# SGE job that extracts summary statistics from chromosome-level
# pangenomes generated with minigraph-cactus.
#
# Analyses performed:
#
#   1. Variant statistics from VCF files.
#   2. Graph statistics from GFA files.
#   3. Exclusive functional classification by genomic region.
#   4. Variant presence/absence by sample.
#
# Note:
#   The loop variable is named "grp" to avoid conflicts with the
#   SGE-exported environment variable "group", which may contain
#   the system group identifier.
#
# Usage:
#   qsub pangenome_stats.sh
# =============================================================

set -euo pipefail

# =========================
# ENVIRONMENT
# =========================

module load bcftools/1.22
module load bedtools/2.28.0

# =========================
# CONFIGURATION
# =========================

OUTPUT_BASE="/Storage/data1/isabella.gallego/MAESTRIA/cactus_chroms"

GFF_FILE="/Storage/data1/isabella.gallego/MAESTRIA/data/genomes/MrorB3.groups.gff3"

STATS_DIR="/Storage/data1/isabella.gallego/MAESTRIA/cactus_chroms/pangenome_stats_cactus_chroms"

PANGENOME_GROUPS=(
    group1
    group2
    group3
    group4
    group5
    group6
    group7
    group8
    group9
    group10
    group11
)

# =========================
# SETUP
# =========================

mkdir -p \
    "${STATS_DIR}/vcf" \
    "${STATS_DIR}/gfa" \
    "${STATS_DIR}/functional"

# =========================
# RUN INFORMATION
# =========================

echo "============================================================="
echo "PANGENOME STATISTICS"
echo "Start          : $(date '+%Y-%m-%d %H:%M:%S')"
echo "Job ID         : ${JOB_ID:-not_available}"
echo "Node           : $(hostname)"
echo "OUTPUT_BASE    : ${OUTPUT_BASE}"
echo "GFF_FILE       : ${GFF_FILE}"
echo "STATS_DIR      : ${STATS_DIR}"
echo "bcftools       : $(which bcftools)"
echo "bedtools       : $(which bedtools)"
echo "Group1 VCF     : $(ls "${OUTPUT_BASE}/group1/cactus_group1_v1.vcf.gz" 2>/dev/null && echo YES || echo NO)"
echo "Group1 GFA     : $(ls "${OUTPUT_BASE}/group1/cactus_group1_v1.gfa.gz" 2>/dev/null && echo YES || echo NO)"
echo "GFF available  : $(ls "${GFF_FILE}" 2>/dev/null && echo YES || echo NO)"
echo "============================================================="

# =========================
# 1. VCF STATISTICS
# =========================

echo ""
echo "[1/4] Extracting variant statistics from VCF files..."

VARIANTS_SUMMARY="${STATS_DIR}/vcf/variants_summary.tsv"

printf "group\tSNPs\tINDELs\tSVs\tDEL\tINS\tINV\tDUP\tBND\ttotal\n" \
    > "${VARIANTS_SUMMARY}"

for grp in "${PANGENOME_GROUPS[@]}"; do
    vcf="${OUTPUT_BASE}/${grp}/cactus_${grp}_v1.vcf.gz"

    if [[ ! -f "${vcf}" ]]; then
        echo "  [WARNING] File not found: ${vcf}"
        continue
    fi

    echo "  Processing ${grp}..."

    bcftools stats "${vcf}" \
        > "${STATS_DIR}/vcf/${grp}_bcftools_stats.txt"

    snps=$(
        grep "^SN.*number of SNPs:" \
            "${STATS_DIR}/vcf/${grp}_bcftools_stats.txt" \
        | cut -f4
    )

    indels=$(
        grep "^SN.*number of indels:" \
            "${STATS_DIR}/vcf/${grp}_bcftools_stats.txt" \
        | cut -f4
    )

    # Structural deletion:
    # REF length >= 50 bp and REF is longer than ALT.
    del=$(
        bcftools query \
            -f '%REF\t%ALT\n' \
            "${vcf}" \
        | awk '
            length($1) >= 50 &&
            length($1) > length($2)
        ' \
        | wc -l
    ) || del=0

    # Structural insertion:
    # ALT length >= 50 bp and ALT is longer than REF.
    ins=$(
        bcftools query \
            -f '%REF\t%ALT\n' \
            "${vcf}" \
        | awk '
            length($2) >= 50 &&
            length($2) > length($1)
        ' \
        | wc -l
    ) || ins=0

    inv=0
    dup=0
    bnd=0

    snps=${snps:-0}
    indels=${indels:-0}
    del=${del:-0}
    ins=${ins:-0}
    inv=${inv:-0}
    dup=${dup:-0}
    bnd=${bnd:-0}

    svs=$((del + ins + inv + dup + bnd))
    total=$((snps + indels + svs))

    printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n" \
        "${grp}" \
        "${snps}" \
        "${indels}" \
        "${svs}" \
        "${del}" \
        "${ins}" \
        "${inv}" \
        "${dup}" \
        "${bnd}" \
        "${total}" \
        >> "${VARIANTS_SUMMARY}"

    bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO\n' \
        "${vcf}" \
        > "${STATS_DIR}/vcf/${grp}_variants_raw.tsv"

    echo \
        "    SNPs=${snps} INDELs=${indels} DEL=${del} " \
        "INS=${ins} INV=${inv} DUP=${dup} BND=${bnd}"
done

echo "  [OK] Variant summary: ${VARIANTS_SUMMARY}"

# =========================
# 2. GRAPH STATISTICS
# =========================

echo ""
echo "[2/4] Extracting graph statistics from GFA files..."

GFA_SUMMARY="${STATS_DIR}/gfa/graph_stats.tsv"

printf "group\tnodes\tedges\tpaths\tsnarls\tgfa_size_MB\n" \
    > "${GFA_SUMMARY}"

for grp in "${PANGENOME_GROUPS[@]}"; do
    gfa="${OUTPUT_BASE}/${grp}/cactus_${grp}_v1.gfa.gz"

    if [[ ! -f "${gfa}" ]]; then
        echo "  [WARNING] File not found: ${gfa}"
        continue
    fi

    echo "  Processing ${grp}..."

    size_mb=$(
        du -m "${gfa}" \
        | cut -f1
    )

    nodes=$(
        zcat "${gfa}" \
        | awk '$1 == "S"' \
        | wc -l
    )

    edges=$(
        zcat "${gfa}" \
        | awk '$1 == "L"' \
        | wc -l
    )

    paths=$(
        zcat "${gfa}" \
        | awk '$1 == "P" || $1 == "W"' \
        | wc -l
    )

    snarls_file="${OUTPUT_BASE}/${grp}/cactus_${grp}_v1.snarls"

    if [[ -f "${snarls_file}" ]]; then
        snarls=$(
            bcftools query \
                -f '%CHROM\t%POS\n' \
                "${OUTPUT_BASE}/${grp}/cactus_${grp}_v1.vcf.gz" \
                2>/dev/null \
            | sort -u \
            | wc -l \
            || echo "NA"
        )
    else
        snarls="NA"
    fi

    printf "%s\t%d\t%d\t%d\t%s\t%d\n" \
        "${grp}" \
        "${nodes}" \
        "${edges}" \
        "${paths}" \
        "${snarls}" \
        "${size_mb}" \
        >> "${GFA_SUMMARY}"

    echo \
        "    nodes=${nodes} edges=${edges} paths=${paths} " \
        "snarls=${snarls} size=${size_mb}MB"
done

echo "  [OK] Graph statistics: ${GFA_SUMMARY}"

# =========================
# 3. EXCLUSIVE FUNCTIONAL ANNOTATION
# =========================

echo ""
echo "[3/4] Performing exclusive functional classification..."

GFF_GROUPS="/Storage/data1/isabella.gallego/MAESTRIA/data/genomes/MrorB3.groups.gff3"

FUNCTIONAL_SUMMARY="${STATS_DIR}/functional/variants_by_region_exclusive.tsv"

printf "group\tCDS\tUTR\tintron\tintergenic\ttotal\n" \
    > "${FUNCTIONAL_SUMMARY}"

FEATURE_BED="${STATS_DIR}/functional/features.bed"

awk '
BEGIN {
    OFS = "\t"
}

!/^#/ && NF >= 8 {
    chrom = $1

    gsub(
        /^MrorB3_/,
        "Mror_B3_",
        chrom
    )

    print chrom, $4 - 1, $5, $3
}
' "${GFF_GROUPS}" \
| sort -k1,1 -k2,2n \
> "${FEATURE_BED}"

echo "  Features in BED: $(wc -l < "${FEATURE_BED}")"

for grp in "${PANGENOME_GROUPS[@]}"; do
    vcf="${OUTPUT_BASE}/${grp}/cactus_${grp}_v1.vcf.gz"

    [[ ! -f "${vcf}" ]] && continue

    echo "  Processing ${grp}..."

    variants_bed="${STATS_DIR}/functional/${grp}_variants.bed"

    bcftools query \
        -f '%CHROM\t%POS0\t%POS\t%TYPE\n' \
        "${vcf}" \
        > "${variants_bed}"

    total=$(
        wc -l < "${variants_bed}"
    )

    # ---------------------------------------------------------
    # CDS
    # ---------------------------------------------------------

    awk '$4 == "CDS"' \
        "${FEATURE_BED}" \
    | bedtools intersect \
        -a "${variants_bed}" \
        -b stdin \
        -u \
        > "${STATS_DIR}/functional/${grp}_CDS.bed"

    cds=$(
        wc -l < "${STATS_DIR}/functional/${grp}_CDS.bed"
    )

    # ---------------------------------------------------------
    # UTR, excluding CDS
    # ---------------------------------------------------------

    awk '
        $4 == "three_prime_UTR" ||
        $4 == "five_prime_UTR"
    ' "${FEATURE_BED}" \
    | bedtools intersect \
        -a "${variants_bed}" \
        -b stdin \
        -u \
    | bedtools intersect \
        -a stdin \
        -b "${STATS_DIR}/functional/${grp}_CDS.bed" \
        -v \
        > "${STATS_DIR}/functional/${grp}_UTR.bed"

    utr=$(
        wc -l < "${STATS_DIR}/functional/${grp}_UTR.bed"
    )

    # ---------------------------------------------------------
    # INTRON, excluding CDS and UTR
    # ---------------------------------------------------------

    awk '$4 == "intron"' \
        "${FEATURE_BED}" \
    | bedtools intersect \
        -a "${variants_bed}" \
        -b stdin \
        -u \
    | bedtools intersect \
        -a stdin \
        -b "${STATS_DIR}/functional/${grp}_CDS.bed" \
        -v \
    | bedtools intersect \
        -a stdin \
        -b "${STATS_DIR}/functional/${grp}_UTR.bed" \
        -v \
        > "${STATS_DIR}/functional/${grp}_intron.bed"

    intron=$(
        wc -l < "${STATS_DIR}/functional/${grp}_intron.bed"
    )

    # ---------------------------------------------------------
    # INTERGENIC
    # ---------------------------------------------------------

    cat \
        "${STATS_DIR}/functional/${grp}_CDS.bed" \
        "${STATS_DIR}/functional/${grp}_UTR.bed" \
        "${STATS_DIR}/functional/${grp}_intron.bed" \
        > "${STATS_DIR}/functional/${grp}_genic.bed"

    intergenic=$(
        bedtools intersect \
            -a "${variants_bed}" \
            -b "${STATS_DIR}/functional/${grp}_genic.bed" \
            -v \
        | wc -l
    )

    printf "%s\t%d\t%d\t%d\t%d\t%d\n" \
        "${grp}" \
        "${cds}" \
        "${utr}" \
        "${intron}" \
        "${intergenic}" \
        "${total}" \
        >> "${FUNCTIONAL_SUMMARY}"

    echo \
        "    CDS=${cds} UTR=${utr} intron=${intron} " \
        "intergenic=${intergenic}"
done

echo "  [OK] Exclusive functional table: ${FUNCTIONAL_SUMMARY}"

# =========================
# 4. SAMPLE-LEVEL PRESENCE/ABSENCE
# =========================

echo ""
echo "[4/4] Calculating variant presence/absence by sample..."

PAV_SUMMARY="${STATS_DIR}/vcf/pav_summary.tsv"

printf "group\tsample\tprivate_variants\tshared_variants\n" \
    > "${PAV_SUMMARY}"

for grp in "${PANGENOME_GROUPS[@]}"; do
    vcf="${OUTPUT_BASE}/${grp}/cactus_${grp}_v1.vcf.gz"

    [[ ! -f "${vcf}" ]] && continue

    echo "  Processing ${grp}..."

    samples=$(
        bcftools query \
            -l "${vcf}" \
            2>/dev/null
    )

    for sample in ${samples}; do
        shared=$(
            bcftools view \
                -s "${sample}" \
                "${vcf}" \
                2>/dev/null \
            | bcftools view \
                -e 'GT="0/0" || GT="./."' \
                2>/dev/null \
            | grep -c "^[^#]"
        ) || shared=0

        private=$(
            bcftools view \
                -s "${sample}" \
                "${vcf}" \
                2>/dev/null \
            | bcftools view \
                -e 'GT="0/0" || GT="./."' \
                2>/dev/null \
            | bcftools view \
                -i 'AC=1' \
                2>/dev/null \
            | grep -c "^[^#]"
        ) || private=0

        printf "%s\t%s\t%s\t%s\n" \
            "${grp}" \
            "${sample}" \
            "${private}" \
            "${shared}" \
            >> "${PAV_SUMMARY}"

        echo \
            "    ${sample}: shared=${shared} private=${private}"
    done
done

echo "  [OK] Presence/absence summary: ${PAV_SUMMARY}"

# =========================
# FINAL REPORT
# =========================

echo ""
echo "============================================================="
echo "PANGENOME STATISTICS COMPLETED"
echo "End          : $(date '+%Y-%m-%d %H:%M:%S')"
echo "Results      : ${STATS_DIR}/"
echo ""
echo "Generated tables:"
echo "  vcf/variants_summary.tsv"
echo "      Variant counts by type."
echo ""
echo "  vcf/pav_summary.tsv"
echo "      Variant presence/absence by sample."
echo ""
echo "  gfa/graph_stats.tsv"
echo "      Graph nodes, edges, paths, snarls, and file size."
echo ""
echo "  functional/variants_by_region_exclusive.tsv"
echo "      Exclusive variant distribution by genomic region."
echo "============================================================="
```
