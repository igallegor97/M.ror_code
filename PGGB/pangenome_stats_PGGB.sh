#!/bin/bash
#$ -N pggb_stats
#$ -q all.q
#$ -cwd
#$ -pe smp 4
#$ -V
#$ -o /Storage/data1/isabella.gallego/MAESTRIA/code/minigraph_cactus/chroms_cactus/logs/pggb_stats_$JOB_ID.out
#$ -e /Storage/data1/isabella.gallego/MAESTRIA/code/minigraph_cactus/chroms_cactus/logs/pggb_stats_$JOB_ID.err

# =============================================================================
# pangenome_stats_pggb.sh — Job SGE
# Extrae estadísticas de pangenomas generados con PGGB
# =============================================================================

module load bcftools/1.22
module load bedtools/2.28.0

set -euo pipefail

PGGB_BASE="/Storage/data1/isabella.gallego/MAESTRIA/PGGB_results_pacbio/pggb_partitioned_results"
GFF="/Storage/data1/isabella.gallego/MAESTRIA/data/genomes/MrorB3.groups.gff3"
STATSDIR="/Storage/data1/isabella.gallego/MAESTRIA/PGGB_results_pacbio/pangenome_stats_pggb"
PREFIX="all_pacbio_pansn.fasta.bf3285f"

mkdir -p "${STATSDIR}/vcf" "${STATSDIR}/gfa" "${STATSDIR}/functional"

# Comunidades 0-9
P_COMMUNITIES=(0 1 2 3 4 5 6 7 8 9)

echo "========================================"
echo " Extrayendo estadísticas PGGB"
echo " Output: ${STATSDIR}"
echo " bcftools: $(which bcftools)"
echo " bedtools: $(which bedtools)"
echo "========================================"

# ── 1. ESTADÍSTICAS DEL VCF ───────────────────────────────────────────────────
echo ""
echo "[1/4] Estadísticas de variantes (VCF)..."

VARIANTS_SUMMARY="${STATSDIR}/vcf/variants_summary.tsv"
printf "community\tSNPs\tINDELs\tSVs\tDEL\tINS\tINV\tDUP\tBND\ttotal\n" > "${VARIANTS_SUMMARY}"

for com in "${P_COMMUNITIES[@]}"; do
    comdir="${PGGB_BASE}/${PREFIX}.community.${com}"
    vcf="${comdir}/variants.vcf"

    if [[ ! -f "${vcf}" ]]; then
        echo "  ⚠ No encontrado: ${vcf}"
        continue
    fi

    echo "  Procesando community.${com}..."

    # bcftools stats necesita VCF indexado o sin comprimir — usar directamente
    bcftools stats "${vcf}" > "${STATSDIR}/vcf/community${com}_bcftools_stats.txt"

    snps=$(grep   "^SN.*number of SNPs:"   "${STATSDIR}/vcf/community${com}_bcftools_stats.txt" | cut -f4)
    indels=$(grep "^SN.*number of indels:" "${STATSDIR}/vcf/community${com}_bcftools_stats.txt" | cut -f4)

    # PGGB tampoco usa SVTYPE — clasificar por tamaño REF/ALT
    del=$(bcftools query -f '%REF\t%ALT\n' "${vcf}" | \
        awk 'length($1)>=50 && length($1)>length($2)' | wc -l) || del=0
    ins=$(bcftools query -f '%REF\t%ALT\n' "${vcf}" | \
        awk 'length($2)>=50 && length($2)>length($1)' | wc -l) || ins=0
    inv=0; dup=0; bnd=0

    snps=${snps:-0}; indels=${indels:-0}
    del=${del:-0}; ins=${ins:-0}
    svs=$((del + ins))
    total=$((snps + indels + svs))

    printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n" \
        "community.${com}" "${snps}" "${indels}" "${svs}" \
        "${del}" "${ins}" "${inv}" "${dup}" "${bnd}" "${total}" \
        >> "${VARIANTS_SUMMARY}"

    # Raw variants
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO\n' "${vcf}" \
        > "${STATSDIR}/vcf/community${com}_variants_raw.tsv"

    echo "    SNPs=${snps} INDELs=${indels} DEL=${del} INS=${ins} total=${total}"
done
echo "  ✓ Resumen de variantes: ${VARIANTS_SUMMARY}"

# ── 2. ESTADÍSTICAS DEL GRAFO (GFA / odgi) ───────────────────────────────────
echo ""
echo "[2/4] Estadísticas del grafo..."

GFA_SUMMARY="${STATSDIR}/gfa/graph_stats.tsv"
printf "community\tnodes\tedges\tpaths\tlength_bp\tgfa_size_MB\n" > "${GFA_SUMMARY}"

for com in "${P_COMMUNITIES[@]}"; do
    comdir="${PGGB_BASE}/${PREFIX}.community.${com}"
    gfa=$(ls "${comdir}"/*.smooth.final.gfa 2>/dev/null | head -1)
    odgi_stats="${comdir}/odgi_stats.txt"

    if [[ ! -f "${gfa}" ]]; then
        echo "  ⚠ No encontrado GFA en: ${comdir}"
        continue
    fi

    echo "  Procesando community.${com}..."

    size_mb=$(du -m "${gfa}" | cut -f1)

    # Parsear GFA directamente
    nodes=$(awk '$1=="S"' "${gfa}" | wc -l)
    edges=$(awk '$1=="L"' "${gfa}" | wc -l)
    paths=$(awk '$1=="P" || $1=="W"' "${gfa}" | wc -l)

    # Longitud total del grafo desde odgi_stats.txt si existe
    length_bp="NA"
    if [[ -f "${odgi_stats}" ]]; then
        length_bp=$(grep -i "length" "${odgi_stats}" | head -1 | awk '{print $NF}' || echo "NA")
    fi

    printf "%s\t%d\t%d\t%d\t%s\t%d\n" \
        "community.${com}" "${nodes}" "${edges}" "${paths}" "${length_bp}" "${size_mb}" \
        >> "${GFA_SUMMARY}"

    echo "    nodes=${nodes} edges=${edges} paths=${paths} length=${length_bp} size=${size_mb}MB"
done
echo "  ✓ Estadísticas del grafo: ${GFA_SUMMARY}"

# ── 3. ANOTACIÓN FUNCIONAL (VCF x GFF) ───────────────────────────────────────
echo ""
echo "[3/4] Intersección de variantes con anotación funcional..."

FUNC_SUMMARY="${STATSDIR}/functional/variants_by_region.tsv"
printf "community\tfeature_type\tvariant_count\n" > "${FUNC_SUMMARY}"

FEATURE_BED="${STATSDIR}/functional/features.bed"

# Convertir GFF a BED
# CHROM en VCF PGGB: B3#0#Group1
# CHROM en GFF:      MrorB3_Group1
# → normalizar GFF a B3#0#Group1
awk '
BEGIN{OFS="\t"}
!/^#/ && NF>=8 {
    chrom=$1
    # MrorB3_Group1 -> B3#0#Group1
    gsub(/^MrorB3_/, "", chrom)
    chrom = "B3#0#" chrom
    s=$4-1; e=$5+0
    if(s>e){tmp=s; s=e; e=tmp}
    if(s<0) s=0
    print chrom, s, e, $3
}
' "${GFF}" | sort -k1,1 -k2,2n > "${FEATURE_BED}"

echo "  Features en BED: $(wc -l < ${FEATURE_BED})"
echo "  Ejemplo de CHROM en BED: $(head -1 ${FEATURE_BED} | cut -f1)"

for com in "${P_COMMUNITIES[@]}"; do
    comdir="${PGGB_BASE}/${PREFIX}.community.${com}"
    vcf="${comdir}/variants.vcf"
    [[ ! -f "${vcf}" ]] && continue

    echo "  Procesando community.${com}..."

    vcf_bed="${STATSDIR}/functional/community${com}_variants.bed"

    # Convertir VCF a BED
    bcftools query -f '%CHROM\t%POS0\t%POS\t%TYPE\n' "${vcf}" > "${vcf_bed}"

    echo "    Variantes en BED: $(wc -l < ${vcf_bed})"
    echo "    Ejemplo CHROM en VCF BED: $(head -1 ${vcf_bed} | cut -f1)"

    for feature in gene mRNA CDS exon intron three_prime_UTR five_prime_UTR; do
        feature_count=$(
            awk -v f="${feature}" '$4==f' "${FEATURE_BED}" | \
            bedtools intersect -a "${vcf_bed}" -b stdin -wa 2>/dev/null | \
            wc -l
        ) || feature_count=0
        printf "%s\t%s\t%d\n" "community.${com}" "${feature}" "${feature_count}" >> "${FUNC_SUMMARY}"
    done

    intergenic=$(
        bedtools intersect -a "${vcf_bed}" -b "${FEATURE_BED}" -v 2>/dev/null | wc -l
    ) || intergenic=0
    printf "%s\tintergenic\t%d\n" "community.${com}" "${intergenic}" >> "${FUNC_SUMMARY}"

    echo "    ✓ Anotación funcional completada"
done
echo "  ✓ Distribución funcional: ${FUNC_SUMMARY}"

# ── 4. PRESENCIA/AUSENCIA POR MUESTRA (PAV) ───────────────────────────────────
echo ""
echo "[4/4] Presencia/ausencia de variantes por muestra..."

PAV_SUMMARY="${STATSDIR}/vcf/pav_summary.tsv"
printf "community\tsample\tprivate_variants\tshared_variants\n" > "${PAV_SUMMARY}"

for com in "${P_COMMUNITIES[@]}"; do
    comdir="${PGGB_BASE}/${PREFIX}.community.${com}"
    vcf="${comdir}/variants.vcf"
    [[ ! -f "${vcf}" ]] && continue

    echo "  Procesando community.${com}..."
    samples=$(bcftools query -l "${vcf}" 2>/dev/null)

    for sample in ${samples}; do
        shared=$(bcftools view -s "${sample}" "${vcf}" 2>/dev/null | \
            bcftools view -e 'GT="0/0" || GT="./."' 2>/dev/null | \
            grep -c "^[^#]") || shared=0

        private=$(bcftools view -s "${sample}" "${vcf}" 2>/dev/null | \
            bcftools view -e 'GT="0/0" || GT="./."' 2>/dev/null | \
            bcftools view -i 'AC=1' 2>/dev/null | \
            grep -c "^[^#]") || private=0

        printf "%s\t%s\t%s\t%s\n" \
            "community.${com}" "${sample}" "${private}" "${shared}" >> "${PAV_SUMMARY}"
        echo "    ${sample}: shared=${shared} private=${private}"
    done
done
echo "  ✓ PAV summary: ${PAV_SUMMARY}"

echo ""
echo "========================================"
echo " ✓ Todas las estadísticas generadas en:"
echo "   ${STATSDIR}/"
echo "   vcf/variants_summary.tsv          → conteos por tipo"
echo "   vcf/pav_summary.tsv               → presencia/ausencia por muestra"
echo "   gfa/graph_stats.tsv               → nodes, edges, paths, length, size"
echo "   functional/variants_by_region.tsv → distribución por región génica"
echo "========================================"
