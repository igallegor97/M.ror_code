#!/bin/bash

#$ -N verify_VCFs
#$ -q all.q
#$ -cwd
#$ -pe smp 8
#$ -V
#$ -o /Storage/data1/isabella.gallego/MAESTRIA/code/PGGB_24_genomes/logs/verification_vcfs_$JOB_ID.out
#$ -e /Storage/data1/isabella.gallego/MAESTRIA/code/PGGB_24_genomes/logs/verification_vcfs_$JOB_ID.err

BASE="/Storage/data1/isabella.gallego/MAESTRIA/PGGB_results_24G/pggb_partitioned_results_24G"

SUMMARY="${BASE}/vcf_reference_chromosome_communities.tsv"

{
    printf "community\tn_vcf_samples\treference\ttotal_samples_after_reference\tn_records\tn_snps\tindex_ok\tvcf_ok\n"

    for vcf in \
        "${BASE}"/all_genomes_pansn.fasta.bf3285f.community.*/variants.vcf.gz
    do
        [[ -f "${vcf}" ]] || continue

        community="$(
            basename "$(dirname "${vcf}")" \
            | grep -o 'community\.[0-9]\+'
        )"

        reference="$(
            awk -F '\t' \
                -v community="${community}" \
                'NR > 1 && $1 == community {print $2}' \
                "${SUMMARY}"
        )"

        n_samples="$(
            bcftools query -l "${vcf}" \
            | wc -l
        )"

        n_records="$(
            bcftools view -H "${vcf}" \
            | wc -l
        )"

        n_snps="$(
            bcftools view -H -v snps "${vcf}" \
            | wc -l
        )"

        if bcftools view -h "${vcf}" \
            | grep -q '^#CHROM'; then
            vcf_ok="YES"
        else
            vcf_ok="NO"
        fi

        if [[ -f "${vcf}.tbi" || -f "${vcf}.csi" ]]; then
            index_ok="YES"
        else
            index_ok="NO"
        fi

        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "${community}" \
            "${n_samples}" \
            "${reference}" \
            "$((n_samples + 1))" \
            "${n_records}" \
            "${n_snps}" \
            "${index_ok}" \
            "${vcf_ok}"

    done
} | sort -V > "${BASE}/chromosome_vcf_qc.tsv"
