#!/usr/bin/env bash
module load GATK
module load samtools

BFILE=${1}
FILESTEM=$(basename ${BFILE})

# VARS
REFDIR="/dat/${USER}/.hg38"

for file in ${REFDIR}/hg38{_chr.map,_chr.reorder.map,.dict,.fa,.fa.fai}; do
    MISSING=0
    if [[ -f ${file} ]]; then
        echo ${file} missing! Did you run prepare-refdata.sh ?
        let MISSING=${MISSING}+1
    fi
done
if [[ ${MISSING} > 0 ]]; then
    echo "Exiting due to missing ref data"; exit 1 
fi


ln -s ${REFDIR}/hg38* .


# Get GENOTYPE fingerprints from plink bfile.bed
plink --bfile ${BFILE} \
    --recode vcf-iid \
    --output-chr chrM \
    --keep-allele-order \
    --out ${FILESTEM}.fingerprints

echo "Done generating genotype fingerprints for ${FILESTEM}"
exit 0

# Get RNA fingerprints from 10X bam
module load GATK
gatk ExtractFingerprint  \
        --INPUT "${bam}" \
        --HAPLOTYPE_MAP ${HAPMAP} \
        --OUTPUT "${sample}.10Xrna.vcf" \
        --REFERENCE_SEQUENCE 'hg38.fa'

# Fix sample name

function fixVCFsampleName() {
    local in_vcf=${1}
    local sample=${2}
    local before=$(awk -vOFS='\t' '$1 ~ /^#CHROM/ {print $9,$10}' ${in_vcf})
    local after="FORMAT\t${sample}"
    sed "0,/${before}$/{s/${before}$/${after}/}" ${in_vcf} > ${in_vcf}.new && \
    mv ${in_vcf}.new ${in_vcf}
}


fixVCFsampleName HBCC_1385.10Xrna.vcf ${sample}
