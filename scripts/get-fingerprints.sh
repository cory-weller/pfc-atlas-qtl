#!/usr/bin/env bash
#SBATCH --time 2:00:00
#SBATCH --mem 24g
#SBATCH --gres lscratch:100
#SBATCH --partition quick,norm
#SBATCH --cpus-per-task 2

STARTDIR='/data/CARD_singlecell/users/wellerca/pfc-atlas-qtl'


SAMPLEFILE=${STARTDIR}/samples.tsv

TMPDIR=/lscratch/${SLURM_JOB_ID}
HAPMAP="hg38_chr.reorder.map"



mkdir -p ${TMPDIR} && cd ${TMPDIR}

# Use positional argument only if not an array job
if [[ ! -z ${SLURM_ARRAY_TASK_ID} ]]; then
    N=${SLURM_ARRAY_TASK_ID}
fi

sample=$(sed -n ${N}p ${SAMPLEFILE} | cut -f 1)
sample=${sample/-ARC/}

COHORT=$(sed -n ${N}p ${SAMPLEFILE} | cut -f 11)

if [[ "${COHORT}" == 'HBCC' ]]; then
    sample=HBCC_${sample}
fi

bam=$(sed -n ${N}p ${SAMPLEFILE} | cut -f 12)
bfile="/data/CARD_singlecell/users/wellerca/pfc-atlas-qtl/data/genotypes/${COHORT}-forQTL"

# Create bam with sample ID forced to be $sample
module load samtools
if [[ "${COHORT}" == 'HBCC' ]]; then
    samtools addreplacerg -r "SM:${sample}" -r "ID:${sample}" -o ${sample}.bam ${bam}
else
    ln -s ${bam} ${sample}.bam
fi


module load plink/1.9

plink --bfile ${bfile} \
    --recode vcf-iid \
    --output-chr chrM \
    --keep-allele-order \
    --keep <(echo -e "${sample}\t${sample}") \
    --out ${sample} || { echo "Sample ${sample} not in ${bfile}"; exit 1; }



module load GATK
# Link reference files
ln -s ${STARTDIR}/hg38* .


# Get fingerprints from 10X RNA-seq VCF
gatk ExtractFingerprint  \
        --INPUT "${bam}" \
        --HAPLOTYPE_MAP ${HAPMAP} \
        --OUTPUT "${sample}.10Xrna.vcf" \
        --REFERENCE_SEQUENCE 'hg38.fa' || exit 1

# Fix misnamed 10Xrna.vcf by replacing 10
for file in *10Xrna.vcf; do
    mv ${file} ${file}.tmp
    sample=${file%.10Xrna.vcf}
    awk -v id=${sample} -vOFS='\t' '$1 ~ /^##/ {print}; $1 ~ /^#CHROM/ {print $1,$2,$3,$4,$5,$6,$7,$8,$9,id}; $1 ~ /chr/ {print}' ${file}.tmp > ${file}
    rm ${file}.tmp
done



# Recreate vcf index
rm ${sample}.10Xrna.vcf.idx
gatk IndexFeatureFile \
    -I "${sample}.10Xrna.vcf"

# Fix flipped ref/alt allele in genotypes VCF
module load R/4.3   # Requires data.table
Rscript ${STARTDIR}/scripts/reorder-vcf.R ${sample} # generates ${sample}.genotype-fingerprints-noheader.vcf

# Give the newly generated vcf header to the genotype vcf
grep '^##' "${sample}.10Xrna.vcf" | cat - ${sample}.genotype-fingerprints-noheader.vcf > "${sample}.genotype-fingerprints.vcf"

# Check fingerprints
gatk CheckFingerprint  \
        --INPUT "${sample}.10Xrna.vcf" \
        --GENOTYPES "${sample}.genotype-fingerprints.vcf" \
        --HAPLOTYPE_MAP ${HAPMAP} \
        --OUTPUT "${sample}.fingerprinting" || exit 1


zip ${sample}-fingerprinting.zip ${sample}.10Xrna.vcf ${sample}.genotype-fingerprints.vcf ${sample}.fingerprinting*

cp ${sample}-fingerprinting.zip ${STARTDIR}/fingerprints/