#!/usr/bin/env bash
#SBATCH --time 8:00:00
#SBATCH --mem 60g
#SBATCH --gres lscratch:100
#SBATCH --partition norm
#SBATCH --ntasks 12
#SBATCH --cpus-per-task 1


COHORT=${1}
STARTDIR='/data/CARD_singlecell/users/wellerca/pfc-atlas-qtl'


TMPDIR=/lscratch/${SLURM_JOB_ID}
HAPMAP="hg38_chr.reorder.map"


mkdir -p ${TMPDIR} && cd ${TMPDIR}



cp ${STARTDIR}/${HAPMAP} .

cp ${STARTDIR}/fingerprints/${COHORT}/*.zip .
for file in *.zip; do unzip $file; done
module load bcftools
for file in *.vcf; do bgzip ${file} && tabix -p vcf ${file}.gz; done

# Generate file list of vcfs to merge for bcftools argument
ls *rna.vcf.gz > files.txt
bcftools merge --file-list files.txt --output ${COHORT}-merged-rna.vcf
bgzip ${COHORT}-merged-rna.vcf && tabix -p vcf ${COHORT}-merged-rna.vcf.gz

module load GATK
all_samples=($(ls *10Xrna.vcf.gz | cat | cut -d '.' -f 1))
for sample in ${all_samples[@]}; do
    gatk CrosscheckFingerprints   \
            --INPUT "${sample}.genotype-fingerprints.vcf.gz" \
            --SECOND_INPUT ${COHORT}-merged-rna.vcf.gz \
            --HAPLOTYPE_MAP ${HAPMAP} \
            --CROSSCHECK_BY SAMPLE \
            --CROSSCHECK_MODE CHECK_ALL_OTHERS \
            --OUTPUT "${sample}.checkothers.txt"

    # gatk CrosscheckFingerprints   \
    #         --INPUT "${sample}.genotype-fingerprints.vcf.gz" \
    #         --SECOND_INPUT ${COHORT}-merged-rna.vcf.gz \
    #         --HAPLOTYPE_MAP ${HAPMAP} \
    #         --CROSSCHECK_BY SAMPLE \
    #         --CROSSCHECK_MODE CHECK_SAME_SAMPLE \
    #         --OUTPUT "${sample}.selfcheck.txt"

    awk 'NR>6' ${sample}.checkothers.txt | sed '/^$/d' > ${sample}.crosscheck.txt
done
