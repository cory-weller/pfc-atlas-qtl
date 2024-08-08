#!/usr/bin/env bash
#SBATCH --time 8:00:00
#SBATCH --mem 60g
#SBATCH --gres lscratch:100
#SBATCH --partition norm
#SBATCH --ntasks 12
#SBATCH --cpus-per-task 1

STARTDIR='/data/CARD_singlecell/users/wellerca/pfc-atlas-qtl'


TMPDIR=/lscratch/${SLURM_JOB_ID}
#TMPDIR=/dev/shm/${SLURM_JOB_ID}
HAPMAP="hg38_chr.reorder.map"


mkdir -p ${TMPDIR} && cd ${TMPDIR}



cp ${STARTDIR}/${HAPMAP} .

if [[ ! -f ${STARTDIR}/fingerprints/fingerprints-merged-rna.vcf ]]; then
    cp ${STARTDIR}/fingerprints/*.zip .
    for file in *.zip; do unzip $file; done
    module load bcftools
    for file in *.vcf; do bgzip ${file} && tabix -p vcf ${file}.gz; done

    # Generate file list of vcfs to merge for bcftools argument
    ls *rna.vcf.gz > files.txt
    bcftools merge --file-list files.txt --output fingerprints-merged-rna.vcf
    bgzip fingerprints-merged-rna.vcf && tabix -p vcf fingerprints-merged-rna.vcf.gz
else
    cp ${STARTDIR}/fingerprints/fingerprints-merged-rna.vcf* .
fi

# module load GATK
# all_samples=($(ls *10Xrna.vcf.gz | cat | cut -d '.' -f 1))
# for sample in ${all_samples[@]}; do
#     gatk CrosscheckFingerprints   \
#             --INPUT "${sample}.genotype-fingerprints.vcf.gz" \
#             --SECOND_INPUT fingerprints-merged-rna.vcf.gz \
#             --HAPLOTYPE_MAP ${HAPMAP} \
#             --CROSSCHECK_BY SAMPLE \
#             --CROSSCHECK_MODE CHECK_ALL_OTHERS \
#             --OUTPUT "${sample}.checkothers.txt"

#     awk 'NR>6' ${sample}.checkothers.txt | sed '/^$/d' > ${sample}.crosscheck.txt
# done


function crosscheck() {
    local sample=${1}
    gatk CrosscheckFingerprints   \
        --INPUT "${sample}.genotype-fingerprints.vcf.gz" \
        --SECOND_INPUT fingerprints-merged-rna.vcf.gz \
        --HAPLOTYPE_MAP hg38_chr.reorder.map \
        --CROSSCHECK_BY SAMPLE \
        --CROSSCHECK_MODE CHECK_ALL_OTHERS \
        --OUTPUT "${sample}.checkothers.txt"

    awk 'NR>7' ${sample}.checkothers.txt | sed '/^$/d' > ${sample}.crosscheck.txt
}

export -f crosscheck

module load GATK
parallel -j 12 crosscheck {1} ::: ${all_samples[@]}

echo -e 'LEFT_GROUP_VALUE\tRIGHT_GROUP_VALUE\tRESULT\tDATA_TYPE\tLOD_SCORE\tLOD_SCORE_TUMOR_NORMAL\tLOD_SCORE_NORMAL_TUMOR\tLEFT_RUN_BARCODE\tLEFT_LANE\tLEFT_MOLECULAR_BARCODE_SEQUENCE\tLEFT_LIBRARY\tLEFT_SAMPLE\tLEFT_FILE\tRIGHT_RUN_BARCODE\tRIGHT_LANE\tRIGHT_MOLECULAR_BARCODE_SEQUENCE\tRIGHT_LIBRARY\tRIGHT_SAMPLE\tRIGHT_FILE' > cross-check_results.txt
cat *crosscheck* >> cross-check_results.txt
cp cross-check_results.txt ${STARTDIR}/fingerprints/crosscheck-results.txt

