#!/usr/bin/env bash
#SBATCH --time 1:00:00
#SBATCH --mem 50g
#SBATCH --gres lscratch:100
#SBATCH --partition quick,norm
#SBATCH --cpus-per-task 2

module load GATK/4.6
module load samtools/1.19
module load R/4.3

if [[ ${1} == '' ]]; then
    echo "ERROR: sample file required as first (and only) argument!"
    exit 1
fi

REFDIR="/data/${USER}/.hg38"
STARTDIR=$(realpath ${PWD})
SAMPLEFILE=$(realpath ${1})

for file in ${REFDIR}/hg38{_chr.map,_chr.reorder.map,.dict,.fa,.fa.fai}; do
    MISSING=0
    if [[ ! -f ${file} ]]; then
        echo ${file} missing! Did you run prepare-refdata.sh ?
        let MISSING=${MISSING}+1
    fi
done
if [[ ${MISSING} > 0 ]]; then
    echo "Exiting due to missing ref data"; exit 1 
fi



TMPDIR=/lscratch/${SLURM_JOB_ID}
mkdir -p ${TMPDIR}
ln -s ${REFDIR}/hg38* ${TMPDIR}
cd ${TMPDIR}


# Use positional argument only if not an array job
if [[ ! -z ${SLURM_ARRAY_TASK_ID} ]]; then
    N=${SLURM_ARRAY_TASK_ID}
fi

sample=$(sed -n ${N}p ${SAMPLEFILE} | awk '{print $1}')
bam=$(sed -n ${N}p ${SAMPLEFILE} | awk '{print $2}')


# Get fingerprints from 10X RNA-seq VCF
gatk ExtractFingerprint  \
        --INPUT "${bam}" \
        --HAPLOTYPE_MAP hg38_chr.reorder.map \
        --OUTPUT "${sample}.10Xrna.vcf" \
        --REFERENCE_SEQUENCE 'hg38.fa' || exit 1

# Fix misnamed 10Xrna.vcf by replacing field 10]
file="${sample}.10Xrna.vcf"
mv ${file} ${file}.tmp
awk -v id=${sample} -vOFS='\t' '$1 ~ /^##/ {print}; $1 ~ /^#CHROM/ {print $1,$2,$3,$4,$5,$6,$7,$8,$9,id}; $1 ~ /chr/ {print}' ${file}.tmp > ${file}
rm ${file}.tmp

bgzip ${file} && tabix -p vcf ${file}.gz

# Export fingerprints vcf and tabix index to permanent directory
mv ${file}.gz ${STARTDIR}/fingerprints/rna
mv ${file}.gz.tbi ${STARTDIR}/fingerprints/rna


