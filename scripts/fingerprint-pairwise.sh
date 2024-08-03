#!/usr/bin/env bash
#SBATCH --time 8:00:00
#SBATCH --mem 60g
#SBATCH --gres lscratch:100
#SBATCH --partition norm
#SBATCH --ntasks 12
#SBATCH --cpus-per-task 1

NTHREADS=12
COHORT=${1}

STARTDIR=${PWD}
FINGERPRINTDIR="${STARTDIR}/fingerprints/${COHORT}"
HAPMAP="${STARTDIR}/hg38_chr.reorder.map"

# Work in local scratch directory
TMPDIR=/lscratch/${SLURM_JOB_ID}
mkdir -p ${TMPDIR} && cd ${TMPDIR}


# Unzip all .zip files containing RNA and WGS fingerprints
for file in $(find ${FINGERPRINTDIR}); do
    unzip ${file}
done

module load GATK

# Fix two-line fingerprint VCFs by removing -ARC
for file in *genotype-fingerprints.vcf; do
    mv ${file} ${file}.tmp
    grep -v '\-ARC' ${file}.tmp > ${file}
    rm ${file}.tmp
done

# Fix misnamed 10Xrna.vcf by replacing 10
for file in *10Xrna.vcf; do
    mv ${file} ${file}.tmp
    sample=${file%.10Xrna.vcf}
    awk -v id=${sample} -vOFS='\t' '$1 ~ /^##/ {print}; $1 ~ /^#CHROM/ {print $1,$2,$3,$4,$5,$6,$7,$8,$9,id}; $1 ~ /chr/ {print}' ${file}.tmp > ${file}
    rm ${file}.tmp
done

# files=($(ls *10Xrna.vcf | cut -d '.' -f 1))
# N=${#files[@]}
# let N=$N-1
# for i in $(seq 1 $N); do
#     for j in $(seq $i $N); do
#         sample1=${files["$i"]}
#         sample2=${files["$j"]}
#         rna_alias=$(grep '#CHROM' ${sample1}.10Xrna.vcf | awk '{print $10}')
#         geno_alias=$(grep '#CHROM' ${sample2}.genotype-fingerprints.vcf | awk '{print $10}')
#         gatk CheckFingerprint  \
#                 --INPUT "${sample1}.10Xrna.vcf" \
#                 --OBSERVED_SAMPLE_ALIAS ${sample1} \
#                 --EXPECTED_SAMPLE_ALIAS ${sample2} \
#                 --GENOTYPES "${sample2}.genotype-fingerprints.vcf" \
#                 --HAPLOTYPE_MAP ${HAPMAP} \
#                 --OUTPUT "${sample1}_x_${sample2}.crosscheck"
#     done
# done

function crosscheck() {
    local HAPMAP=${1}
    local sample1=${2}
    local sample2=${3}
    gatk CheckFingerprint  \
            --INPUT "${sample1}.10Xrna.vcf" \
            --OBSERVED_SAMPLE_ALIAS ${sample1} \
            --EXPECTED_SAMPLE_ALIAS ${sample2} \
            --GENOTYPES "${sample2}.genotype-fingerprints.vcf" \
            --HAPLOTYPE_MAP ${HAPMAP} \
            --OUTPUT "${sample1}_x_${sample2}.crosscheck"
}

export -f crosscheck
files=($(ls *10Xrna.vcf | cut -d '.' -f 1))

parallel -j ${NTHREADS} crosscheck {1} {2} {3} ::: ${HAPMAP} ::: ${files[@]} ::: ${files[@]}

mkdir ${COHORT}
mv *crosscheck* ${COHORT}
zip -r NABEC-crosscheck-fingerprints.zip ${COHORT}

mv NABEC-crosscheck-fingerprints.zip ${FINGERPRINTDIR}