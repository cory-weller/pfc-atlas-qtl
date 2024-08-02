#!/usr/bin/env bash
#SBATCH --time 1:00:00
#SBATCH --mem 24g
#SBATCH --gres lscratch:100
#SBATCH --partition quick,norm
#SBATCH --cpus-per-task 2

COHORT=${1}

STARTDIR='/data/CARD_singlecell/users/wellerca/pfc-atlas-qtl'

bamlist="${STARTDIR}/${COHORT}_bam_file_paths.txt"


# Retrieve and reorder haplotype map

if [[ ! -f 'hg38_chr.map' ]]; then
    wget -P fingerprints https://github.com/naumanjaved/fingerprint_maps/raw/master/map_files/hg38_chr.map
fi

# Build new map
function reorderMap() {
    local fn=${1}
    local out_fn=${fn/map/reorder.map}
    echo "$fn being reordered to $out_fn"
    grep '^@SQ' hg38_chr.map  > header.txt
    grep 'chr1\s' header.txt > ${out_fn}
    grep 'chr10\s' header.txt >> ${out_fn}
    grep 'chr11\s' header.txt >> ${out_fn}
    grep 'chr12\s' header.txt >> ${out_fn}
    grep 'chr13\s' header.txt >> ${out_fn}
    grep 'chr14\s' header.txt >> ${out_fn}
    grep 'chr15\s' header.txt >> ${out_fn}
    grep 'chr16\s' header.txt >> ${out_fn}
    grep 'chr17\s' header.txt >> ${out_fn}
    grep 'chr18\s' header.txt >> ${out_fn}
    grep 'chr19\s' header.txt >> ${out_fn}
    grep 'chr2\s' header.txt >> ${out_fn}
    grep 'chr20\s' header.txt >> ${out_fn}
    grep 'chr21\s' header.txt >> ${out_fn}
    grep 'chr22\s' header.txt >> ${out_fn}
    grep 'chr3\s' header.txt >> ${out_fn}
    grep 'chr4\s' header.txt >> ${out_fn}
    grep 'chr5\s' header.txt >> ${out_fn}
    grep 'chr6\s' header.txt >> ${out_fn}
    grep 'chr7\s' header.txt >> ${out_fn}
    grep 'chr8\s' header.txt >> ${out_fn}
    grep 'chr9\s' header.txt >> ${out_fn}
    grep '#' hg38_chr.map >> ${out_fn}
    grep '^chr1\s' hg38_chr.map >> ${out_fn}
    grep '^chr10\s' hg38_chr.map >> ${out_fn}
    grep '^chr11\s' hg38_chr.map >> ${out_fn}
    grep '^chr12\s' hg38_chr.map >> ${out_fn}
    grep '^chr13\s' hg38_chr.map >> ${out_fn}
    grep '^chr14\s' hg38_chr.map >> ${out_fn}
    grep '^chr15\s' hg38_chr.map >> ${out_fn}
    grep '^chr16\s' hg38_chr.map >> ${out_fn}
    grep '^chr17\s' hg38_chr.map >> ${out_fn}
    grep '^chr18\s' hg38_chr.map >> ${out_fn}
    grep '^chr19\s' hg38_chr.map >> ${out_fn}
    grep '^chr2\s' hg38_chr.map >> ${out_fn}
    grep '^chr20\s' hg38_chr.map >> ${out_fn}
    grep '^chr21\s' hg38_chr.map >> ${out_fn}
    grep '^chr22\s' hg38_chr.map >> ${out_fn}
    grep '^chr3\s' hg38_chr.map >> ${out_fn}
    grep '^chr4\s' hg38_chr.map >> ${out_fn}
    grep '^chr5\s' hg38_chr.map >> ${out_fn}
    grep '^chr6\s' hg38_chr.map >> ${out_fn}
    grep '^chr7\s' hg38_chr.map >> ${out_fn}
    grep '^chr8\s' hg38_chr.map >> ${out_fn}
    grep '^chr9\s' hg38_chr.map >> ${out_fn}
}

if [[ ! -f 'fingerprints/hg38_chr.reorder.map' ]]; then
    reorderMap fingerprints/hg38_chr.map
fi



SAMPLEFILE=samples.txt
TMPDIR=/lscratch/${SLURM_JOB_ID}
HAPMAP="${STARTDIR}/hg38_chr.reorder.map"

mkdir -p ${TMPDIR} && cd ${TMPDIR}

paste <(sed 's@/.*Multiome/@@g' ${bamlist}  | sed 's@/.*.bam$@@g') ${bamlist} > samples.txt

N=${SLURM_ARRAY_TASK_ID}

sample=$(sed -n ${N}p ${SAMPLEFILE} | cut -f 1)
sample=${sample/-ARC/}
bam=$(sed -n ${N}p ${SAMPLEFILE} | cut -f 2)
bfile="/data/CARD_singlecell/users/wellerca/pfc-atlas-qtl/data/genotypes/${COHORT}-forQTL"

cp ~/pfc-atlas-qtl/hg38_chr.reorder.map .

module load plink/1.9

plink --bfile ${bfile} \
    --recode vcf-iid \
    --output-chr chrM \
    --keep-allele-order \
    --keep <(echo -e "${sample}\t${sample}") \
    --out ${sample} || exit 1





module load GATK
REF='/fdb/cellranger-arc/refdata-cellranger-arc-GRCh38-2020-A/fasta/genome.fa'
ln -s ${REF} hg38.fa


gatk CreateSequenceDictionary \
    R=hg38.fa 

ml samtools
samtools faidx hg38.fa

# Get fingerprints from 10X RNA-seq VCF
gatk ExtractFingerprint  \
        --INPUT "${bam}" \
        --HAPLOTYPE_MAP ${HAPMAP} \
        --OUTPUT "${sample}.10Xrna.vcf" \
        --REFERENCE_SEQUENCE 'hg38.fa' || exit 1

# Fix flipped ref/alt allele in genotypes VCF
module load R/4.3   # Requires data.table
Rscript ${STARTDIR}/scripts/reorder-vcf.R ${sample} # generates ${sample}.genotype-fingerprints-noheader.vcf

# Give the newly generated vcf header to the genotype vcf
grep '^#' "${sample}.10Xrna.vcf" | cat - ${sample}.genotype-fingerprints-noheader.vcf > "${sample}.genotype-fingerprints.vcf"

# Check fingerprints
gatk CheckFingerprint  \
        --INPUT "${sample}.10Xrna.vcf" \
        --GENOTYPES "${sample}.genotype-fingerprints.vcf" \
        --HAPLOTYPE_MAP ${HAPMAP} \
        --OUTPUT "${sample}.fingerprinting" || exit 1

zip ${sample}-fingerprinting.zip ${sample}.10Xrna.vcf ${sample}.genotype-fingerprints.vcf ${sample}.fingerprinting*

cp ${sample}-fingerprinting.zip ${STARTDIR}/fingerprints