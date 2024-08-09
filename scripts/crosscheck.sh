#!/usr/bin/env bash
#SBATCH --time 1:30:00
#SBATCH --mem 60g
#SBATCH --partition quick
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4

module load GATK

gatk CrosscheckFingerprints   \
    --INPUT fingerprints/dna/dna-merged.vcf.gz \
    --SECOND_INPUT fingerprints/rna/rna-merged.vcf.gz \
    --HAPLOTYPE_MAP hg38_chr.reorder.map \
    --CROSSCHECK_BY SAMPLE \
    --CROSSCHECK_MODE CHECK_ALL_OTHERS \
    --OUTPUT fingerprints/crosscheck/all-pairs.txt

