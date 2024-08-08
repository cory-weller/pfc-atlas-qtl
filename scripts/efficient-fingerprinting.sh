#!/usr/bin/env bash
module load GATK
module load samtools
module load R/4.3


# VARS
REFDIR="/data/${USER}/.hg38"

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
cp fingerprints/fingerprints-merged-rna.vcf.gz ${TMPDIR}
ln -s ${REFDIR}/hg38* ${TMPDIR}
cd ${TMPDIR}



echo '/data/CARD_singlecell/users/wellerca/pfc-atlas-qtl/data/genotypes/NABEC-forQTL' > bfiles.txt
echo '/data/CARD_singlecell/users/wellerca/pfc-atlas-qtl/data/genotypes/HBCC' >> bfiles.txt





module load plink/1.9

# Extract genotype vcf information from plink files
while read BFILE; do
    FILESTEM=$(basename ${BFILE})

    # Get GENOTYPE fingerprints from plink bfile.bed
    plink --bfile ${BFILE} \
        --recode vcf-iid \
        --output-chr chrM \
        --keep-allele-order \
        --out ${FILESTEM}.genotypes



    header_length=$(grep '^#' ${FILESTEM}.genotypes.vcf | wc -l)

    awk '$1 ~ /^chr/ {print $1,$2}' ${FILESTEM}.genotypes.vcf > vcf.loci
    zcat fingerprints-merged-rna.vcf.gz | awk '$1 ~ /^chr/ {print $1,$2}' > map.loci

    Rscript - <<EOF
library(data.table)
vcf.loci <- fread('vcf.loci')
vcf.loci[, idx := 1:.N]
map.loci <- fread('map.loci')
setkey(vcf.loci, V1, V2)
setkey(map.loci, V1, V2)
dat <- merge(vcf.loci, map.loci)
indices <- as.character(dat[order(idx), idx])
writeLines(indices, con='indices.txt')
EOF


    # Extract desired rows from genotype vcf
    grep '^#CHROM' ${FILESTEM}.genotypes.vcf > genotypes.common.vcf
    awk 'FNR==NR{h[$1]; c++; next} FNR in h{print; if (!--c) exit}' \
        <(awk -v s=${header_length} '{print $1+s}' indices.txt) \
        ${FILESTEM}.genotypes.vcf >> genotypes.common.vcf

    # From here, reorder ref and alt alleles across all samples

    Rscript ~/pfc-atlas-qtl/scripts/harmonize-vcf-alleles.R ${FILESTEM}   # Creates ${FILESTEM}.genotype.fingerprints.noheader.vcf
    zgrep '^##' fingerprints-merged-rna.vcf.gz > ${FILESTEM}.genotype-fingerprints.vcf
    cat ${FILESTEM}.genotype-fingerprints.noheader.vcf >> ${FILESTEM}.genotype-fingerprints.vcf
    rm indices.txt vcf.loci map.loci genotypes.common.vcf
    #bgzip ${FILESTEM}.genotype-fingerprints.vcf && tabix -p vcf ${FILESTEM}.genotype-fingerprints.vcf.gz
done < bfiles.txt

bcftools merge -o merged-genos.vcf *genotype-fingerprints.vcf.gz 

```
R
library(data.table)
library(foreach)

files <- list.files(pattern='*genotype-fingerprints.vcf')

o <- foreach(file=files) %do% {
    fread(file, skip='#CHROM')
}

for(dt in o) {
    setkey(dt, '#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT')
}

fillNA <- function(DT, x) {
    # replaces all NA values with x in data.table DT
    for (j in seq_len(ncol(DT)))
        set(DT,which(is.na(DT[[j]])),j,x)
}

out <- do.call(merge, c(o, all=TRUE))
fillNA(out, './.')

fwrite(out, file='merged-genotype-noheader.vcf', quote=F, row.names=F, col.names=T, sep='\t')
```


zgrep '^##' fingerprints-merged-rna.vcf.gz > merged.genotype-fingerprints.vcf
cat merged-genotype-noheader.vcf >> merged.genotype-fingerprints.vcf
bgzip merged.genotype-fingerprints.vcf && tabix -p vcf merged.genotype-fingerprints.vcf.gz

module load GATK

gatk CrosscheckFingerprints   \
    --INPUT merged.genotype-fingerprints.vcf.gz \
    --SECOND_INPUT fingerprints-merged-rna.vcf.gz \
    --HAPLOTYPE_MAP hg38_chr.reorder.map \
    --CROSSCHECK_BY SAMPLE \
    --CROSSCHECK_MODE CHECK_ALL_OTHERS \
    --OUTPUT "all-crosscheck.txt"

























echo "Done generating genotype fingerprints for ${FILESTEM}"
exit 0

# Get RNA fingerprints from 10X bam files and merge
module load GATK
gatk ExtractFingerprint  \
        --INPUT "${bam}" \
        --HAPLOTYPE_MAP hg38_chr.reorder.map \
        --OUTPUT "${sample}.10Xrna.vcf" \
        --REFERENCE_SEQUENCE 'hg38.fa'

bgzip and index
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
