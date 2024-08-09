#!/usr/bin/env bash
#SBATCH --time 1:00:00
#SBATCH --mem 100g
#SBATCH --gres lscratch:100
#SBATCH --partition quick,norm
#SBATCH --cpus-per-task 4

TMPDIR=/lscratch/${SLURM_JOB_ID}
REFDIR="/data/${USER}/.hg38"
OUTDIR=$(realpath fingerprints/dna)
RNAVCF=$(realpath fingerprints/rna/rna-merged.vcf.gz)
STARTDIR=$(realpath ${PWD})
mkdir -p ${TMPDIR}

ln -s ${REFDIR}/hg38* ${TMPDIR}
cd ${TMPDIR}



echo '/data/CARD_singlecell/users/wellerca/pfc-atlas-qtl/data/genotypes/NABEC-forQTL' > bfiles.txt
echo '/data/CARD_singlecell/users/wellerca/pfc-atlas-qtl/data/genotypes/HBCC' >> bfiles.txt





module load plink/1.9
module load R/4.3

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
    zcat ${RNAVCF} | awk '$1 ~ /^chr/ {print $1,$2}' > map.loci

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

    Rscript ${STARTDIR}/scripts/harmonize-vcf-alleles.R ${FILESTEM} ${RNAVCF}  # Creates ${FILESTEM}.genotype.fingerprints.noheader.vcf
    zgrep '^##' ${RNAVCF} > ${FILESTEM}.genotype-fingerprints.vcf
    cat ${FILESTEM}.genotype-fingerprints.noheader.vcf >> ${FILESTEM}.genotype-fingerprints.vcf
    rm indices.txt vcf.loci map.loci genotypes.common.vcf
done < bfiles.txt

module load samtools

Rscript - <<EOF
library(data.table)
library(foreach)

files <- list.files(pattern='*genotype-fingerprints.vcf$')

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
EOF

zgrep '^##' ${RNAVCF} > dna-merged.vcf
cat merged-genotype-noheader.vcf >> dna-merged.vcf
bgzip dna-merged.vcf && tabix -p vcf dna-merged.vcf.gz

mv dna-merged.vcf.gz ${OUTDIR}
mv dna-merged.vcf.gz.tbi ${OUTDIR}