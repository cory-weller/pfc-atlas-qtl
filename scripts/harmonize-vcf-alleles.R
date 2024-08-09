#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
FILESTEM <- args[1]
RNAVCF <- args[2]

library(data.table)
vcf <- fread('genotypes.common.vcf', header=T)
samples <- colnames(vcf[,-c(1:9)])
setnames(vcf, '#CHROM','CHROM')

rnaloci <- fread(cmd=paste0('zcat ', RNAVCF, " | awk '$1 ~ /^chr/ {print $1,$2,$4,$5}'"))
#rnaloci <- fread(cmd="awk '$1 ~ /^chr/ {print $1,$2,$4,$5}' hg38_chr.reorder.map")
# Commented out line was only to test cases where flipping is necessary

setnames(rnaloci, c('CHROM','POS','REF','ALT'))
chrom.order <- unique(rnaloci[,CHROM])

vcf.info <- vcf[, c('CHROM','POS','REF','ALT')]
setkey(vcf.info, CHROM, POS)
setkey(rnaloci, CHROM, POS)
setkey(vcf, CHROM, POS)
dat <- merge(vcf.info, rnaloci)
setnames(dat, 'REF.y', 'REF.rna')
setnames(dat, 'ALT.y', 'ALT.rna')
setnames(dat, 'REF.x','REF.dna')
setnames(dat, 'ALT.x', 'ALT.dna')

dat[REF.rna == REF.dna & ALT.dna == ALT.rna, flip := FALSE]

dat[REF.rna == ALT.dna & ALT.rna == REF.dna, flip := TRUE]
dat <- dat[! is.na(flip)]   # Exclude other variants that don't line up  (likely indels)
# 1/1 to #/#
# 0/0 to 1/1
# #/# to 0/0
vcf <- vcf[dat] # Merge in REF/ALT/flip cols
vcf.keeporder <- copy(vcf[flip==FALSE])
vcf.reorder <- copy(vcf[flip==TRUE])
vcf.reorder.infocols <- copy(vcf.reorder[, 1:9])
vcf.reorder.infocols.order <- copy(colnames(vcf.reorder.infocols))
setnames(vcf.reorder.infocols, c('REF','ALT'), c('ALT','REF'))

vcf.reorder <- vcf.reorder[, lapply(.SD, gsub, pattern='0', replacement='#'), .SDcols=samples]
vcf.reorder <- vcf.reorder[, lapply(.SD, gsub, pattern='1', replacement='0'), .SDcols=samples]
vcf.reorder <- vcf.reorder[, lapply(.SD, gsub, pattern='#', replacement='1'), .SDcols=samples]
vcf.reorder <- cbind(vcf.reorder.infocols, vcf.reorder)
setcolorder(vcf.reorder, vcf.reorder.infocols.order)

vcf.reorder[, ID := paste(CHROM, POS, REF, ALT, sep=':')]
vcf.final <- rbindlist(list(vcf.reorder, vcf.keeporder[, .SD, .SDcols=colnames(vcf.reorder)]))
vcf.final[, CHROM := factor(CHROM, levels=chrom.order)]

setkey(vcf.final, CHROM, POS)
setnames(vcf.final, 'CHROM', '#CHROM')
fwrite(vcf.final, file=paste0(FILESTEM, '.genotype-fingerprints.noheader.vcf'), quote=F, row.names=F, col.names=T, sep='\t')