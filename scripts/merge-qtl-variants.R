#!/usr/bin/env Rscript
library(data.table)
library(foreach)

# Import data

nabec <- fread('qtl_analysis/nabec.dosage.vcf', skip='#CHROM')
hbcc <- fread('qtl_analysis//hbcc.dosage.vcf', skip='#CHROM')

# Rename columns

nabec.samples.old <- grep('^[UKS].*[0-9]$', colnames(nabec), value=T)
hbcc.samples.old <- grep('^HBCC.*[0-9]$', colnames(hbcc), value=T)

nabec.samples <- gsub('_[KSU].*[0-9]$', '', nabec.samples.old)
hbcc.samples <- gsub('_HBCC.*[0-9]$', '', hbcc.samples.old)

setnames(hbcc, hbcc.samples.old, hbcc.samples)
setnames(nabec, nabec.samples.old, nabec.samples)
setnames(hbcc, '#CHROM', 'CHROM')
setnames(nabec, '#CHROM', 'CHROM')

hbcc[, c('QUAL','FILTER','INFO','FORMAT') := NULL]
nabec[, c('QUAL','FILTER','INFO','FORMAT') := NULL]

setkey(hbcc, CHROM, POS, ID)
setkey(nabec, CHROM, POS, ID)

dat.merge <- merge(hbcc[,1:5], nabec[,1:5], all=T)

hbcc.long <- melt(hbcc, measure.vars=hbcc.samples, value.name='Dosage', variable.name='Sample')
nabec.long <- melt(nabec, measure.vars=nabec.samples, value.name='Dosage', variable.name='Sample')
combined.long <- rbindlist(list(hbcc.long, nabec.long))

fwrite(combined.long, file='data/alt-allele-dosage.tsv', sep='\t')