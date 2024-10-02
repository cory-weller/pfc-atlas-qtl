#!/usr/bin/env Rscript

# R/4.3.2

library(data.table)

setwd('data')

tss <- fread('TSS_plusminus1000.bed')
setnames(tss, c('chr','start','stop','refseq'))
fwrite(tss, file='TSS-blacklist.bed', row.names=F, col.names=T, sep='\t', quote=F)

boyle <- fread('hg38-blacklist.v2.bed.gz', col.names=c('chr','start','stop','reason'))
fwrite(boyle, file='boyle-blacklist.bed', row.names=F, col.names=T, sep='\t', quote=F)

tss[, reason := paste0('TSS-', refseq)]
tss[, refseq := NULL]

both <- rbind(tss, boyle)
setkey(both, chr, start, stop)
fwrite(both, file='boyle-plus-TSS-blacklist.bed', row.names=F, col.names=T, sep='\t', quote=F)
