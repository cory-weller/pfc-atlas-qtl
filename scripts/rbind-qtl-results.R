#!/usr/bin/env Rscript

library(data.table)
library(foreach)

o <- foreach(file=list.files('output', pattern='*assoc.txt.gz', full.names=TRUE), .combine='rbind') %do% {
    dat <- fread(file)
    split_name <- unlist(strsplit(basename(file), split='_|\\.'))
    cohort <- split_name[1]
    celltype <- split_name[2]
    mode <- split_name[3]
    cbind(dat, cohort, celltype, mode)
}

o[, c('variant_chr','variant_pos','variant_ref','variant_alt') := tstrsplit(variant_id, split=':')]
o[, variant_pos := as.numeric(variant_pos)]

# Exclude 4 NA rows
o <- o[!is.na(variant_id)][variant_id != '.']

fwrite(o, file='cis-QTL-combined.tsv', row.names=F, col.names=T, sep='\t', quote=F)
fwrite(o[mode=='rna'], file='cis-eQTL.tsv', row.names=F, col.names=T, sep='\t', quote=F)
fwrite(o[mode=='atac'], file='cis-caQTL.tsv', row.names=F, col.names=T, sep='\t', quote=F)