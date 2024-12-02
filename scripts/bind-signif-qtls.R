#!/usr/bin/env Rscript

library(data.table)
library(foreach)


# Do mean


# Get interaction ones
o <- foreach(method=c('mean','sum'), .combine='rbind') %do% {
    files <- list.files(paste0('/data/CARD_singlecell/users/wellerca/pfc-atlas-qtl/QTL-output/',method, '-nominal-20241124'), full.names=T, recursive=T, pattern='*.tsv.gz')
    interaction_files <- grep('Age', files, value=T)
    noninteraction_files <- grep('Age', invert=T, files, value=T)
    foreach(file=interaction_files, .combine='rbind') %do% {
        splitname <- unlist(strsplit(file, split='/'))[9]
        splitdir <- unlist(strsplit(splitname, split='-'))
        cohort <- splitdir[1]
        celltype <- splitdir[2]
        mode <- splitdir[3]
        interaction <- splitdir[5]
        dat <- fread(file)
        dat <- dat[pval_gi_bh < 0.05]
        dat[, 'cohort' := cohort]
        dat[, 'celltype' := celltype]
        dat[, 'mode' := mode]
        dat[, 'interaction' := interaction]
        dat[, 'pseudobulk_method' := method]
        return(dat)
    }
}

fwrite(o, file='/data/CARD_singlecell/users/wellerca/pfc-atlas-qtl/significant-interaction-qtls.tsv.gz', quote=F, row.names=F, col.names=T, sep='\t')
rm(o)

# Get noninteraction files
# Get interaction ones
o <- foreach(method=c('mean','sum'), .combine='rbind') %do% {
    files <- list.files(paste0('/data/CARD_singlecell/users/wellerca/pfc-atlas-qtl/QTL-output/',method, '-nominal-20241124'), full.names=T, recursive=T, pattern='*.tsv.gz')
    interaction_files <- grep('Age', files, value=T)
    noninteraction_files <- grep('Age', invert=T, files, value=T)
    foreach(file=noninteraction_files, .combine='rbind') %do% {
        splitname <- unlist(strsplit(file, split='/'))[9]
        splitdir <- unlist(strsplit(splitname, split='-'))
        cohort <- splitdir[1]
        celltype <- splitdir[2]
        mode <- splitdir[3]
        interaction <- splitdir[5]
        dat <- fread(file)
        dat <- dat[pval_nominal_bh < 0.05]
        dat[, 'cohort' := cohort]
        dat[, 'celltype' := celltype]
        dat[, 'mode' := mode]
        dat[, 'interaction' := interaction]
        dat[, 'pseudobulk_method' := method]
        return(dat)
    }
}

fwrite(o, file='/data/CARD_singlecell/users/wellerca/pfc-atlas-qtl/significant-noninteraction-qtls.tsv.gz', quote=F, row.names=F, col.names=T, sep='\t')
rm(o)