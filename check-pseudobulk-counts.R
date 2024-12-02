#!/usr/bin/env Rscript


library(data.table)
library(foreach)

setwd('/data/CARD_singlecell/cortex_subtype/output')

get_data <- function(mode, celltype, cohort) {
    filename <- paste0(mode, '/', celltype, '/pseudobulk/', cohort, '_cpm_sum_pseudobulk_model_counts.csv')
    dat <- fread(filename)
    dat <- melt(dat, measure.vars=colnames(dat)[2:length(colnames(dat))], variable.name='sample', value.name='pseudobulk_cpm_sum')
    dat[, 'mode' := mode]
    dat[, 'celltype' := celltype]
    dat[, 'cohort' := cohort]
}

get_celltype <- function(celltype) {
    foreach(mode='atac', .combine='rbind') %do% {
        foreach(cohort=c('hbcc','nabec'), .combine='rbind') %do% {
            dat <- get_data(mode, celltype, cohort)
        }
    }
}

o <- get_celltype('VC')

library(ggplot2)

ggplot(o, aes(x=sample, y=pseudobulk_cpm_sum)) + geom_violin() + facet_grid(.~cohort, scales='free_x')

o.quantiles <- o[, list(quantile(pseudobulk_cpm_sum)), by=list(cohort,sample)]
o.quantiles[, q := c('min','q1','med','q3','max'), by=sample]
o.quantiles[, q := factor(q, levels=rev(c('min','q1','med','q3','max')))]

ggplot(o.quantiles, aes(x=sample, y=V1)) + geom_point() + facet_grid(q~cohort, scales='free')
