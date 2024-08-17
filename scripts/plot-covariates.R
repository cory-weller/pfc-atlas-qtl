#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(ggthemes)
library(viridis)

covs <- fread('data/covariates.csv')

nabec.pcs <- fread('data/genotypes/NABEC.pruned_pc.txt')
hbcc.pcs <- fread('data/genotypes/HBCC.pruned_pc.txt')

hbcc.covs <- covs[cohort == 'hbcc']
hbcc.covs[, sample := gsub('-ARC$', '', sample,)]
hbcc.covs[, sample := paste0('HBCC_', sample)]
nabec.covs <- covs[cohort == 'nabec']
nabec.covs[, sample := gsub('-ARC$', '', sample,)]

nabec <- merge(nabec.covs, nabec.pcs, by.x='sample', by.y='FID')
hbcc <- merge(hbcc.covs, hbcc.pcs, by.x='sample', by.y='FID')

pcs <- rbindlist(list(nabec, hbcc))

covs_to_plot <- c('Sex','Age','PMI','Ancestry','Brain_bank','Sequencing','batch','LibraryPrep','Homogenization')

pcs[, Sex := factor(Sex)]
pcs[, Ancestry := factor(Ancestry)]
pcs[, Brain_bank := factor(Brain_bank)]
pcs[, Sequencing := factor(Sequencing)]
pcs[, batch := factor(batch)]
pcs[, Homogenization := factor(Homogenization)]
pcs[, LibraryPrep := factor(LibraryPrep)]



plot_cov <- function(DT, .cohort, .cov) {
    dat <- DT[cohort==.cohort]

    g <- ggplot(dat, aes(x=PC1, y=PC2, color=.data[[.cov]])) +
        geom_point() +
        theme_few()
    if(! is.factor(dat[[.cov]])) {
        g <- g + scale_color_viridis()
    }
    outfile <- paste0('plots/', toupper(.cohort), '-', .cov, '.png')
    ggsave(g, file=outfile, width=12, heigh=12, units='cm')
}

for(i in covs_to_plot) {
    for(cohort in c('nabec','hbcc')) {
        plot_cov(pcs, cohort, i)
    }
}