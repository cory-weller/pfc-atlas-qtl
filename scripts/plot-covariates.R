#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(ggthemes)
library(viridis)
library(ggrepel)

covs <- fread('samples.tsv')

nabec.pcs <- fread('genotypes/NABEC_polarized.pruned_pc.txt')
hbcc.pcs <- fread('genotypes/HBCC_polarized.pruned_pc.txt')

hbcc.covs <- covs[cohort == 'HBCC']
hbcc.covs[, sample := gsub('-ARC$', '', sample,)]
hbcc.covs[, sample := paste0('HBCC_', sample)]
nabec.covs <- covs[cohort == 'NABEC']
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
    dat[PC1 > 0.2 | PC2 > 0.25 | PC2 < -0.5, 'txt' := IID]


    g <- ggplot(dat, aes(x=PC1, y=PC2, color=.data[[.cov]])) +
        geom_point() +
        theme_few()
    if(! is.factor(dat[[.cov]])) {
        g <- g + scale_color_viridis()
    }
    # Save SVGs
    outfile <- paste0('PCA-plots/', toupper(.cohort), '-', .cov, '.svg')
    if(.cohort == 'HBCC') {
        g <- g + geom_text_repel(data=dat, size =1, segment.size=0.2, color='black', max.overlaps=20, min.segment.length = unit(0, 'lines'), aes(label=txt))
    }
    ggsave(g, file=outfile, width=12, heigh=12, units='cm')

    # Save PNGs
    outfile <- paste0('PCA-plots/', toupper(.cohort), '-', .cov, '.png')
    if(.cohort == 'HBCC') {
        g <- g + geom_text_repel(data=dat, size =1, segment.size=0.2, color='black', max.overlaps=20, min.segment.length = unit(0, 'lines'), aes(label=txt))
    }
    ggsave(g, file=outfile, width=12, heigh=12, units='cm')
}

for(i in covs_to_plot) {
    for(cohort in c('NABEC','HBCC')) {
        plot_cov(pcs, cohort, i)
    }
}