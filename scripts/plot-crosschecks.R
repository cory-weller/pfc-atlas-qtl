#!/usr/bin/env Rscript

library(data.table)
library(foreach)
library(ggplot2)
library(ggthemes)
library(viridis)

cohort <- commandArgs(trailingOnly=TRUE)[1]

dat <- foreach(file=list.files(pattern='*crosscheck*'), .combine='rbind') %do% {
    fread(file)
}

dat[, 'ID_N' := tstrsplit(RIGHT_GROUP_VALUE, split='_')[2]]
dat[, ID_N := as.numeric(ID_N)]
setkey(dat, ID_N)
samples <- unique(dat$RIGHT_GROUP_VALUE)
dat[, ID_N := NULL]
cols_to_remove <- c('_RUN_BARCODE','_LANE','_MOLECULAR_BARCODE_SEQUENCE','_LIBRARY')
dat[, paste0('LEFT', cols_to_remove) := NULL]
dat[, paste0('RIGHT', cols_to_remove) := NULL]

dat[, LEFT_GROUP_VALUE := factor(LEFT_GROUP_VALUE, levels=samples)]
dat[, RIGHT_GROUP_VALUE := factor(RIGHT_GROUP_VALUE, levels=samples)]
fwrite(dat, file=paste0('~/pfc-atlas-qtl/fingerprints/',cohort,'-crosscheck.tsv'), quote=F, row.names=F, col.names=T, sep='\t')
n_x <- length(unique(dat[, LEFT_GROUP_VALUE]))
n_y <- length(unique(dat[, RIGHT_GROUP_VALUE]))
g <- ggplot(dat, aes(x=LEFT_GROUP_VALUE, y=RIGHT_GROUP_VALUE, fill=LOD_SCORE)) +
    geom_tile() +
    theme_few() +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    scale_fill_viridis() +
    labs(x='Sample Genotype', y='Sample snRNA')


ggsave(g, file=paste0('~/pfc-atlas-qtl/fingerprints/',cohort,'-crosscheck.png'), width=10+n_x/2, height=n_y/2, units='cm')

