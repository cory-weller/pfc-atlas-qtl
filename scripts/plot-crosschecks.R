#!/usr/bin/env Rscript

library(data.table)
library(foreach)
library(ggplot2)
library(ggthemes)
library(viridis)

dat <- fread('crosscheck-results.txt')

dat[RIGHT_SAMPLE %like% 'HBCC', 'cohort' := 'HBCC']
dat[! RIGHT_SAMPLE %like% 'HBCC', 'cohort' := 'NABEC']

NABEC <- dat[cohort == 'NABEC']
HBCC <- dat[cohort == 'HBCC']

NABEC[RIGHT_GROUP_VALUE %like% 'UMARY|KEN', c('BrainBank','ID_N') := tstrsplit(RIGHT_GROUP_VALUE, split='-')]
NABEC[RIGHT_GROUP_VALUE %like% 'SH', c('BrainBank','ID_N1','ID_N2') := tstrsplit(RIGHT_GROUP_VALUE, split='-')]
NABEC[RIGHT_GROUP_VALUE %like% 'SH', 'ID_N' := paste0(ID_N1, ID_N2)]
NABEC[, ID_N := as.numeric(ID_N)]
NABEC[, c('ID_N1','ID_N2') := NULL]
setkey(NABEC, BrainBank, ID_N)
NABEC.samples <- unique(NABEC$RIGHT_GROUP_VALUE)

HBCC[, 'BrainBank' := 'HBCC']
HBCC[, 'ID_N' := tstrsplit(RIGHT_GROUP_VALUE, split='_')[2]]
HBCC[, ID_N := as.numeric(ID_N)]
setkey(HBCC, ID_N)
HBCC.samples <- unique(HBCC$RIGHT_GROUP_VALUE)

samples <- c(NABEC.samples, HBCC.samples)

dat <- rbindlist(list(NABEC, HBCC))

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


ggsave(g, file=paste0('~/pfc-atlas-qtl/fingerprints/crosscheck.png'), width=10+n_x/3, height=n_y/3, units='cm', limitsize=FALSE)

