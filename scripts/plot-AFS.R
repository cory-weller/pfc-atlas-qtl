#!/usr/bin/env bash

library(data.table)
library(ggplot2)
library(ggthemes)

setwd('/vf/users/CARD_singlecell/users/wellerca/pfc-atlas-qtl/genotypes')

dat <- fread('HBCC_autosomes.traw')

sample_cols <- colnames(dat)[-c(1,2,3,4,5,6)]
n_samples <- length(sample_cols)

dat[, 'alt_count' := apply(.SD, 1, function(x) sum(x, na.rm=T)), .SDcols=sample_cols]
dat[, 'na_samples' := apply(.SD, 1, function(x) sum(is.na(x))), .SDcols=sample_cols]

dat[, total_alleles := 2*n_samples - 2*na_samples]
dat[, ref_count := total_alleles - alt_count]
dat[, af := alt_count / total_alleles]
dat[af > 0.5, maf := 1-af]
dat[af <= 0.5, maf := af]


HBCC_maf <- copy(dat$maf)

dat <- fread('NABEC_autosomes.traw')
gc()


sample_cols <- colnames(dat)[-c(1,2,3,4,5,6)]
n_samples <- length(sample_cols)
dat[, 'alt_count' := apply(.SD, 1, function(x) sum(x, na.rm=T)), .SDcols=sample_cols]
dat[, 'na_samples' := apply(.SD, 1, function(x) sum(is.na(x))), .SDcols=sample_cols]
dat[, total_alleles := 2*n_samples - 2*na_samples]
dat[, ref_count := total_alleles - alt_count]
dat[, af := alt_count / total_alleles]
dat[af > 0.5, maf := 1-af]
dat[af <= 0.5, maf := af]
NABEC_maf <- copy(dat$maf)


nabec <- data.table(table(cut(NABEC_maf, breaks=seq(0,1,0.01))))
setnames(nabec, c('bin','N'))
nabec[, 'bin' := gsub('\\(','',bin)]
nabec[, 'bin' := gsub('\\]','',bin)]
nabec[, 'bin' := tstrsplit(bin, split=',')[2]]
nabec[, cohort := 'NABEC']
nabec[, 'Alleles (Percent)' := N/sum(N)]


hbcc <- data.table(table(cut(HBCC_maf, breaks=seq(0,1,0.01))))
setnames(hbcc, c('bin','N'))
hbcc[, 'bin' := gsub('\\(','',bin)]
hbcc[, 'bin' := gsub('\\]','',bin)]
hbcc[, 'bin' := tstrsplit(bin, split=',')[2]]
hbcc[, cohort := 'HBCC']
hbcc[, 'Alleles (Percent)' := N/sum(N)]

dat <- rbindlist(list(hbcc, nabec))
dat[, bin := as.numeric(bin)]
dat.long <- melt(dat, measure.vars=c('N','Alleles (Percent)'))
dat.long[variable=='N', variable := 'Alleles (N)']

options(scipen=999)

g <- ggplot(dat.long[bin %between% c(0.045,0.5)], aes(x=bin, y=value, color=cohort)) +
    geom_line() +
    facet_grid(variable~., scales='free_y', switch='y') +
    scale_x_continuous(breaks=seq(0,0.5, 0.05), limits=c(0,0.5)) +
    labs(x='Minor allele frequency (1% bins)', y='') +
    theme_few() +
    theme(strip.text.y = element_text(angle = 180), 
      axis.title.y = element_text(vjust = -10),  # value by experiment
      strip.placement = "outside")

ggsave(g, file='PFC-HBCC-NABEC-AFS.png')
ggsave(g, file='PFC-HBCC-NABEC-AFS.svg')