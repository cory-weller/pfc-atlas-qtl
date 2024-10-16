#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(ggthemes)
library(viridis)

# Import covs
covariates <- fread('samples.tsv')
covariates[Brain_bank == 'HBCC', cohort := 'HBCC']
covariates[Brain_bank != 'HBCC', cohort := 'NABEC']
covariates[cohort == 'HBCC', sample := paste0('HBCC_',sample)]
covariates[, 'FID' := gsub('-ARC', '', sample)]

# Import PCs from king
nabec.geno.pcs <- fread('genotypes/NABEC_polarized.pruned_pc.txt')
hbcc.geno.pcs <- fread('genotypes/HBCC_polarized.pruned_pc.txt')

# Combine dtables
hbcc <- merge(hbcc.geno.pcs, covariates, by.x='FID', by.y='FID')
nabec <- merge(nabec.geno.pcs, covariates, by.x='FID', by.y='FID')
dat <- rbindlist(list(hbcc, nabec))

# Convert to factor when necessary
dat[, Sex := factor(Sex)]
dat[, Homogenization := factor(Homogenization)]
dat[, Brain_bank := factor(Brain_bank)]
dat[, Ancestry := factor(Ancestry)]
dat[, LibraryPrep := factor(LibraryPrep)]
dat[, Sequencing := factor(Sequencing)]

# Save as tabular file
fwrite(dat, file='data/all-covariates.tsv', quote=F, row.names=F, col.names=T, sep='\t')

# Remove samples determined to be PC outliers for HBCC
exclude_hbcc <- fread('genotypes/HBCC_to_remove.txt')$FID
dat <- dat[! FID %in% exclude_hbcc]
dat[, Sex := toupper(Sex)]

# Convert sex to dichotomous factor
dat[Sex=='MALE', Sex := 0]
dat[Sex=='FEMALE', Sex := 1]

# Generate covariate files
desired_covs <- c(paste0('PC',1:10), 'Sex', 'Age', 'PMI')
hbcc <- dat[cohort == 'HBCC'][, .SD, .SDcols=c('FID',desired_covs)]
nabec <- dat[cohort == 'NABEC'][, .SD, .SDcols=c('FID',desired_covs)]
fwrite(hbcc, file='data/HBCC-covariates.tsv', quote=F, row.names=F, col.names=T, sep='\t')
fwrite(nabec, file='data/NABEC-covariates.tsv', quote=F, row.names=F, col.names=T, sep='\t')

# # Save interaction (age)
## DEPRECATED: Instead covariates are pulled out of the covariates file, specified by --covariates
## when running intersect-files.R

# fwrite(hbcc[, .SD, .SDcols=c('FID','Age')], file='data/HBCC-interaction.tsv', quote=F, row.names=F, col.names=T, sep='\t')
# fwrite(nabec[, .SD, .SDcols=c('FID','Age')], file='data/NABEC-interaction.tsv', quote=F, row.names=F, col.names=T, sep='\t')

# hbcc.t <- t(hbcc[, !c('Age','FID')])
# colnames(hbcc.t) <- hbcc$FID
# write.csv(hbcc.t, file='data/HBCC-covariates.txt', quote=F)


# nabec <- dat[cohort == 'NABEC'][, .SD, .SDcols=c('FID', paste0('PC',1:10), 'Sex', 'Age', 'PMI')]
# # Save interaction (age)
# 
# nabec.t <- t(nabec[, !c('Age','FID')])
# colnames(nabec.t) <- nabec$FID
# write.csv(nabec.t, file='data/NABEC-covariates.txt', quote=F)