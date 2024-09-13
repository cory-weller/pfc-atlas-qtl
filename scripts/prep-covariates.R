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
nabec.geno.pcs <- fread('data/genotypes/NABEC.pruned_pc.txt')
hbcc.geno.pcs <- fread('data/genotypes/HBCC.pruned_pc.txt')

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

# Define plotting function
save_feature_plot <- function(feature, DT, .cohort) {
    DT <- DT[cohort == .cohort]
    fn <- paste0('QC-plots/', cohort, '-', feature, '.png')
    g <- ggplot(DT, aes(x=PC1, y=PC2, color=get(feature))) + geom_point(shape=21) +
    labs(color=feature, title=cohort) +
    theme_few()

    if(! is.factor(DT[[feature]])) {
        g <- g + scale_color_viridis(limits=c(20,100))
    }

    ggsave(g, file=fn, width=15, height=12, units='cm')
}

# Get feature plots
features <- c('Age','Ancestry','Sex','Homogenization','LibraryPrep','Sequencing')
cohorts <- c('NABEC','HBCC')

for(cohort in cohorts) {
    for(feature in features) {
        save_feature_plot(feature, DT=dat, .cohort=cohort)
    }
}




exclude_hbcc <- fread('data/HBCC-remove.txt')$FID
dat <- dat[! FID %in% exclude_hbcc]
dat[, Sex := toupper(Sex)]

# Convert sex to dichotomous factor
dat[Sex=='MALE', Sex := 0]
dat[Sex=='FEMALE', Sex := 1]


hbcc <- dat[cohort == 'HBCC'][, .SD, .SDcols=c('FID', paste0('PC',1:10), 'Sex', 'Age', 'PMI')]
# Save interaction (age)
fwrite(hbcc[, .SD, .SDcols=c('FID','Age')], file='data/HBCC-interaction.tsv', quote=F, row.names=F, col.names=T, sep='\t')
hbcc.t <- t(hbcc[, !c('Age','FID')])
colnames(hbcc.t) <- hbcc$FID
write.csv(hbcc.t, file='data/HBCC-covariates.txt', quote=F)


nabec <- dat[cohort == 'NABEC'][, .SD, .SDcols=c('FID', paste0('PC',1:10), 'Sex', 'Age', 'PMI')]
# Save interaction (age)
fwrite(nabec[, .SD, .SDcols=c('FID','Age')], file='data/NABEC-interaction.tsv', quote=F, row.names=F, col.names=T, sep='\t')
nabec.t <- t(nabec[, !c('Age','FID')])
colnames(nabec.t) <- nabec$FID
write.csv(nabec.t, file='data/NABEC-covariates.txt', quote=F)