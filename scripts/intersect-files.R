#!/usr/bin/env Rscript

library(data.table)

args <- commandArgs(trailingOnly=TRUE)
cohort <- args[1]
mode <- args[2]
celltype <- args[3]


counts_fn <- paste0('data/',mode,'/',celltype,'-counts.bed')
subset_counts_fn <- paste0('tensorqtl-subsets/', cohort, '-',mode,'-',celltype,'-counts.bed')
covariates_fn <- paste0('data/',cohort,'-covariates.txt')
subset_covariates_fn <- paste0('tensorqtl-subsets/', cohort, '-', mode, '-', celltype, '-covariates.txt')
interaction_fn <- paste0('data/',cohort,'-interaction.tsv')
subset_interaction_fn <- paste0('tensorqtl-subsets/', cohort, '-', mode, '-', celltype, '-interaction.tsv')
counts <- fread(counts_fn)
covariates <- fread(covariates_fn)
interaction <- fread(interaction_fn)

# Get counts samples
if(cohort == 'HBCC') {
    counts.samples <- grep('^HBCC', colnames(counts), value=TRUE)
} else if (cohort == 'NABEC') {
    counts.samples <- c(grep('^KEN-', colnames(counts), value=TRUE),
                    grep('^SH-', colnames(counts), value=TRUE),
                    grep('^UMARY-', colnames(counts), value=TRUE))
}

# Get covariates samples
covariates.samples <- colnames(covariates)[-1]

# Get interaction samples
interaction.samples <- interaction$FID

# Get intersection
sampleset <- intersect(interaction.samples, intersect(counts.samples, covariates.samples))

# Write intersection
samples.DT <- data.table(FID=sampleset, IID=sampleset)

fwrite(samples.DT, file=paste0('tensorqtl-subsets/', cohort, '-', mode, '-', celltype, '.txt'),
row.names=F, col.names=T, sep='\t', quote=F)

# Subset counts
counts_subset <- counts[, .SD, .SDcols=c('#chr','start','end','phenotype_id',sampleset)]
fwrite(counts_subset, file=subset_counts_fn, quote=F, row.names=F, col.names=T, sep='\t')

# Subste covariates
covariates <- covariates[, .SD, .SDcols=c('V1', sampleset)]
covariates_subset <- as.data.frame(covariates[, .SD, .SDcols=sampleset])
rownames(covariates_subset) <- covariates$V1
fwrite(covariates_subset, file=subset_covariates_fn, quote=F, sep='\t', col.names=T, row.names=T)

# Subset interaction
setkey(interaction, FID)
interaction <- interaction[sampleset]
fwrite(interaction[FID %in% sampleset], file=subset_interaction_fn, quote=F, row.names=F, col.names=T, sep='\t')


stopifnot(identical(interaction$FID, colnames(covariates)[-1]))
stopifnot(identical(colnames(covariates)[-1], colnames(counts_subset)[5:ncol(counts_subset)]))