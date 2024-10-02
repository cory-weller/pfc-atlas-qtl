#!/usr/bin/env Rscript


silent <- function(pkgs) {
    for(pkg in pkgs) {
        suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    }
}

silent(c('argparse','data.table','GenomicRanges'))

# args <- commandArgs(trailingOnly=TRUE)
debug_args <- c(  '--cohort', 'HBCC',
            '--mode', 'atac',
            '--celltype', 'Astro',
            '--projdir', '/home/wellerca/pfc-atlas-qtl',
            '--exclude', 'data/TSS-blacklist.bed' 
)


parser <- ArgumentParser()

parser$add_argument("--cohort", type='character', nargs=1, required=TRUE, dest="cohort", 
                    choices=c('HBCC','NABEC')
                    )

parser$add_argument("--mode", type='character', nargs=1, required=TRUE, dest="mode", 
                    choices=c('rna','atac')
                    )

parser$add_argument("--celltype", type='character', nargs=1, required=TRUE, dest="celltype", 
                    choices=c('Oligo', 'Astro', 'ExN', 'InN', 'MG', 'OPC', 'VC')
                    )

parser$add_argument("--projdir", type='character', nargs=1, required=TRUE, dest="projdir", 
                    help="Top level directory of project"
                    )

parser$add_argument("--exclude", type='character', nargs=1, default=NULL, required=FALSE, 
                    dest="exclude", help="bed file containing feature blacklists with cols chr, start, stop"
                    )


if(length(commandArgs(trailingOnly=TRUE)) == 0) {
    args <- parser$parse_args(args=debug_args)
} else {
    args <- parser$parse_args()
}

# For compatibility, use single-word variable names from args
cohort <- args$cohort
mode <- args$mode
celltype <- args$celltype
projdir <- args$projdir
exclude <- args$exclude



counts_fn <- paste0(projdir, '/QTL-pseudobulk-counts/',mode,'-', cohort, '-', celltype,'-counts.bed')
subset_counts_fn <- 'pseudobulk-counts.bed'
covariates_fn <- paste0(projdir, '/data/',cohort,'-covariates.tsv')
subset_covariates_fn <- paste0('covariates.txt')
interaction_fn <- paste0(projdir, '/data/',cohort,'-interaction.tsv')
subset_interaction_fn <- paste0('interaction.tsv')

counts <- fread(counts_fn)
covariates <- fread(covariates_fn)
interaction <- fread(interaction_fn)

# Get feature blacklist
## prepend project dir if relative path given
if(! is.null(exclude)) {
    if(! startsWith(exclude, '/')) {
        exclude <- paste0(projdir, '/', exclude)
    }
    # Import blacklist bed file and convert to GRanges
    exclude <- fread(exclude)
    exclude <- GRanges(exclude)

    # Collect pseudobulk features from counts table and convert to GRanges
    # If ATAC data, get range from phenotype_id
    if (mode == 'atac') {
        features <- copy(counts[, .SD, .SDcols=c('phenotype_id')])
        features[, c('chr','start','stop') := tstrsplit(phenotype_id, split=':|-')]
        setcolorder(features, c('chr','start','stop','phenotype_id'))
    } else if(mode=='rna') {
        # When RNA
        features <- copy(counts[, .SD, .SDcols=c('#chr','start','end','phenotype_id')])
        setnames(features, '#chr','chr')
        setnames(features, 'end','stop')
    }
    features <- GRanges(features)

    # Get overlapping regions between counts features and blacklist
    hits <- GenomicRanges::findOverlaps(features, exclude)
    overlaps <- pintersect(Pairs(features, exclude, hits=hits))
    overlaps <- as.data.table(overlaps)

    to_exclude <- overlaps[width > 9, phenotype_id]
    counts <- counts[! phenotype_id %in% to_exclude]
}


# Get counts samples
if(cohort == 'HBCC') {
    counts.samples <- grep('^HBCC', colnames(counts), value=TRUE)
} else if (cohort == 'NABEC') {
    counts.samples <- c(grep('^KEN-', colnames(counts), value=TRUE),
                    grep('^SH-', colnames(counts), value=TRUE),
                    grep('^UMARY-', colnames(counts), value=TRUE))
}

# Get covariates samples
covariates.samples <- covariates$FID

# Get interaction samples
interaction.samples <- interaction$FID

# Get intersection
sampleset <- intersect(counts.samples, covariates.samples)
sampleset <- intersect(sampleset, interaction.samples)

# Write intersection
samples.DT <- data.table(FID=sampleset, IID=sampleset)
setkey(samples.DT)
samples.DT <- samples.DT[sampleset]

setkey(covariates, FID)
covariates <- covariates[sampleset]

setkey(interaction, FID)
interaction <- interaction[sampleset]


fwrite(samples.DT, file='samples.txt',row.names=F, col.names=T, sep='\t', quote=F)

# Subset counts
counts_subset <- counts[, .SD, .SDcols=c('#chr','start','end','phenotype_id',sampleset)]
fwrite(counts_subset, file=subset_counts_fn, quote=F, row.names=F, col.names=T, sep='\t')

# Subset and order covariate table sampleIDs the same as all other files
setkey(covariates, FID)
#covariates[, FID := NULL]
covariates_cols <- colnames(covariates)
covariates <- t(as.data.frame(covariates))
rownames(covariates) <- c('', rownames(covariates)[2:nrow(covariates)])
write.table(covariates, file='covariates.txt', quote=F, col.names=F, sep='\t')


# Subset interaction
fwrite(interaction, file=subset_interaction_fn, quote=F, row.names=F, col.names=T, sep='\t')


stopifnot(identical(covariates[1,], colnames(counts_subset)[5:ncol(counts_subset)]))