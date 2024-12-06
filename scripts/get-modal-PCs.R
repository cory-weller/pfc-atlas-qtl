#!/usr/bin/env Rscript

# Script for generating principal components to serve as regression covariates during QTL analysis
# from expression (RNA) or chromatin accessibility (ATAC) data


library(data.table)


rename_header <- function(DT, cohort) {
    ids <- colnames(DT)[-1]
    newnames <- gsub('-ARC', '', ids)
    if(cohort == 'hbcc') {
        newnames <- gsub('^', 'HBCC_', newnames)
    } 
    return(c('feature', newnames))
}

get_pcs <- function(mode, celltype, cohort) {
    filename <- paste0('/data/CARD_singlecell/cortex_subtype/output/', mode, '/', celltype, '/pseudobulk/', cohort, '_cpm_sum_pseudobulk_model_counts.csv')
    dat <- fread(filename)
    setnames(dat, rename_header(dat, cohort))
    pcs <- prcomp(dat[, c(-1)])
    modal_pcs <- as.data.table(pcs$rotation)[,1:20]
    modal_pcs[, sample := rownames(pcs$rotation)]
    return(modal_pcs[, .SD, .SDcols=c('sample', paste0('PC',1:20))])
}


for(mode in c('atac','rna')) {
    for(celltype in c('Astro', 'ExN', 'InN', 'MG', 'Oligo', 'OPC', 'VC')) {
        for(cohort in c('hbcc','nabec')) {
            x <- get_pcs(mode, celltype, cohort)
            fwrite(x, file=paste0('data/', cohort, '-', celltype, '-', mode, '-PCs.tsv'), quote=F, row.names=F, col.names=T, sep='\t')
        }
    }
}




