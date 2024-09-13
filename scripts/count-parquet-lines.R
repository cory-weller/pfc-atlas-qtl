#!/usr/bin/env Rscript

library(arrow)
library(data.table)
library(foreach)
library(doMC)
registerDoMC(cores=4)

cohorts <- c('HBCC','NABEC')
celltypes <- c('ExN','InN','Astro','MG','Oligo','OPC','VC')
chrs <- paste0('chr',1:22)
modes <- c('rna')


import_parquet <- function(cohort, celltype, mode, chr) {
    fn <- paste0('output/', cohort, '_', celltype, '_', mode, '.cis_qtl_pairs.', chr, '.parquet')
    dat <- arrow::read_parquet(fn)
    setDT(dat)
    N <- nrow(dat)
    data.table(cohort, celltype, mode, chr, N)
}

combos <- CJ('cohort'=cohorts, 'celltype'=celltypes, 'mode'=modes, 'chr'=chrs)

dat <- foreach('cohort'=combos$cohort, 'celltype'=combos$celltype, 'mode'=combos$mode, 'chr'=combos$chr, .combine='rbind') %dopar% {
    dt <- import_parquet(cohort, celltype, mode, chr)
    return(dt)
}



