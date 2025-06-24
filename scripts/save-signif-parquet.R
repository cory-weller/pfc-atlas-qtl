#!/usr/bin/env Rscript

library(arrow)
library(data.table)
library(foreach)
library(doMC)
registerDoMC(cores=12)

args <- commandArgs(trailingOnly=TRUE)

dir_name <- args[1]
batch <- basename(dir_name)

output_dir <- 'QTL'

get_bonferroni <- function(N_tests, n_celltypes, n_terms, alpha=0.05) {
    # Returns bonferroni-adjusted p-value threshold at a given alpha
    threshold <- alpha / ( N_tests / n_celltypes / n_terms)
    return(threshold)
}

count_parquet <- function(fn) {
    # Accepts a parquet file name (from tensorQTL) and outputs a data.table with the cols
    # cohort    |   ('HBCC','NABEC')
    # celltype  |   ('ExN','InN','Astro','MG','Oligo','OPC','VC')
    # mode      |   ('atac','rna'
    # chr       |   ('chr1' - 'chr22') i.e. autosomes only
    # N         |   integer, N rows (tests done)
    N <- nrow(arrow::read_parquet(fn))
    s <- strsplit(basename(fn), split='_|\\.')
    s <- unlist(s)
    cohort <- s[1]
    celltype <- s[2]
    mode <- s[3]
    chr <- s[7]
    if(chr == 'chrX') {
        return(NULL)
    }
    out <- data.table(cohort, celltype, mode, chr, N)
    return(out)
}

check_ntests <- function(dir) {
    files <- list.files(dir, pattern='*.parquet', full.names=T)
    o <- foreach(file=files, .combine='rbind') %dopar% {
        count_parquet(file)
    }
    return(o)
}

files <- data.table(fn = list.files(dir_name, pattern='*.parquet', recursive=T, full.names=T))
files[, idx := 1:.N]
files[, c('cohort','celltype','mode'):= tstrsplit(basename(dirname(fn)), split='-'), by=idx]
files[, chr := gsub('cis-nominal.cis_qtl_pairs.', '', gsub('.parquet', '', basename(fn)))]
setkey(files, cohort, celltype, mode, chr)
files[, grp := rleid(cohort, celltype, mode)]

grp_ids <- unique(files$grp)


# For each row of `parquet_counts`, import the parquet file
# Return only rows with significant Genetic 
import_parquet <- function(fn) {
    
    splitstring <- unlist(strsplit(basename(dirname(fn)), split='-'))
    cohort <- splitstring[1]
    celltype <- splitstring[2]
    mode <- splitstring[3]
    chr <- gsub('.parquet', '', gsub('cis-nominal.cis_qtl_pairs.','', basename(fn)))
    
    dat <- arrow::read_parquet(fn)
    setDT(dat)
    # dat[, GenoAssoc := FALSE]
    # dat[, Interaction := FALSE]
    # dat[, AgeAssoc := FALSE]
    dat[, 'celltype' := celltype]
    dat[, 'cohort' := cohort]
    dat[, 'mode' := mode]
    dat[, 'chr' := chr]
    dat <- dat[!is.na(af)]
    #dat <- dat[! is.na(pval_bh)][]
    return(dat[])
}

if(batch %like% 'rna') {
    o2 <- foreach(filename=files$fn, .combine='rbind', .errorhandling='remove') %do% {
        o <- import_parquet(filename)
        o[, pval_bonferroni := p.adjust(pval_nominal, method='bonferroni')]
        o[, pval_BH := p.adjust(pval_nominal, method='BH')]
        o[pval_bonferroni <= 0.05][]
    }
    fwrite(o2[pval_bonferroni < 0.05], file=paste0('QTL/', batch, '-signif-bonferroni.csv'), sep=',', quote=F, row.names=F, col.names=T)
    fwrite(o2[pval_BH < 0.05], file=paste0('QTL/', batch, '-signif-BH.csv'), sep=',', quote=F, row.names=F, col.names=T)
} else if(batch %like% 'atac') {
    for(filename in files$fn) {
        o <- import_parquet(filename)
        out_fn <- gsub('cis-nominal.cis_qtl_pairs.', '', filename)
        out_fn <- gsub('.parquet', '.csv', out_fn)         
        fwrite(o, file=out_fn, quote=F, row.names=F, col.names=T, sep=',')
    }
}

# o <- o[, .SD, .SDcols=c('pval_nominal')][order(pval_nominal)]
# out_fn <- paste0('QTL/', cohort, '-', celltype, '-', mode, '-nominal_p.txt')
# fwrite(o, file=out_fn, sep=',', quote=F, row.names=F, col.names=F)



quit()
