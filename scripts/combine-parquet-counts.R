#!/usr/bin/env Rscript

library(arrow)
library(data.table)
library(foreach)
library(doMC)
registerDoMC(cores=12)

# cohorts <- c('HBCC','NABEC')
# celltypes <- c('ExN','InN')
# chrs <- paste0('chr',1:22)
# modes <- c('rna', 'atac')
# combos <- CJ('cohort'=cohorts, 'celltype'=celltypes, 'mode'=modes, 'chr'=chrs)

output_dir <- 'QTL'
summary_stats_fn <- 'QTL/cis-QTLs-bonferroni.tsv'

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

files <- data.table(fn = list.files('QTL', pattern='*.parquet', recursive=T, full.names=T))
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



for (i in grp_ids) {
    print(i)
    o <- foreach(filename=files[grp==i, fn], .combine='rbind', .errorhandling='remove') %do% {
        import_parquet(filename)
    }
    o[, 'pval_bh' := p.adjust(pval_nominal, method='BH')]
    cohort <- files[grp==i][1,cohort]
    celltype <- files[grp==i][1,celltype]
    mode <- files[grp==i][1,mode]
    
    out_fn <- paste0('QTL/', cohort, '-', celltype, '-', mode, '-ALL.csv')
    fwrite(o, file=out_fn, sep=',', quote=F, row.names=F, col.names=T)
    
    o <- o[pval_bh < 0.05]
    out_fn <- paste0('QTL/', cohort, '-', celltype, '-', mode, '-BH05.csv')
    fwrite(o, file=out_fn, sep=',', quote=F, row.names=F, col.names=T)
    rm(o)
    gc()
}



# Get N counts
files <- list.files('QTL', pattern='*-ALL.csv', full.names=TRUE, recursive=TRUE)

ntests <- foreach(file=files, .combine='rbind') %do% {
    fread(cmd=paste0('wc -l ', file))
}

setnames(ntests, c('N','fn'))
ntests[, 'tmp' := gsub('QTL/', '', fn)]
ntests[, 'tmp' := gsub('-ALL.csv', '', tmp)]
ntests[, c('cohort','celltype','mode') := tstrsplit(tmp, split='-')]
ntests[, 'tmp' := NULL]
ntests[, 'fn' := NULL]
ntests[, sum(N), by=cohort]






dat <- fread('QTL/HBCC-Astro-rna-ALL.csv', select='pval_nominal')
dat2 <- fread('QTL/HBCC-ExN-rna-ALL.csv', select='pval_nominal')
dat3<- fread('QTL/HBCC-InN-rna-ALL.csv', select='pval_nominal')
dat4 <- fread('QTL/HBCC-MG-rna-ALL.csv', select='pval_nominal')
dat5 <- fread('QTL/HBCC-Oligo-rna-ALL.csv', select='pval_nominal')
dat6 <- fread('QTL/HBCC-OPC-rna-ALL.csv', select='pval_nominal')
dat7 <- fread('QTL/HBCC-VC-rna-ALL.csv', select='pval_nominal')





# Iterate over parquet files in the 'output' directory
parquet_counts <- check_ntests(dir='tensorqtl_sumScvi_cpm')

# Get bonferroni-adjusted p value threshold for ATAC QTLs
atac_threshold <-   get_bonferroni(N_tests=sum(parquet_counts[mode=='atac', N]),
                        n_celltypes=length(celltypes),
                        n_terms=2,
                        alpha=0.05
                    )

# Get bonferroni-adjusted p value threshold for RNA QTLs
rna_threshold <-    get_bonferroni(N_tests=sum(parquet_counts[mode=='rna', N]), 
                            n_celltypes=length(celltypes), 
                            n_terms=2, 
                            alpha=0.05
                    )



fwrite(o[cohort=='HBCC'], file='hbcc_scvi_sum.tsv', quote=F, row.names=F, col.names=T, sep='\t')
fwrite(o[cohort=='NABEC'], file='nabec_scvi_sum.tsv', quote=F, row.names=F, col.names=T, sep='\t')

nabec <- o[cohort=='NABEC']
hbcc <- o[cohort=='HBCC']
nabec[, 'p.adjust' := p.adjust(pval_g, method='BH')]
hbcc[, 'p.adjust' := p.adjust(pval_g, method='BH')]

dat <- rbindlist(list(
nabec[p.adjust < 0.05],
hbcc[p.adjust < 0.05]
))

dat <- fread(summary_stats_fn)


geno_assocs <- dat[GenoAssoc==TRUE, .N, by=list(celltype, cohort, mode)][order(celltype, cohort, mode)]
interact_assoc <- dat[Interaction==TRUE, .N, by=list(celltype, cohort, mode)][order(celltype, cohort, mode)]
all_assoc <- dat[Interaction==TRUE & GenoAssoc == TRUE]

# Get variants for which associations are found
if(!dir.exists('qtl_analysis')) {dir.create('qtl_analysis')}
cat(unique(sort(dat$variant_id)), file='qtl_analysis/variant_ids.txt', sep='\n')


