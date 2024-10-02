#!/usr/bin/env Rscript

library(arrow)
library(data.table)
library(foreach)
library(doMC)
registerDoMC(cores=4)


args <- commandArgs(trailingOnly=TRUE)
tensorqtl_output_dir <- 'QTL-scvicounts-polarized-interaction'
celltypes <- c('Astro','ExN','InN','MG','Oligo','OPC','VC')
filenames <- list.files(tensorqtl_output_dir, pattern='*.parquet', full.names=T, recursive=T)

# # fn is shorthand for filename
# dat <- foreach(fn=filenames, .combine='rbind', .errorhandling='remove') %do% {
#     dat.tmp <- fread(fn)
#     fn_split <- unlist(strsplit(basename(fn), split='_|\\.'))
#     cohort <- fn_split[1]
#     celltype <- fn_split[2]
#     mode <- fn_split[3]
#     return(cbind(dat.tmp, cohort, celltype, mode))
# }

#  fwrite(dat, file=paste0(output_dir, '/all_top_significant_hits.tsv'), quote=F, row.names=F, col.names=T, sep='\t')

summary_stats_fn <- paste0(tensorqtl_output_dir, '/cis-QTLs-Bonferroni.tsv')

get_bonferroni <- function(N_tests, n_celltypes, n_terms, alpha=0.5) {
    # Returns bonferroni-adjusted p-value threshold at a given alpha
    threshold <- alpha / ( N_tests / n_celltypes / n_terms)
    return(threshold)
}

read_parquet <- function(fn) {
    # Accepts a parquet file name (from tensorQTL) and outputs a data.table with the cols
    # cohort    |   ('HBCC','NABEC')
    # celltype  |   ('ExN','InN','Astro','MG','Oligo','OPC','VC')
    # mode      |   ('atac','rna'
    # chr       |   ('chr1' - 'chr22') i.e. autosomes only
    # N         |   integer, N rows (tests done)
    dat <- arrow::read_parquet(fn)
    setDT(dat)
    N <- nrow(dat)
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

check_ntests <- function(files) {
    o <- foreach(file=files, .combine='rbind') %dopar% {
        count_parquet(file)
    }
    return(o)
}

# Iterate over parquet files in the 'output' directory
parquet_counts <- check_ntests(files=filenames)
fwrite(parquet_counts, file=paste0(tensorqtl_output_dir, '/parquet_counts.tsv'), quote=F, row.names=F, col.names=T, sep='\t')

# Get bonferroni-adjusted p value threshold for ATAC QTLs
hbcc_threshold <-   get_bonferroni(N_tests=sum(parquet_counts[cohort=='HBCC', N]),
                        n_celltypes=length(celltypes),
                        n_terms=2,
                        alpha=0.1
                    )

# Get bonferroni-adjusted p value threshold for RNA QTLs
nabec_threshold <-    get_bonferroni(N_tests=sum(parquet_counts[cohort=='NABEC', N]), 
                            n_celltypes=length(celltypes), 
                            n_terms=2, 
                            alpha=0.1
                    )

rna_threshold <-   get_bonferroni(N_tests=sum(parquet_counts[mode=='rna', N]),
                        n_celltypes=length(celltypes),
                        n_terms=2,
                        alpha=0.1
                    )

# Get bonferroni-adjusted p value threshold for RNA QTLs
atac_threshold <-    get_bonferroni(N_tests=sum(parquet_counts[mode=='atac', N]), 
                            n_celltypes=length(celltypes), 
                            n_terms=2, 
                            alpha=0.1
                    )

# For each row of `parquet_counts`, import the parquet file
# Return only rows with significant Genetic 
import_parquet <- function(cohort, celltype, mode, chr, threshold) {
    fn <- paste0('output/', cohort, '_', celltype, '_', mode, '.cis_qtl_pairs.', chr, '.parquet')
    dat <- arrow::read_parquet(fn)
    setDT(dat)

}

import_parquet <- function(fn, atac_thresh, rna_thresh) {
    dat <- arrow::read_parquet(fn)
    setDT(dat)
    s <- strsplit(basename(fn), split='_|\\.')
    s <- unlist(s)
    cohort <- s[1]
    celltype <- s[2]
    mode <- s[3]
    chr <- s[7]
    if(mode == 'rna') {
        threshold <- rna_thresh
    } else if(mode == 'atac') {
        threshold <- atac_thresh
    }

    dat[, GenoAssoc := FALSE]
    dat[, Interaction := FALSE]
    dat[, AgeAssoc := FALSE]
    dat[pval_gi < threshold, Interaction := TRUE]
    dat[pval_g < threshold, GenoAssoc := TRUE]
    dat[, 'celltype' := celltype]
    dat[, 'cohort' := cohort]
    dat[, 'mode' := mode]
    dat[, 'chr' := chr]
    return(dat[Interaction == TRUE][])
}




if(!file.exists(summary_stats_fn)) {
   o <- foreach(file=filenames) %dopar% {
       
        import_parquet(fn=file,
                        atac_thresh=atac_threshold,
                        rna_thresh=rna_threshold)


    }   
    o <- rbindlist(o)
    # merge in gene IDs
    fwrite(o, file=summary_stats_fn, sep='\t')
    rm(o)
    gc()
}


dat <- fread(summary_stats_fn)


# Get variants for which associations are found
if(!dir.exists('qtl_analysis')) {dir.create('qtl_analysis')}
cat(unique(sort(dat$variant_id)), file='qtl_analysis/variant_ids.txt', sep='\n')
