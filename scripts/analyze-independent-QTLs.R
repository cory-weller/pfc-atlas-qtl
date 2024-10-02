#!/usr/bin/env Rscript

#library(arrow)
library(data.table)
library(foreach)
library(doMC)
registerDoMC(cores=4)

output_dir <- 'tensorqtl_permutations'

filenames <- list.files(output_dir, pattern='*independent_qtl.txt.gz', full.names=T)

# fn is shorthand for filename
dat <- foreach(fn=filenames, .combine='rbind', .errorhandling='remove') %do% {
    dat.tmp <- fread(fn)
    fn_split <- unlist(strsplit(basename(fn), split='_|\\.'))
    cohort <- fn_split[1]
    celltype <- fn_split[2]
    mode <- fn_split[3]
    return(cbind(dat.tmp, cohort, celltype, mode))
}

fwrite(dat, file=paste0(output_dir, '/all_top_significant_hits.tsv'), quote=F, row.names=F, col.names=T, sep='\t')

out_dir <- 'tensorqtl_sumScvi_cpm'
summary_stats_fn <- paste0(out_dir, '/cis-QTLs-bonferroni.tsv')
if(! dir.exists(out_dir) ) { dir.create(out_dir) }

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

# Iterate over parquet files in the 'output' directory
parquet_counts <- check_ntests(dir='output')
fwrite(parquet_counts, file=paste0(out_dir, '/parquet_counts.tsv'), quote=F, row.names=F, col.names=T, sep='\t')

# Get bonferroni-adjusted p value threshold for ATAC QTLs
hbcc_threshold <-   get_bonferroni(N_tests=sum(parquet_counts[cohort=='HBCC', N]),
                        n_celltypes=length(celltypes),
                        n_terms=2,
                        alpha=0.05
                    )

# Get bonferroni-adjusted p value threshold for RNA QTLs
nabec_threshold <-    get_bonferroni(N_tests=sum(parquet_counts[cohort=='NABEC', N]), 
                            n_celltypes=length(celltypes), 
                            n_terms=2, 
                            alpha=0.05
                    )

# For each row of `parquet_counts`, import the parquet file
# Return only rows with significant Genetic 
import_parquet <- function(cohort, celltype, mode, chr, threshold) {
    fn <- paste0('output/', cohort, '_', celltype, '_', mode, '.cis_qtl_pairs.', chr, '.parquet')
    dat <- arrow::read_parquet(fn)
    setDT(dat)
    dat[, GenoAssoc := FALSE]
    dat[, Interaction := FALSE]
    dat[, AgeAssoc := FALSE]
    dat[pval_gi < threshold, Interaction := TRUE]
    dat[pval_g < threshold, GenoAssoc := TRUE]
    dat[, 'celltype' := celltype]
    dat[, 'cohort' := cohort]
    dat[, 'mode' := mode]
    dat[, 'chr' := chr]
    return(dat[GenoAssoc==TRUE | Interaction == TRUE][])
}




if(!file.exists(summary_stats_fn)) {
   o <- foreach(i=1:nrow(parquet_counts)) %dopar% {
        row_i <- parquet_counts[i,]
        cohort_i <- row_i[['cohort']]
        celltype_i <- row_i[['celltype']]
        mode_i <- row_i[['mode']]
        chr_i <- row_i[['chr']]
        if(cohort_i == 'HBCC') {
            p_threshold_i <- hbcc_threshold
        } else if(cohort_i == 'NABEC') {
            p_threshold_i <- nabec_threshold
        }
        import_parquet(cohort = cohort_i,
                        celltype = celltype_i,
                        mode = mode_i,
                        chr = chr_i,
                        threshold = p_threshold_i
        )

    }   
    o <- rbindlist(o)
    # merge in gene IDs
    fwrite(o, file=summary_stats_fn, sep='\t')
    rm(o)
    gc()
}

dat <- fread(summary_stats_fn)


geno_assocs <- dat[GenoAssoc==TRUE, .N, by=list(celltype, cohort, mode)][order(celltype, cohort, mode)]
interact_assoc <- dat[Interaction==TRUE, .N, by=list(celltype, cohort, mode)][order(celltype, cohort, mode)]
all_assoc <- dat[Interaction==TRUE & GenoAssoc == TRUE]

# Get variants for which associations are found
if(!dir.exists('qtl_analysis')) {dir.create('qtl_analysis')}
cat(unique(sort(dat$variant_id)), file='qtl_analysis/variant_ids.txt', sep='\n')


