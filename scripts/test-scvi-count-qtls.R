#!/usr/bin/env Rscript

library(arrow)
library(data.table)
library(foreach)
library(doMC)
registerDoMC(cores=4)

cohorts <- c('HBCC','NABEC')
celltypes <- c('ExN','InN')
chrs <- paste0('chr',1:22)
modes <- c('rna')
combos <- CJ('cohort'=cohorts, 'celltype'=celltypes, 'mode'=modes, 'chr'=chrs)

#summary_stats_fn <- 'cis-QTLs-bonferroni.tsv'

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

# For each row of `parquet_counts`, import the parquet file
# Return only rows with significant Genetic 
import_parquet <- function(cohort, celltype, mode, chr) {
    fn <- paste0('tensorqtl_sumScvi_cpm/', cohort, '_', celltype, '_', mode, '.cis_qtl_pairs.', chr, '.parquet')
    dat <- arrow::read_parquet(fn)
    setDT(dat)
    dat[, GenoAssoc := FALSE]
    dat[, Interaction := FALSE]
    dat[, AgeAssoc := FALSE]
    dat[, 'celltype' := celltype]
    dat[, 'cohort' := cohort]
    dat[, 'mode' := mode]
    dat[, 'chr' := chr]
    return(dat[])
}



if(!file.exists(summary_stats_fn)) {
   o <- foreach(i=1:nrow(combos)) %dopar% {
        row_i <- combos[i,]
        cohort_i <- row_i[['cohort']]
        celltype_i <- row_i[['celltype']]
        mode_i <- row_i[['mode']]
        chr_i <- row_i[['chr']]

        import_parquet(cohort = cohort_i,
                        celltype = celltype_i,
                        mode = mode_i,
                        chr = chr_i
        )

    }   
    o <- rbindlist(o)
    # merge in gene IDs
    fwrite(o, file=summary_stats_fn, sep='\t')
    rm(o)
    gc()
}



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


