#!/usr/bin/env Rscript

library(arrow)
library(data.table)
library(foreach)
library(doMC)
registerDoMC(cores=4)

library(argparse)

# args <- commandArgs(trailingOnly=TRUE)
debug_args <- c(
            '--batch', 'QTL-output/scvicounts-interaction-oct15',
            '--method', 'interaction'
)

parser <- ArgumentParser()


# parser$add_argument("--cohort", type='character', nargs=1, required=TRUE, dest="cohort", 
#                     choices=c('HBCC','NABEC')
#                     )

# parser$add_argument("--mode", type='character', nargs=1, required=TRUE, dest="mode", 
#                     choices=c('rna','atac')
#                     )

# parser$add_argument("--celltype", type='character', nargs=1, required=TRUE, dest="celltype", 
#                     choices=c('Oligo', 'Astro', 'ExN', 'InN', 'MG', 'OPC', 'VC')
#                     )


parser$add_argument("--batchdir", type='character', nargs=1, required=TRUE, dest="batch", 
                    help="Directory containing tensorQTL parquet files"
                    )
parser$add_argument("--method", type='character', nargs=1, required=TRUE, dest="method", 
                    help="'nominal' or 'interaction'"
                    )

if(length(commandArgs(trailingOnly=TRUE)) == 0) {
    args <- parser$parse_args(args=debug_args)
} else {
    args <- parser$parse_args()
}

batch <- args$batch

# cohort <- args$cohort
# mode <- args$mode
# celltype <- args$celltype


cohorts <- c('HBCC','NABEC')
celltypes <- c('ExN','InN','Astro','MG','Oligo','OPC','VC')
chrs <- paste0('chr',1:22)
modes <- c('rna','atac')
#combos <- CJ('cohort'=cohorts, 'celltype'=celltypes, 'mode'=modes, 'chr'=chrs)

summary_stats_fn <- paste0(batch, '/cis-QTLs-bonferroni.tsv')

# if(! dir.exists(batch_dir) ) { dir.create(batch_dir) }

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

check_ntests <- function(filenames) {
    o <- foreach(file=filenames, .combine='rbind') %dopar% {
        count_parquet(file)
    }
    return(o)
}

filenames <- list.files(batch, pattern='*.parquet', full.names=T, recursive=T)


# Iterate over parquet files in the 'output' directory
parquet_counts <- check_ntests(filenames)

fwrite(parquet_counts, file=paste0(batch, '/parquet_counts.tsv'), quote=F, row.names=F, col.names=T, sep='\t')

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
NABEC_atac_threshold <-    get_bonferroni(N_tests=sum(parquet_counts[mode=='atac' & cohort=='NABEC', N]), 
                            n_celltypes=length(celltypes), 
                            n_terms=4, 
                            alpha=0.05
                    )
NABEC_rna_threshold <-    get_bonferroni(N_tests=sum(parquet_counts[mode=='rna' & cohort=='NABEC', N]), 
                            n_celltypes=length(celltypes), 
                            n_terms=4, 
                            alpha=0.05
                    )
HBCC_atac_threshold <-    get_bonferroni(N_tests=sum(parquet_counts[mode=='atac' & cohort=='HBCC', N]), 
                            n_celltypes=length(celltypes), 
                            n_terms=4, 
                            alpha=0.05
                    )
HBCC_rna_threshold <-    get_bonferroni(N_tests=sum(parquet_counts[mode=='rna' & cohort=='HBCC', N]), 
                            n_celltypes=length(celltypes), 
                            n_terms=4, 
                            alpha=0.05
                    )

# For each row of `parquet_counts`, import the parquet file
# Return only rows with significant Genetic 

thresholds <- CJ('COHORT'=c('HBCC','NABEC'), 'MODE'=c('rna','atac'))
thresholds[COHORT=='HBCC' &  MODE == 'rna', threshold := HBCC_rna_threshold]
thresholds[COHORT=='HBCC' &  MODE == 'atac', threshold := HBCC_atac_threshold]
thresholds[COHORT=='NABEC' & MODE == 'rna', threshold := NABEC_rna_threshold]
thresholds[COHORT=='NABEC' & MODE == 'atac', threshold := NABEC_atac_threshold]

import_parquet <- function(fn, method, signif_thresholds) {
    dat <- arrow::read_parquet(fn)
    setDT(dat)
    s <- strsplit(basename(fn), split='_|\\.')
    s <- unlist(s)
    cohort <- s[1]
    celltype <- s[2]
    mode <- s[3]
    chr <- s[7]
    threshold <- signif_thresholds[COHORT==cohort & MODE==mode, threshold]
    if(mode == 'rna') {
        if(cohort=='NABEC') { threshold <- NABEC_rna_threshold } else
        if(cohort=='HBCC') { threshold <- HBCC_rna_threshold } 
    } else if(mode == 'atac') {
        if(cohort=='NABEC') { threshold <- NABEC_atac_threshold } else
        if(cohort=='HBCC') { threshold <- HBCC_atac_threshold }
    }
    dat[, 'celltype' := celltype]
    dat[, 'cohort' := cohort]
    dat[, 'mode' := mode]
    dat[, 'chr' := chr]
    if(method=='nominal') {
        return(dat[pval_nominal <= threshold][])
    }else if(method=='interaction') {
        dat[, GenoAssoc := FALSE]
        dat[, Interaction := FALSE]
        dat[pval_gi < threshold, Interaction := TRUE]
        dat[pval_g < threshold, GenoAssoc := TRUE]

        return(dat[Interaction == TRUE | GenoAssoc==TRUE][])
    } else {
        cat('ERROR: Not yet implemented')
    }
}



if(!file.exists(summary_stats_fn)) {
   o <- foreach(file=filenames) %dopar% {
        import_parquet(fn=file,
                        method=args$method,
                        signif_thresholds=thresholds)
    }
    o <- rbindlist(o)
    # merge in gene IDs
    fwrite(o, file=summary_stats_fn, sep='\t')
    rm(o)
    gc()
}


qtls <- fread(summary_stats_fn)
qtls[, .N, by=list(celltype,cohort,mode)]

# if(!file.exists(summary_stats_fn)) {
#    o <- foreach(i=1:nrow(parquet_counts)) %dopar% {
#         row_i <- parquet_counts[i,]
#         cohort_i <- row_i[['cohort']]
#         celltype_i <- row_i[['celltype']]
#         mode_i <- row_i[['mode']]
#         chr_i <- row_i[['chr']]
#         if(cohort_i == 'HBCC') {
#             p_threshold_i <- hbcc_threshold
#         } else if(cohort_i == 'NABEC') {
#             p_threshold_i <- nabec_threshold
#         }
#         import_parquet(cohort = cohort_i,
#                         celltype = celltype_i,
#                         mode = mode_i,
#                         chr = chr_i,
#                         threshold = p_threshold_i
#         )

#     }   
#     o <- rbindlist(o)
#     # merge in gene IDs
#     fwrite(o, file=summary_stats_fn, sep='\t')
#     rm(o)
#     gc()
# }

# dat <- fread(summary_stats_fn)


# geno_assocs <- dat[GenoAssoc==TRUE, .N, by=list(celltype, cohort, mode)][order(celltype, cohort, mode)]
# interact_assoc <- dat[Interaction==TRUE, .N, by=list(celltype, cohort, mode)][order(celltype, cohort, mode)]
# all_assoc <- dat[Interaction==TRUE & GenoAssoc == TRUE]

# # Get variants for which associations are found
# if(!dir.exists('qtl_analysis')) {dir.create('qtl_analysis')}
# cat(unique(sort(dat$variant_id)), file='qtl_analysis/variant_ids.txt', sep='\n')


