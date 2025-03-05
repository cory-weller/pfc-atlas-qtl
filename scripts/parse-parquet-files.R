#!/usr/bin/env Rscript

library(arrow)
library(data.table)
library(foreach)
library(doMC)
library(ggplot2)
library(ggthemes)
registerDoMC(cores=12)


method <- Sys.getenv('QTLMETHOD')
N <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

fdr <- 0.05

qtldir <- paste0('QTL-output/', method, '-nominal-20241210')
qtldirs <- list.dirs(qtldir, recursive=FALSE)
chosendir <- qtldirs[N]

signif_qtls <- paste0(chosendir, '/significant.tsv.gz')
lambda_file <- gsub('significant.tsv.gz', 'lambda.txt', signif_qtls)


files <- list.files(chosendir, recursive=TRUE, include.dirs=TRUE, pattern='*.parquet', full.names=T)
myvars <- unlist(strsplit(basename(chosendir), split='-'))
cohort <- myvars[1]
celltype <- myvars[2]
mode <- myvars[3]
interaction <- myvars[5]

get_lambda <- function(ps) {
    chisq <- qchisq(1 - ps, 1)
    lambda <- median(chisq) / qchisq(0.5, 1)
    lambda
}


# Iterate over parquet files and bind together
o <- foreach(fn=files, .combine='rbind') %do% {
                dat <- arrow::read_parquet(fn)
                setDT(dat)
                return(dat)
}

# Adjust p-vals and filter
if(interaction == 'None') {
    o[, pval_nominal_bh := p.adjust(pval_nominal, method='BH')]
    mylambda <- get_lambda(o$pval_nominal)
    o <- o[pval_nominal_bh < fdr]
} else {
    o[, pval_gi_bh := p.adjust(pval_gi, method='BH')]
    mylambda <- get_lambda(o$pval_gi)
    o <- o[pval_gi_bh < fdr]
}

fwrite(o, file=signif_qtls, quote=F, row.names=F, col.names=T, sep='\t')
writeLines(as.character(mylambda), con=lambda_file)



quit(status=0)

