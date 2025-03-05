#!/usr/bin/env Rscript

library(data.table)
library(foreach)

# testline

o <- foreach(method=c('mean','sum'), .combine='rbind') %do% {

    files <- list.files(paste0('QTL-output/',method,'-nominal-20241124/'), recursive=T, full.names=T, pattern='*.tsv.gz')
        foreach(file=files, .combine='rbind') %do% {
            splitname <- strsplit(file, split='/')[[1]][4]
            splitname <- strsplit(splitname, split='-')[[1]]
            cohort <- splitname[1]
            celltype <- splitname[2]
            mode <- splitname[3]
            interaction <- splitname[5]
            if(interaction=='None') {
                dat <- fread(file)
            } else {
                return(NULL)
            }
            dat[, 'cohort' := cohort]
            dat[, 'celltype' := celltype]
            dat[, 'mode' := mode]
            dat[, 'interaction' := interaction]
            dat[, 'method' := method]
    }
}
dcast(o[, .N, by=list(cohort, celltype, mode, interaction, method)], cohort+celltype+mode+interaction~method, value.var='N')


method <- 'sum'
files <- list.files(paste0('QTL-output/',method,'-nominal-20241124/'), recursive=T, full.names=T, pattern='*.tsv.gz')

o.interaction <- foreach(file=files, .combine='rbind') %do% {
    splitname <- strsplit(file, split='/')[[1]][4]
    splitname <- strsplit(splitname, split='-')[[1]]
    cohort <- splitname[1]
    celltype <- splitname[2]
    mode <- splitname[3]
    interaction <- splitname[5]
    if(interaction=='Age') {
        dat <- fread(file)
    } else {
        return(NULL)
    }
    dat[, 'cohort' := cohort]
    dat[, 'celltype' := celltype]
    dat[, 'mode' := mode]
    dat[, 'interaction' := interaction]
}

o.interaction[, .N, by=list(cohort, celltype, mode, interaction)]


## Look at Allele frequency spectrum

library(data.table)
library(ggplot2)
library(foreach)



# HBCC



# NABEC

o <- foreach(cohort=c('HBCC','NABEC'), .combine='rbind') %do% {
    fname <- paste0('genotypes/', cohort, '-alt-dosage.traw')
    dat <- fread(fname)
    sids <- colnames(dat)[7:length(colnames(dat))]
    dat[, nRef := apply(.SD, 1, function(x) 2*sum(x==0, na.rm=T)+sum(x==1, na.rm=T)), .SDcols=sids]
    dat[, nAlt := apply(.SD, 1, function(x) 2*sum(x==2, na.rm=T)+sum(x==1, na.rm=T)), .SDcols=sids]
    dat[, p := nRef / ( nRef + nAlt)]
    set(dat, , sids, NULL)
    dat[p > 0.5, p := 1 - p]
    dat <- dat[, .SD, .SDcols=c('nRef','nAlt','p')]
    dat[, 'cohort' := cohort][]
    return(dat)
}

o <- o[p > 0.05]

ggplot(o, aes(x=p)) + geom_histogram(bins=20) + facet_wrap(~cohort, scales='free_y') + labs(y='count', x='minor allele frequency')

ggplot(o, aes(x=cohort, y=p)) + geom_violin()