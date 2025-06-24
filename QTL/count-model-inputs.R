#!/usr/bin/env Rscript

library(data.table)

files <- list.files(pattern='model_summary.txt', full.names=T, recursive=T)

library(foreach)

o <- foreach(file=files, .combine='rbind') %do% {
    dat <- fread(file, header=F)
    dat[,V1 := gsub(':','',V1)]
    dat <- as.data.table(t(dat))
    setnames(dat, as.character(dat[1,]))
    dat <- dat[2,]
    batch <- basename(dirname(file))
    cohort <- unlist(strsplit(batch, split='-'))[1]
    celltype <- unlist(strsplit(batch, split='-'))[2]
    mode <- unlist(strsplit(batch, split='-'))[3]
    dat[, 'cohort' := cohort]
    dat[, 'celltype' := celltype]
    dat[, 'mode' := mode]
    return(dat[])
}

o[, Samples := as.numeric(Samples)]
o[, Features := as.numeric(Features)]
o[, Variants := as.numeric(Variants)]
setkey(o, cohort, mode, celltype)

fwrite(o, file='all-model-inputs.csv', sep=',', quote=F, row.names=F, col.names=T)
