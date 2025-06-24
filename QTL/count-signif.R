#!/usr/bin/env Rscript

library(data.table)

files <- list.files(pattern='pval-correction-mapping.txt', full.names=T, recursive=T)

library(foreach)

o <- foreach(file=files, .combine='rbind') %do% {
    dat <- fread(file, header=T)
    batch <- basename(dirname(file))
    cohort <- unlist(strsplit(batch, split='-'))[1]
    celltype <- unlist(strsplit(batch, split='-'))[2]
    mode <- unlist(strsplit(batch, split='-'))[3]
    dat[, 'cohort' := cohort]
    dat[, 'celltype' := celltype]
    dat[, 'mode' := mode]
    return(dat[])
}

# Number passing BH correction
o1 <- o[, .N, by=list(cohort, celltype, mode)]
o1[, method := 'BH']
    
    
# Number passing Bonferroni Correction
o2 <- o[p_Bonferroni < 0.05, .N, by=list(cohort, celltype, mode)]
o2[, method := 'Bonferroni']

dat <- rbindlist(list(o1, o2))

options(scipen=999)

library(ggplot2)
library(ggrepel)
ggplot(dat, aes(x=celltype, fill=cohort, y=N)) + geom_bar(position='dodge', stat='identity') +
    facet_grid(method~., scales='free_y') +
    labs(title='p-value corrected significant caQTLs') +
    geom_label_repel(aes(label=N))