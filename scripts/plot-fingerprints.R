#!/usr/bin/env Rscript

library(data.table)
library(foreach)
library(ggplot2)
library(ggthemes)

o <- foreach(file=list.files('fingerprints', pattern='*.zip', full.names=T), .combine='rbind') %do% {
    sample <- strsplit(basename(file), split='-fingerprinting')[[1]][1]
    unzip_cmd <- paste('unzip -p', file, paste0(sample, '.fingerprinting.fingerprinting_summary_metrics'))
    dat <- fread(cmd=unzip_cmd, skip='READ_GROUP')
}

o <- o[order(LOD_EXPECTED_SAMPLE)]
o[, SAMPLE := factor(SAMPLE, levels=o$SAMPLE)]

g <- ggplot(o, aes(x=SAMPLE, y=LOD_EXPECTED_SAMPLE)) +
    geom_point() +
    theme_few() +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    geom_hline(yintercept=0, linetype='dashed', color='gray') +
    labs(x='NABEC Sample', y='LOD (more positive is better)', title='Fingerprinting of snRNA x WGS Genotypes')

ggsave(g, file='plots/NABEC-fingerprinting.png', width=55, height=12, units='cm')