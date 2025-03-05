#!/usr/bin/env Rscript

library(data.table)
library(foreach)
library(doMC)
library(parallel)
registerDoMC(cores=4)

library(ggplot2)
library(ggthemes)

setwd('/data/CARD_singlecell/users/wellerca/pfc-atlas-qtl')

# Combine significant tables

date <- '20241210'

top_output_dir <- paste0('QTL-output/sum-nominal-', date)
qc_output_dir <- paste0(top_output_dir, '/QC-plots')
tensorqtl_output_dir <- paste0('QTL-output/sum-nominal-', date, '/tensorqtl/')

files <- list.files(tensorqtl_output_dir, pattern='significant.tsv.gz', full.names=T, recursive=T)
qtl_files <- grep(files, pattern='-None/', value=T)
interaction_files <- grep(files, pattern='-Age/', value=T)

qtls_filename <- paste0(top_output_dir, '/significant-QTLs.tsv')
interaction_qtls_filename <- paste0(top_output_dir, '/significant-interaction-QTLs.tsv')

if(!file.exists(qtls_filename)) {
    o <- foreach(file=qtl_files, .combine='rbind') %do% {
        dat <- fread(file)
        params <- basename(gsub('/significant.tsv.gz', '', file))
        params <- unlist(strsplit(params, split='-'))
        cohort <- params[1]
        celltype <- params[2]
        mode <- params[3]
        interaction <- params[5]
        dat[, 'cohort' := cohort]
        dat[, 'celltype' := celltype]
        dat[, 'mode' := mode]
        dat[, 'interaction' := interaction]
        return(dat)
    }

    fwrite(o, file=qtls_filename, quote=F, row.names=F, col.names=T, sep='\t')
}

if(!file.exists(interaction_qtls_filename)) {
    o <- foreach(file=interaction_files, .combine='rbind') %do% {
        dat <- fread(file)
        params <- basename(gsub('/significant.tsv.gz', '', file))
        params <- unlist(strsplit(params, split='-'))
        cohort <- params[1]
        celltype <- params[2]
        mode <- params[3]
        interaction <- params[5]
        dat[, 'cohort' := cohort]
        dat[, 'celltype' := celltype]
        dat[, 'mode' := mode]
        dat[, 'interaction' := interaction]
        return(dat)
    }

    fwrite(o, file=interaction_qtls_filename, quote=F, row.names=F, col.names=T, sep='\t')
}

interaction_qtls <- fread(interaction_qtls_filename)
qtls <- fread(qtls_filename)

qtl_counts <- rbindlist(list(qtls[, .N, by=list(celltype, cohort, mode, interaction)],
            interaction_qtls[, .N, by=list(celltype, cohort, mode, interaction)]))


qtl_counts[celltype == 'ExN', celltype := 'Excitatory Neuron']
qtl_counts[celltype == 'InN', celltype := 'Inhibitory Neuron']
qtl_counts[celltype == 'MG', celltype := 'Microglia']
qtl_counts[celltype == 'Oligo', celltype := 'Oligodendrocyte']
qtl_counts[celltype == 'OPC', celltype := 'OPC']
qtl_counts[celltype == 'VC', celltype := 'Vascular Cell']
qtl_counts[celltype == 'Astro', celltype := 'Astrocyte']
qtl_counts[, celltype := factor(celltype, levels=c('Excitatory Neuron','Inhibitory Neuron','Astrocyte','Microglia','Oligodendrocyte','OPC','Vascular Cell'))]

qtl_counts[, cohort_total := sum(N), by=list(cohort, mode, interaction)]
qtl_counts[, Proportion := N/cohort_total]

plt.colors <- c(
    'Oligodendrocyte' = '#945247ff',
    'Excitatory Neuron' = '#ff7500ff',
    'Inhibitory Neuron' = '#00a162ff',
    'Astrocyte' = '#0078b8ff',
    'Microglia' = '#ea0016ff',
    'OPC' = '#b735ffff',
    'Vascular Cell' = '#f26ec5ff'    
    )
    
plot_nqtls <- function(DT, clrs) {
    ggplot(DT, aes(x=cohort, y=N, fill=celltype)) +
    geom_bar(stat='identity') +
    theme_few(14) +
    labs(fill='Cell Type', x='Cohort', y='Number', title='Significant QTLs (BH correction)') +
    facet_grid(mode~interaction, scales='free_y', labeller=labeller(.rows=label_both,.cols=label_both)) +
    guides(fill = guide_legend()) +
    scale_fill_manual(values=clrs)
}
options(scipen = 999) 
g1 <- plot_nqtls(qtl_counts, plt.colors)
ggsave(g1, file=paste0(qc_output_dir, '/N-QTLs.png'), width=20, height=20, units='cm')
ggsave(g1, file=paste0(qc_output_dir, '/N-QTLs.svg'), width=20, height=20, units='cm')

plot_propqtls <- function(DT, clrs) {
    ggplot(DT, aes(x=cohort, y=Proportion, fill=celltype)) +
    geom_bar(stat='identity') +
    theme_few(14) +
    labs(fill='Cell Type', x='Cohort', y='Proportion', title='Significant QTLs (BH correction)') +
    facet_grid(mode~interaction, scales='free_y', labeller=labeller(.rows=label_both,.cols=label_both)) +
    guides(fill = guide_legend()) +
    scale_fill_manual(values=clrs)
}

g2 <- plot_propqtls(qtl_counts, plt.colors)
ggsave(g2, file=paste0(qc_output_dir, '/Proportion-QTLs.png'), width=20, height=20, units='cm')
ggsave(g2, file=paste0(qc_output_dir, '/Proportion-QTLs.svg'), width=20, height=20, units='cm')


# Plot lambda values
files <- list.files(tensorqtl_output_dir, pattern='lambda.txt', full.names=T, recursive=T)
lambdas <- foreach(file=files, .combine='rbind') %do% {
    lambda <- as.numeric(readLines(file))
    params <- basename(gsub('/lambda.txt', '', file))
    params <- unlist(strsplit(params, split='-'))
    cohort <- params[1]
    celltype <- params[2]
    mode <- params[3]
    interaction <- params[5]
    return(data.table(cohort, celltype, mode, interaction, lambda))
}


lambdas[celltype == 'ExN', celltype := 'Excitatory Neuron']
lambdas[celltype == 'InN', celltype := 'Inhibitory Neuron']
lambdas[celltype == 'MG', celltype := 'Microglia']
lambdas[celltype == 'Oligo', celltype := 'Oligodendrocyte']
lambdas[celltype == 'OPC', celltype := 'OPC']
lambdas[celltype == 'VC', celltype := 'Vascular Cell']
lambdas[celltype == 'Astro', celltype := 'Astrocyte']
lambdas[, celltype := factor(celltype, levels=c('Excitatory Neuron','Inhibitory Neuron','Astrocyte','Microglia','Oligodendrocyte','OPC','Vascular Cell'))]




g <- ggplot(lambdas, aes(x=cohort, fill=celltype, y=lambda)) + 
    geom_bar(color='black', stat='identity', position=position_dodge2(preserve='single', padding=0)) + 
    facet_grid(mode~interaction) +
    geom_hline(yintercept=1, linetype='dashed', alpha=0.6) +
    theme_few(12) +
    scale_fill_manual(values=plt.colors)

ggsave(g, file=paste0(qc_output_dir, '/lambda-gc.png'), width=25, height=18, units='cm')
ggsave(g, file=paste0(qc_output_dir, '/lambda-gc.svg'), width=25, height=18, units='cm')
