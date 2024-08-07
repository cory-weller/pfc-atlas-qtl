#!/usr/bin/env Rscript

library(data.table)

snakemake_samplefile <- '/data/CARD_singlecell/brain_atlas_rna/input/samples.csv'

dat <- fread(snakemake_samplefile)

dat[, 'bampath' := paste0('/data/CARD_singlecell/Brain_atlas/', toupper(cohort), '_multiome/', batch, '/Multiome/', sample, '/outs/gex_possorted_bam.bam')]
dat[, cohort := toupper(cohort)]

fwrite(dat, file='samples.tsv', quote=F, row.names=F, col.names=T, sep='\t')

