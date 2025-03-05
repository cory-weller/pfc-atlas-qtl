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
qtls_filename <- paste0(top_output_dir, '/significant-QTLs.tsv')
interaction_qtls_filename <- paste0(top_output_dir, '/significant-interaction-QTLs.tsv')

if(!file.exists(qtls_filename) | !file.exists(interaction_qtls_filename)) {
   cat("ERROR: Summarized QTL files not yet generated, did you run plot-N-qtl-lambda.R ?\n")
   exit(1)
}


signif_interaction <- fread(interaction_qtls_filename)
signif_noninteraction <- fread(qtls_filename)

rna_variants <- unique(signif_noninteraction[mode=='rna', variant_id])
atac_variants <- unique(signif_noninteraction[mode=='atac', variant_id])
dual_mode_variants <- intersect(rna_variants, atac_variants)
rna_only_variants <- setdiff(rna_variants, atac_variants)
atac_only_variants <- setdiff(atac_variants, rna_variants)

length(dual_mode_variants)
length(rna_only_variants)
length(atac_only_variants)

signif_noninteraction[variant_id %in% dual_mode_variants, 'dual_QTL_any' := TRUE]
signif_noninteraction[is.na(dual_QTL_any), dual_QTL_any := FALSE]

for(.celltype in unique(signif_noninteraction$celltype)) {
    dual_mode_variants <- intersect(signif_noninteraction[mode=='rna' & celltype == .celltype, variant_id],
                                signif_noninteraction[mode=='atac' & celltype == .celltype, variant_id])
    dual_per_celltype <- copy(signif_noninteraction[celltype == .celltype & variant_id %in% dual_mode_variants])[order(variant_id)]
    fwrite(dual_per_celltype, file=paste0(top_output_dir, '/dualQTL-', .celltype, '.tsv'), quote=F, row.names=F, col.names=T, sep='\t')
    #signif_noninteraction[, paste0('dual_QTL_', .celltype) := FALSE]
    #signif_noninteraction[celltype == .celltype & variant_id %in% dual_mode_variants, paste0('dual_QTL_', .celltype) := TRUE][]
}




eQTL <- signif_interaction[mode=='rna' & variant_id %in% rna_only_variants]
caQTL <- signif_interaction[mode=='atac' & variant_id %in% atac_only_variants]



fwrite(eQTL, file=paste0(top_output_dir, '/eQTL-only.tsv'), quote=F, row.names=F, col.names=T, sep='\t')
fwrite(caQTL, file=paste0(top_output_dir, '/caQTL-only.tsv'), quote=F, row.names=F, col.names=T, sep='\t')
