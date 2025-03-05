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

plt.colors <- c(
    'Oligodendrocyte' = '#945247ff',
    'Excitatory Neuron' = '#ff7500ff',
    'Inhibitory Neuron' = '#00a162ff',
    'Astrocyte' = '#0078b8ff',
    'Microglia' = '#ea0016ff',
    'OPC' = '#b735ffff',
    'Vascular Cell' = '#f26ec5ff'    
    )
    
featurefile <- paste0('QTL-output/sum-nominal-', date, '/cpm_sum_counts.tsv')

# Get all ATAC feature IDs and generate combined set
files <- list.files('/data/CARD_singlecell/cortex_subtype/output/atac', recursive=T, pattern='*_cpm_sum_pseudobulk_model_counts.csv', full.names=T)
atac_features <- foreach(file=files, .combine='rbind') %do% {
     fread(file, select='peaks')
}

setnames(atac_features, 'feature')
atac_features <- unique(atac_features)
atac_features[, c('chr','start','end') := tstrsplit(feature, split=':')]
atac_features[, start := as.numeric(start)]
atac_features[, end := as.numeric(end)]

library(GenomicRanges)
all_peaks <- makeGRangesFromDataFrame(atac_features, keep.extra.columns=TRUE)
merged_peaks <- reduce(all_peaks)
overlaps <- findOverlaps(all_peaks, merged_peaks)
final_peakset <- cbind(atac_features, as.data.table(merged_peaks[overlaps@to]@ranges))
setnames(final_peakset, c('feature_before','chr','start_before','end_before','start','end','width'))
final_peakset[, 'feature' := paste0(chr, '_', start, '_', end)]
final_peakset[, width_before := end_before - start_before]
final_peakset[, feature_before := gsub(':', '_', feature_before)]

signif_interaction <- fread(interaction_qtls_filename)
signif_noninteraction <- fread(qtls_filename)

signif_atac_features <- unique(sort(c(signif_interaction[mode=='atac'][['phenotype_id']], signif_noninteraction[mode=='atac'][['phenotype_id']])))
signif_rna_features <- unique(sort(c(signif_interaction[mode=='rna'][['phenotype_id']], signif_noninteraction[mode=='rna'][['phenotype_id']])))


# Get set of peaks that overlap with a significant peak
signif_merged_peaks <- final_peakset[feature_before %in% signif_atac_features]$feature
signif_peakset <- final_peakset[feature %in% signif_merged_peaks]
signif_unmerged_peaks <- signif_peakset$feature_before


signif_features <- unique(sort(c(signif_unmerged_peaks, signif_rna_features)))





# Get table of feature counts
if(file.exists(featurefile)) {
    featurecounts <- fread(featurefile)
} else {
    signif_features <- gsub('_', ':', signif_features)

    get_subset <- function(cohort, mode, celltypes=c('ExN','InN','Astro','MG','Oligo','OPC','VC'), signif_features) {
        dat <- foreach(celltype=celltypes) %do% {
            pseudobulk <- fread(paste0('/data/CARD_singlecell/cortex_subtype/output/', mode, '/', celltype, '/pseudobulk/', cohort, '_cpm_sum_pseudobulk_model_counts.csv'))
            # Remove trailing '-ARC'
            colnames(pseudobulk)[1] <- 'feature'
            pseudobulk <- pseudobulk[feature %in% signif_features]
            setnames(pseudobulk, gsub('-ARC', '', colnames(pseudobulk)))
            pseudobulk <- melt(pseudobulk, id.vars='feature', variable.name='sample', value.name='cpm_sum')
            pseudobulk[,'celltype' := celltype]
            pseudobulk[,'cohort' := cohort]
            pseudobulk[,'mode' := mode]
            pseudobulk[sample %like% '^[0-9]', sample := paste0('HBCC_', sample)]
            return(pseudobulk)
        }
        return(rbindlist(dat))
    }

    featurecounts <- foreach(cohort=c('hbcc','nabec'), .combine='rbind') %do% {
        foreach(mode=c('rna','atac'), .combine='rbind') %do% {
                get_subset(cohort=cohort, mode=mode, signif_features=signif_features)
            }
    }
    
    # Add in age info
    agetable <- fread('data/covariates.csv', select=c('sample','Sex','Age'))
    agetable[sample %like% 'UMARY-', sample:=gsub('-ARC','',sample)]
    agetable[sample %like% 'SH-', sample:=gsub('-ARC','',sample)]
    agetable[sample %like% 'KEN-', sample:=gsub('-ARC','',sample)]
    agetable[sample%like% '-ARC', sample := paste0('HBCC_', sample)]
    agetable[sample%like% '-ARC', sample:=gsub('-ARC','',sample)]

    featurecounts <- merge(featurecounts, agetable, by='sample')
    featurecounts[mode=='atac', feature := gsub(':','_',feature)]
    setkey(featurecounts, sample)
    fwrite(featurecounts, file=featurefile, quote=F, row.names=F, col.names=T, sep='\t')
}

# Add in new feature counts
setnames(featurecounts, 'feature', 'feature_before')
featurecounts <- merge(featurecounts, final_peakset[, .SD, .SDcols=c('feature_before','feature')], by='feature_before', all=TRUE)
featurecounts[mode=='rna', feature := feature_before]
#featurecounts[, 'feature_before' := NULL]
featurecounts <- featurecounts[!is.na(sample)]


# Save table of old-to-new feature mapping
feature_old_new <- copy(unique(featurecounts[, .SD, .SDcols=c('feature_before','feature')]))
setkey(feature_old_new, feature_before)

# Add column of new feature ID
signif_noninteraction <- merge(feature_old_new, signif_noninteraction, by.x='feature_before', by.y='phenotype_id')
signif_interaction <- merge(feature_old_new, signif_interaction, by.x='feature_before', by.y='phenotype_id')
 
signif_noninteraction[mode=='atac', c('chr','start','end') := tstrsplit(feature_before, split='_')]
signif_noninteraction[mode=='atac', start := as.numeric(start)]
signif_noninteraction[mode=='atac', end := as.numeric(end)]

signif_interaction[mode=='atac', c('chr','start','end') := tstrsplit(feature_before, split='_')]
signif_interaction[mode=='atac', start := as.numeric(start)]
signif_interaction[mode=='atac', end := as.numeric(end)]

# Load genotypes
genotypes.hbcc <- fread('genotypes/HBCC-alt-dosage.traw')
setkey(genotypes.hbcc, CHR, POS)
genotypes.nabec <- fread('genotypes/NABEC-alt-dosage.traw')
setkey(genotypes.nabec, CHR, POS)

# Aggregate cpm per newly mapped feature
featurecounts <- featurecounts[, list(N_peaks_merged=.N, 'cpm_sum'=sum(cpm_sum, na.rm=T)), by=list(sample, celltype, cohort, mode, Sex, Age, feature)]


# Define plotting function
plot_qtl <- function(hbcc_dt=genotypes.hbcc,
                    nabec_dt=genotypes.nabec,
                    feature_dt=featurecounts,
                        variant_id,
                        feature_name,
                        interaction,
                        outdir=tensorqtl_output_dir,
                        translate_dt=feature_old_new) {

    feature_name <- translate_dt[.(feature_name), feature]
    # Columns to exclude when finding sample IDs
    infocols <- c('CHR','SNP','(C)M','POS','COUNTED','ALT')
    variant_chr <- unlist(strsplit(variant_id, split=':'))[1]
    variant_pos <- unlist(strsplit(variant_id, split=':'))[2]
    variant_pos <- as.numeric(variant_pos)

    # Get hbcc genotypes in long format
    hbcc_sub <- hbcc_dt[CHR==variant_chr & POS == variant_pos]
    hbcc_samples <- setdiff(colnames(hbcc_sub),infocols)
    hbcc_genotypes <- melt(hbcc_sub, measure.vars=setdiff(colnames(genotypes.hbcc), infocols), variable.name='sample',value.name='alt_dosage')
    hbcc_genotypes <- hbcc_genotypes[!is.na(alt_dosage)]

    # get nabec genotypes in long format
    nabec_sub <- nabec_dt[CHR==variant_chr & POS == variant_pos]
    nabec_samples <- setdiff(colnames(nabec_sub),infocols)
    nabec_genotypes <- melt(nabec_sub, measure.vars=setdiff(colnames(genotypes.nabec), infocols), variable.name='sample',value.name='alt_dosage')
    nabec_genotypes <- nabec_genotypes[!is.na(alt_dosage)]

    # Merge genotypes together
    genotypes <- rbindlist(list(hbcc_genotypes, nabec_genotypes))
    genotypes[, alt_dosage := factor(alt_dosage, levels=c(0,1,2))]
    features <- feature_dt[feature==feature_name]

    # Merge genotypes with feature counts
    setkey(genotypes, sample)
    setkey(features, sample)
    dat <- merge(genotypes, features)
    dat[cohort=='hbcc', cohort := 'HBCC']
    dat[cohort=='nabec',cohort := 'NABEC']

    if(interaction=='None') {
        g <- ggplot(data=dat, aes(x=alt_dosage, y=cpm_sum)) +
            geom_boxplot(outliers=FALSE) +
            geom_jitter(alpha=0.2) +
            facet_grid(celltype~cohort) +
            theme_few() +
            labs(x=variant_id, y=feature_name)
    } else if(interaction=='Age') {
        g <- ggplot(data=dat, aes(x=Age, y=cpm_sum, group=alt_dosage, color=alt_dosage)) +
            geom_point() +
            geom_smooth(method='lm', se=FALSE) +
            facet_grid(celltype~cohort) +
            theme_few() +
            labs(title=variant_id, y=feature_name)
    }
    outname <- paste0(outdir, '/dosage-feature-plots/', feature_name, '_', variant_id, '_interaction-', interaction, '.png')
    outname <- gsub(':','-', outname)
    ggsave(g, file=outname, width=15, height=20, units='cm')
    #return(g)
}

cat("Ready to plot!\n")


# Plot (non-interaction) QTLs
plot_qtl(feature_name='CHRM5', variant_id='chr15:33984369:G:A', interaction='None')
plot_qtl(feature_name='DPP10', variant_id='chr2:113698175:A:G', interaction='None')

plot_qtl(feature_name='chr14_33749639_33750157', variant_id='chr14:33119657:G:A', interaction='None')

# Plot (interaction) QTLs

plot_qtl(feature_name='chr5_56374890_56375398', variant_id='chr5:56626084:G:A', interaction='Age')


