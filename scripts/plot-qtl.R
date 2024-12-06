#!/usr/bin/env Rscript

library(data.table)
library(foreach)
library(doMC)
library(parallel)
registerDoMC(cores=4)

library(ggplot2)
library(ggthemes)

setwd('/data/CARD_singlecell/users/wellerca/pfc-atlas-qtl')

featurefile <- 'data/cpm_sum_counts.RDS'

# Get all ATAC features

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

signif_interaction <- fread('significant-interaction-qtls.tsv.gz')
signif_interaction <- signif_interaction[pseudobulk_method=='sum']

signif_noninteraction <- fread('significant-noninteraction-qtls.tsv.gz')
signif_noninteraction <- signif_noninteraction[pseudobulk_method=='sum']

signif_atac_features <- unique(sort(c(signif_interaction[mode=='atac'][['phenotype_id']], signif_noninteraction[mode=='atac'][['phenotype_id']])))
signif_rna_features <- unique(sort(c(signif_interaction[mode=='rna'][['phenotype_id']], signif_noninteraction[mode=='rna'][['phenotype_id']])))


# Get set of peaks that overlap with a significant peak
signif_merged_peaks <- final_peakset[feature_before %in% signif_atac_features]$feature
signif_peakset <- final_peakset[feature %in% signif_merged_peaks]
signif_unmerged_peaks <- signif_peakset$feature_before


signif_features <- unique(sort(c(signif_unmerged_peaks, signif_rna_features)))





# Get table of feature counts
if(file.exists(featurefile)) {
    featurecounts <- readRDS(featurefile)
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
    saveRDS(featurecounts, file=featurefile)
}

# Add in new feature counts
setnames(featurecounts, 'feature', 'feature_before')
featurecounts <- merge(featurecounts, final_peakset[, .SD, .SDcols=c('feature_before','feature')], by='feature_before', all=TRUE)
featurecounts[mode=='rna', feature := feature_before]
featurecounts[, 'feature_before' := NULL]
featurecounts <- featurecounts[!is.na(sample)]

# What happened to feature 'chr6:4134024:4134696'

# Load genotypes
genotypes.hbcc <- fread('genotypes/HBCC-alt-dosage.traw')
setkey(genotypes.hbcc, CHR, POS)
genotypes.nabec <- fread('genotypes/NABEC-alt-dosage.traw')
setkey(genotypes.nabec, CHR, POS)

# Define plotting function
plot_qtl <- function(hbcc_dt=genotypes.hbcc,
                    nabec_dt=genotypes.nabec,
                    feature_dt=featurecounts,
                        variant_id,
                        feature_name,
                        interaction,
                        translate_dt=feature_old_new) {

    feature_name <- feature_old_new[.(feature_name), feature]
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
    outname <- paste0('QTL-plots/', feature_name, '_', variant_id, '_interaction-', interaction, '.png')
    outname <- gsub(':','-', outname)
    ggsave(g, file=outname, width=15, height=20, units='cm')
    #return(g)
}

# Load QTLs
qtls <- fread('significant-noninteraction-qtls.tsv.gz')[pseudobulk_method == 'sum']
interaction_qtls <- fread('significant-interaction-qtls.tsv.gz')[pseudobulk_method == 'sum']

# Save table of old-to-new feature mapping
feature_old_new <- copy(unique(featurecounts[, .SD, .SDcols=c('feature_before','feature')]))

# Add column of new feature ID
qtls <- merge(feature_old_new, qtls, by.x='feature_before', by.y='phenotype_id')
interaction_qtls <- merge(feature_old_new, interaction_qtls, by.x='feature_before', by.y='phenotype_id')
 

qtls[, c('chr','start','end') := tstrsplit(feature_before, split='_')]
qtls[, start := as.numeric(start)]
qtls[, end := as.numeric(end)]

interaction_qtls[, c('chr','start','end') := tstrsplit(feature_before, split='_')]
interaction_qtls[, start := as.numeric(start)]
interaction_qtls[, end := as.numeric(end)]

# Aggregate cpm per newly mapped feature
featurecounts <- featurecounts[, list(N_peaks_merged=.N, 'cpm_sum'=sum(cpm_sum, na.rm=T)), by=list(sample, celltype, cohort, mode, Sex, Age, feature)]




# Plot (non-interaction) QTLs
plot_qtl(feature_name='CHRM5', variant_id='chr15:33984369:G:A', interaction='None')
plot_qtl(feature_name='DPP10', variant_id='chr2:113698175:A:G', interaction='None')


# Plot (age interaction) QTLs
plot_qtl(feature_name='chr6_4134726_4136360', variant_id='chr6:4945769:G:A', interaction='Age')
