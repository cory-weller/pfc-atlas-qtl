#!/usr/bin/env Rscript

library(data.table)
library(foreach)
library(doMC)
library(parallel)
registerDoMC(cores=4)

library(ggplot2)
library(ggthemes)

featurefile <- 'data/cpm_sum_counts.RDS'

# Get table of feature counts
if(file.exists(featurefile)) {
    featurecounts <- readRDS(featurefile)
} else {
    # get list of all significant features
    signif_interaction <- fread('significant-interaction-qtls.tsv.gz')
    signif_interaction <- signif_interaction[['phenotype_id']]
    signif_noninteraction <- fread('significant-noninteraction-qtls.tsv.gz')
    signif_noninteraction <- signif_noninteraction[['phenotype_id']]
    signif_features <- unique(sort(c(signif_interaction, signif_noninteraction)))


    # Subset pseudobulk counts to these features

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




# Load genotypes
genotypes.hbcc <- fread('genotypes/HBCC-alt-dosage.traw')
setkey(genotypes.hbcc, CHR, POS)
genotypes.nabec <- fread('genotypes/NABEC-alt-dosage.traw')
setkey(genotypes.nabec, CHR, POS)

# Define plotting function
plot_qtl <- function(hbcc_dt=genotypes.hbcc,
                    nabec_dt=genotypes.nabec,
                    feature_dt=featurecounts,
                        variant_id, feature_name, interaction) {
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

    genotypes <- rbindlist(list(hbcc_genotypes, nabec_genotypes))
    genotypes[, alt_dosage := factor(alt_dosage, levels=c(0,1,2))]
    features <- feature_dt[feature==feature_name]

    # Merge
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

# Plot (non-interaction) QTLs
 plot_qtl(feature_name='CHRM5', variant_id='chr15:33984369:G:A', interaction='None')
 plot_qtl(feature_name='DPP10', variant_id='chr2:113698175:A:G', interaction='None')


# Plot (age interaction) QTLs
 plot_qtl(feature_name='CHRM5', variant_id='chr15:33984369:G:A', interaction='Age')
