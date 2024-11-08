#!/usr/bin/env Rscript

library(arrow)
library(data.table)
library(foreach)
library(doMC)
registerDoMC(cores=12)

cohorts <- c('hbcc','nabec')
celltypes <- c('ExN','InN','Astro','MG','Oligo','OPC','VC')
chrs <- paste0('chr',1:22)
modes <- c('rna','atac')

qtl_counts_fn <- 'QTL-output/tensorqtl-counts.tsv'
featurecounts_fn <- 'QTL-output/feature-counts.tsv'
if(!file.exists(qtl_counts_fn)) {
    files <- list.files('QTL-output', recursive=TRUE, include.dirs=TRUE, pattern='*.parquet', full.names=T)
    #dirs <- unique(unlist(lapply(files, dirname)))

    files.dt <- data.table(fn=files)
    files.dt[, c('d1','d2','bn') := tstrsplit(fn, split='/')[2:4]]
    files.dt[, 'pseudobulkmethod' := tstrsplit(d1, split='-')[1]]
    files.dt[, c('cohort','celltype','mode','tmp','interaction') := tstrsplit(d2, split='-')]
    files.dt[, tmp := NULL]
    files.dt[, 'chr' := tstrsplit(bn, split='\\.')[3]]
    files.dt[, c('d1','d2','bn') := NULL]

    #   For testing a smaller subset of rows before running 2000+
    #   files.dt <- files.dt[1:20]


    o <- foreach(fn=files.dt$fn, 
                pseudobulkmethod=files.dt$pseudobulkmethod,
                cohort=files.dt$cohort,
                celltype=files.dt$celltype,
                mode=files.dt$mode,
                interaction=files.dt$interaction,
                chr=files.dt$chr,
                .combine='rbind') %do% {
                    N <- nrow(arrow::read_parquet(fn))
                    dat.out <- data.table(pseudobulkmethod, cohort, celltype, mode, interaction, chr, N)
                    return(dat.out)
    }
    fwrite(o, file=qtl_counts_fn, quote=F, row.names=F, col.names=T, sep='\t')
} else {
    o <- fread(qtl_counts_fn)
}
       
o.summary <- o[, list('N'=sum(N), 'threshold'=0.05/sum(N)), by=list(pseudobulkmethod,cohort,interaction,mode)]

# Get significant results for non-interaction model
o.noninteraction <- foreach(i=which(o.summary$interaction == 'None'), .combine='rbind') %do% {
    pseudobulkmethod <- o.summary[[i,'pseudobulkmethod']]
    cohort <- o.summary[[i, 'cohort']]
    mode <- o.summary[[i, 'mode']]
    interaction <- o.summary[[i, 'interaction']]
    threshold <- o.summary[[i, 'threshold']]
    foreach(celltype=celltypes, .combine='rbind') %do% {
        foreach(chr=chrs, .combine='rbind') %do% {
            fn <- paste0('QTL-output/', pseudobulkmethod, '-nominal-20241031/', cohort, '-', celltype, '-', mode, '-interaction-', interaction, '/cis-nominal.cis_qtl_pairs.', chr, '.parquet')
            dat <- arrow::read_parquet(fn)
            setDT(dat)
            if(interaction == 'None') {
                dat <- dat[pval_nominal < threshold]
            } else {
                dat <- dat[pval_gi < threshold]
            }
            data.table(dat, fn, threshold, pseudobulkmethod, cohort, mode, interaction, threshold, celltype, chr)
        }
    }
}

fwrite(o.noninteraction, file='QTL-output/signif-noninteraction-qtls.tsv', quote=F, row.names=F, col.names=T, sep='\t')

# Get significant results for interaction model
o.interaction <- foreach(i=which(o.summary$interaction == 'Age'), .combine='rbind') %do% {
    pseudobulkmethod <- o.summary[[i,'pseudobulkmethod']]
    cohort <- o.summary[[i, 'cohort']]
    mode <- o.summary[[i, 'mode']]
    interaction <- o.summary[[i, 'interaction']]
    threshold <- o.summary[[i, 'threshold']]
    foreach(celltype=celltypes, .combine='rbind') %do% {
        foreach(chr=chrs, .combine='rbind') %do% {
            fn <- paste0('QTL-output/', pseudobulkmethod, '-nominal-20241031/', cohort, '-', celltype, '-', mode, '-interaction-', interaction, '/cis-nominal.cis_qtl_pairs.', chr, '.parquet')
            dat <- arrow::read_parquet(fn)
            setDT(dat)
            if(interaction == 'None') {
                dat <- dat[pval_nominal < threshold]
            } else {
                dat <- dat[pval_gi < threshold]
            }
            data.table(dat, fn, threshold, pseudobulkmethod, cohort, mode, interaction, threshold, celltype, chr)
        }
    }
}

fwrite(o.interaction, file='QTL-output/signif-interaction-qtls.tsv', quote=F, row.names=F, col.names=T, sep='\t')

# Go with sum just because
if(!file.exists(featurecounts_fn)) {
    files  <- list.files('/data/CARD_singlecell/cortex_subtype/', pattern='cpm_sum_log_pseudobulk_model_counts.csv', full.names=T, recursive=T)

    get_subset <- function(cohort, mode, celltypes=c('ExN','InN','Astro','MG','Oligo','OPC','VC')) {
        dat <- foreach(celltype=celltypes) %do% {
            pseudobulk <- fread(paste0('/data/CARD_singlecell/cortex_subtype/output/', mode, '/', celltype, '/pseudobulk/', cohort, '_cpm_sum_log_pseudobulk_model_counts.csv'))
            pseudobulk[,'celltype' := celltype]
            pseudobulk[,'cohort' := cohort]
            pseudobulk[,'mode' := mode]
            return(pseudobulk)
        }
        return(rbindlist(dat, fill=TRUE))
    }

    all_samples <- c()

    hbcc_rna <- get_subset('hbcc','rna')
    hbcc_samples <- grep('-ARC', colnames(hbcc_rna), value=TRUE)
    newnames <- paste0('HBCC_', hbcc_samples)
    newnames <- gsub('-ARC', '', newnames)
    all_samples <- c(all_samples, newnames)
    setnames(hbcc_rna, hbcc_samples, newnames)


    hbcc_atac <- get_subset('hbcc','atac')
    hbcc_samples <- grep('-ARC', colnames(hbcc_atac), value=TRUE)
    newnames <- paste0('HBCC_', hbcc_samples)
    newnames <- gsub('-ARC', '', newnames)
    all_samples <- c(all_samples, newnames)
    setnames(hbcc_atac, hbcc_samples, newnames)

    nabec_rna <- get_subset('nabec','rna')
    nabec_samples <- grep('-ARC', colnames(nabec_rna), value=TRUE)
    newnames <- gsub('-ARC', '', nabec_samples)
    all_samples <- c(all_samples, newnames)
    setnames(nabec_rna, nabec_samples, newnames)

    nabec_atac <- get_subset('nabec','atac')
    nabec_samples <- grep('-ARC', colnames(nabec_atac), value=TRUE)
    newnames <- gsub('-ARC', '', nabec_samples)
    all_samples <- c(all_samples, newnames)
    setnames(nabec_atac, nabec_samples, newnames)

    setnames(nabec_atac, 'peaks', 'feature')
    setnames(hbcc_atac, 'peaks', 'feature')
    setnames(nabec_rna, 'V1', 'feature')
    setnames(hbcc_rna, 'V1', 'feature')

    all_samples <- unique(all_samples)

    dat.counts <- rbindlist(list(hbcc_atac, hbcc_rna, nabec_atac, nabec_rna), fill=TRUE)

    featurecounts <- melt(dat.counts, measure.vars=all_samples, variable.name='sample', value.name='count')

    # Add in age info
    agetable <- fread('data/covariates.csv', select=c('sample','Sex','Age'))
    agetable[sample %like% 'UMARY-', sample:=gsub('-ARC','',sample)]
    agetable[sample %like% 'SH-', sample:=gsub('-ARC','',sample)]
    agetable[sample %like% 'KEN-', sample:=gsub('-ARC','',sample)]
    agetable[sample%like% '-ARC', sample := paste0('HBCC_', sample)]
    agetable[sample%like% '-ARC', sample:=gsub('-ARC','',sample)]

    featurecounts <- merge(featurecounts, agetable, by='sample')
    featurecounts[mode=='atac', feature := gsub(':','_',feature)]

    fwrite(featurecounts, file=featurecounts_fn, quote=F, row.names=F, col.names=T, sep='\t')
} else {
    featurecounts <- fread(featurecounts_fn)
}



# Get genotypes
genotypes.hbcc <- fread('genotypes/HBCC-alt-dosage.traw')
setkey(genotypes.hbcc, CHR, POS)
genotypes.nabec <- fread('genotypes/NABEC-alt-dosage.traw')
setkey(genotypes.nabec, CHR, POS)

plot_qtl <- function(hbcc_dt=genotypes.hbcc, nabec_dt=genotypes.nabec, feature_dt=featurecounts, variant_id, feature_name, method) {
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

    if(method=='Genotype') {
        g <- ggplot(data=dat, aes(x=alt_dosage, y=count)) +
            geom_boxplot() +
            geom_jitter(alpha=0.2) +
            facet_grid(celltype~cohort) +
            theme_few() +
            labs(x=variant_id, y=feature_name)
    } else if(method=='Age') {
        g <- ggplot(data=dat, aes(x=Age, y=count, group=alt_dosage, color=alt_dosage)) +
            geom_point() +
            geom_smooth(method='lm', se=FALSE) +
            facet_grid(celltype~cohort) +
            theme_few() +
            labs(title=variant_id, y=feature_name)
    } else if(method=='Sex') {
        # NYI
        g <- NULL
    }
    outname <- paste0('QTL-plots/', feature_name, '_', variant_id, '_', method, '.png')
    outname <- gsub(':','-', outname)
    ggsave(g, file=outname, width=15, height=20, units='cm')
    #return(g)

}

qtls <- fread('QTL-output/signif-noninteraction-qtls.tsv')
plot_qtl(feature_name='ARL17B', variant_id='chr17:46269852:C:A')

library(ggplot2)
library(ggthemes)

# eQTL
plot_qtl(feature_name='CHRM5', variant_id='chr15:33984369:G:A', method='Genotype')

plot_qtl(feature_name='ARL17B', variant_id='chr17:46269852:C:A', method='Genotype')
plot_qtl(feature_name='AC119673.2', variant_id='chr1:205813590:A:G', method='Genotype')    
plot_qtl(feature_name='HPR', variant_id='chr16:72060792:T:C', method='Genotype')
plot_qtl(feature_name='TXNL4B', variant_id='chr16:72061413:C:T', method='Genotype')
plot_qtl(feature_name='AL591686.2', variant_id='chr1:242147130:T:G', method='Genotype')    
plot_qtl(feature_name='HLA-DMA', variant_id='chr6:32946949:G:C', method='Genotype')
plot_qtl(feature_name='AC093151.8', variant_id='chr1:41148075:G:T', method='Genotype')   
plot_qtl(feature_name='APIP', variant_id='chr11:34884440:T:C', method='Genotype')
plot_qtl(feature_name='PKD1L3', variant_id='chr16:71998478:T:C', method='Genotype')
plot_qtl(feature_name='TOGARAM2', variant_id='chr2:28959915:C:T', method='Genotype')
plot_qtl(feature_name='PPM1K-DT', variant_id='chr4:88285020:T:C', method='Genotype')
plot_qtl(feature_name='LERFS', variant_id='chr9:62890385:C:G', method='Genotype')
plot_qtl(feature_name='PLD5', variant_id='chr1:242148188:A:G', method='Genotype')
plot_qtl(feature_name='PM20D1', variant_id='chr1:205813590:A:G', method='Genotype')
plot_qtl(feature_name='FAM198B-AS1', variant_id='chr4:158171481:A:G', method='Genotype')     


# age-interaction QTL
qtls <- fread('QTL-output/signif-interaction-qtls.tsv')

plot_qtl(feature_name='SMIM12', variant_id='chr1:34681413:C:T', method='Age')
plot_qtl(feature_name='GABRB2', variant_id='chr5:161277005:T:C', method='Age')
plot_qtl(feature_name='GABRB2', variant_id='chr5:161279237:C:G', method='Age')
plot_qtl(feature_name='PTPRM', variant_id='chr18:7469023:C:G', method='Age')
plot_qtl(feature_name='SLC7A11', variant_id='chr4:138217741:C:A', method='Age')
plot_qtl(feature_name='PLD5', variant_id='chr1:242147130:T:G', method='Age')
plot_qtl(feature_name='SLC18A1', variant_id='chr8:20175318:C:A', method='Age')

		
plot_qtl(feature_name='SMIM12',  variant_id='chr1:34681413:C:T', method='Age')
plot_qtl(feature_name='SMIM12',  variant_id='chr1:34681972:T:C', method='Age')
plot_qtl(feature_name='AL158064.1', variant_id='chr13:80024792:G:T', method='Age')
plot_qtl(feature_name='AL137781.1', variant_id='chr13:80024792:G:T', method='Age')
plot_qtl(feature_name='SMIM12',  variant_id='chr1:34681413:C:T', method='Age')
plot_qtl(feature_name='SMIM12',  variant_id='chr1:34681972:T:C', method='Age')
plot_qtl(feature_name='AL158064.1', variant_id='chr13:80024792:G:T', method='Age')
plot_qtl(feature_name='AL137781.1', variant_id='chr13:80024792:G:T', method='Age')
	
plot_qtl(feature_name='PRDM5',variant_id='chr4:120622647:G:C', method='Age')
plot_qtl(feature_name='SMIM12'   ,variant_id='chr1:34681972:T:C', method='Age')
plot_qtl(feature_name='AL137781.1'  ,variant_id='chr13:80024792:G:T', method='Age')
plot_qtl(feature_name='chr1_6601970_6603350'    ,variant_id='chr1:6589168:T:G', method='Age')
plot_qtl(feature_name='chr3_69385114_69386888'   ,variant_id='chr3:69377088:G:A', method='Age')
plot_qtl(feature_name='AL158064.1'  ,variant_id='chr13:80024792:G:T', method='Age')
plot_qtl(eature_name='chr1_241639703_241640734'  ,variant_id='chr1:241728320:G:A', method='Age')
plot_qtl(feature_name='chr19_58058668_58059926'  ,variant_id='chr19:58131851:G:A', method='Age')
plot_qtl(feature_name='chr6_16712779_16713673'   ,variant_id='chr6:16782405:C:T', method='Age')
plot_qtl(feature_name='chr20_461526_463296'    ,variant_id='chr20:402763:G:A', method='Age')
plot_qtl(feature_name='KDM6B'   ,variant_id='chr17:7819345:G:A', method='Age')
plot_qtl(feature_name='chr9_4489408_4491889'    ,variant_id='chr9:4407513:G:A', method='Age')
plot_qtl(feature_name='KDM6B'   ,variant_id='chr17:7819999:G:C', method='Age')
plot_qtl(feature_name='CPEB4'  ,variant_id='chr5:173826777:T:C', method='Age')
plot_qtl(feature_name='CPEB4'  ,variant_id='chr5:173822320:A:G', method='Age')
plot_qtl(feature_name='chr3_69385114_69386888'   ,variant_id='chr3:69379706:T:C', method='Age')
plot_qtl(feature_name='chr3_69385114_69386888'   ,variant_id='chr3:69374758:T:C', method='Age')
plot_qtl(feature_name='chr3_69385114_69386888'   ,variant_id='chr3:69375613:T:C', method='Age')
plot_qtl(feature_name='chr17_20189006_20190001'  ,variant_id='chr17:20275827:T:C', method='Age')
plot_qtl(feature_name='chr17_20189006_20190001'  ,variant_id='chr17:20177795:C:T', method='Age')
plot_qtl(feature_name='chr17_20189006_20190001'  ,variant_id='chr17:20200473:A:G', method='Age')
plot_qtl(feature_name='chr17_20189006_20190001'  ,variant_id='chr17:20250134:A:G', method='Age')
plot_qtl(feature_name='chr17_20189006_20190001'  ,variant_id='chr17:20130902:G:A', method='Age')
plot_qtl(feature_name='chr17_20189006_20190001'  ,variant_id='chr17:20133924:G:T', method='Age')
plot_qtl(feature_name='chr17_20189006_20190001'  ,variant_id='chr17:20245036:T:C', method='Age')
plot_qtl(feature_name='chr17_20189006_20190001'  ,variant_id='chr17:20271402:T:A', method='Age')
plot_qtl(feature_name='chr19_58182018_58183985'  ,variant_id='chr19:58131851:G:A', method='Age')
plot_qtl(feature_name='CPEB4'  ,variant_id='chr5:173826777:T:C', method='Age')
plot_qtl(feature_name='hr12_121132313_121133477' ,variant_id='chr12:121138184:G:A', method='Age')
plot_qtl(feature_name='chr17_20189006_20190001'  ,variant_id='chr17:20219862:G:A', method='Age')
plot_qtl(feature_name='chr17_20189006_20190001'  ,variant_id='chr17:20265746:A:C', method='Age')
plot_qtl(feature_name='chr6_16712779_16713673'   ,variant_id='chr6:16780117:G:A', method='Age')
plot_qtl(feature_name='CPEB4'  ,variant_id='chr5:173828289:G:A', method='Age')
plot_qtl(feature_name='ERC2'   ,variant_id='chr3:55551908:C:T', method='Age')
plot_qtl(feature_name='ERC2'   ,variant_id='chr3:55553610:C:T', method='Age')
plot_qtl(feature_name='ERC2'   ,variant_id='chr3:55574312:C:A', method='Age')
plot_qtl(feature_name='ERC2'   ,variant_id='chr3:55555150:C:T', method='Age')
plot_qtl(feature_name='ERC2'   ,variant_id='chr3:55551511:C:T', method='Age')
plot_qtl(feature_name='ERC2'   ,variant_id='chr3:55553514:C:T', method='Age')
plot_qtl(feature_name='ERC2'   ,variant_id='chr3:55572059:G:A', method='Age')
plot_qtl(feature_name='SMIM12'   ,variant_id='chr1:34681413:C:T', method='Age')
plot_qtl(feature_name='chr8_85462949_85464398'   ,variant_id='chr8:85494902:C:G', method='Age')
plot_qtl(feature_name='chr19_58058668_58059926'  ,variant_id='chr19:58157968:C:T', method='Age')
plot_qtl(feature_name='chr11_61947474_61947989'  ,variant_id='chr11:62032773:T:C', method='Age')
plot_qtl(feature_name='chr17_32274901_32275735'  ,variant_id='chr17:32216259:A:G', method='Age')
plot_qtl(feature_name='chr17_32274901_32275735'  ,variant_id='chr17:32206287:C:T', method='Age')