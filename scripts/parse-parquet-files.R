#!/usr/bin/env Rscript

library(arrow)
library(data.table)
library(foreach)
library(doMC)
library(ggplot2)
library(ggthemes)
registerDoMC(cores=12)


method <- Sys.getenv('QTLMETHOD')
N <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

fdr <- 0.05

qtldir <- paste0('QTL-output/', method, '-nominal-20241124')
qtldirs <- list.dirs(qtldir, recursive=FALSE)
chosendir <- qtldirs[N]

signif_qtls <- paste0(chosendir, '/significant.tsv.gz')


files <- list.files(chosendir, recursive=TRUE, include.dirs=TRUE, pattern='*.parquet', full.names=T)
myvars <- unlist(strsplit(basename(chosendir), split='-'))
cohort <- myvars[1]
celltype <- myvars[2]
mode <- myvars[3]
interaction <- myvars[5]

inflation <- function(ps) {
    chisq <- qchisq(1 - ps, 1)
    lambda <- median(chisq) / qchisq(0.5, 1)
    lambda
}

if(FALSE) {
    # Iterate over parquet files and bind together
    foreach(method=c('sum'), .combine='rbind') %do% {
        qtldir <- paste0('QTL-output/', method, '-nominal-20241124')
        qtldirs <- list.dirs(qtldir, recursive=FALSE)
        foreach(mydir=qtldirs, .combine='rbind') %do% {

            files <- list.files(mydir, recursive=TRUE, include.dirs=TRUE, pattern='*.parquet', full.names=T)
            myvars <- unlist(strsplit(basename(mydir), split='-'))
            cohort <- myvars[1]
            celltype <- myvars[2]
            mode <- myvars[3]
            interaction <- myvars[5]

            o <- foreach(fn=files, .combine='rbind') %do% {
                            dat <- arrow::read_parquet(fn)
                            setDT(dat)
                            samplesize <- floor(nrow(dat) * 1e-4)
                            return(dat[sample(.N, size=samplesize)])
            }

            png(paste0('QQ/', cohort, '-', celltype, '-', mode, '.png'))
            qq(o$pval_nominal)
            dev.off()

            mylambda <- inflation(o$pval_nominal)
            return(data.table(cohort, celltype, method, mode, interaction, 'lambda'=mylambda))
        }   
    }
}

# Iterate over parquet files and bind together
o <- foreach(fn=files, .combine='rbind') %do% {
                dat <- arrow::read_parquet(fn)
                setDT(dat)
                return(dat)
}

# Adjust p-vals and filter
if(interaction == 'None') {
    o[, pval_nominal_bh := p.adjust(pval_nominal, method='BH')]
    o <- o[pval_nominal_bh < fdr]
} else {
    o[, pval_i_bh := p.adjust(pval_i, method='BH')]
    o[, pval_gi_bh := p.adjust(pval_gi, method='BH')]
    o <- o[pval_gi_bh < fdr]
}
fwrite(o, file=signif_qtls, quote=F, row.names=F, col.names=T, sep='\t')


quit(status=0)


# Get genotypes
genotypes.hbcc <- fread('genotypes/HBCC-alt-dosage.traw')
setkey(genotypes.hbcc, CHR, POS)
genotypes.nabec <- fread('genotypes/NABEC-alt-dosage.traw')
setkey(genotypes.nabec, CHR, POS)

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
        g <- ggplot(data=dat, aes(x=alt_dosage, y=count)) +
            geom_boxplot() +
            geom_jitter(alpha=0.2) +
            facet_grid(celltype~cohort) +
            theme_few() +
            labs(x=variant_id, y=feature_name)
    } else if(interaction=='Age') {
        g <- ggplot(data=dat, aes(x=Age, y=count, group=alt_dosage, color=alt_dosage)) +
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

interaction <- 'None'

if(interaction=='None') {
    qtls <- fread('significant-noninteraction-qtls.tsv.gz')
} else if(interaction=='Age') {
    qtls <- fread('significant-interaction-qtls.tsv.gz')
}

if(!file.exists('significant-atac-features.tsv')) {
    qtls <- fread('significant-noninteraction-qtls.tsv.gz')
    atac_noninteraction <- featurecounts[feature %in% unique(qtls[mode=='atac']$phenotype_id)]
    qtls <- fread('significant-interaction-qtls.tsv.gz')
    atac_interaction <- featurecounts[feature %in% unique(qtls[mode=='atac']$phenotype_id)]
    atac_signif_features <- unique(rbindlist(list(atac_noninteraction, atac_interaction)))
    fwrite(atac_signif_features, file='significant-atac-features.tsv', quote=F, row.names=F, col.names=T, sep='\t')
} else {
    atac_signif_features <- fread('significant-atac-features.tsv')
}
plot_qtl(feature_name='ARL17B', variant_id='chr17:46269852:C:A', interaction=interaction)



# # eQTL
# plot_qtl(feature_name='CHRM5', variant_id='chr15:33984369:G:A', interaction='None')

# plot_qtl(feature_name='ARL17B', variant_id='chr17:46269852:C:A', interaction='None')
# plot_qtl(feature_name='AC119673.2', variant_id='chr1:205813590:A:G', interaction='None')    
# plot_qtl(feature_name='HPR', variant_id='chr16:72060792:T:C', interaction='None')
# plot_qtl(feature_name='TXNL4B', variant_id='chr16:72061413:C:T', interaction='None')
# plot_qtl(feature_name='AL591686.2', variant_id='chr1:242147130:T:G', interaction='None')    
# plot_qtl(feature_name='HLA-DMA', variant_id='chr6:32946949:G:C', interaction='None')
# plot_qtl(feature_name='AC093151.8', variant_id='chr1:41148075:G:T', interaction='None')   
# plot_qtl(feature_name='APIP', variant_id='chr11:34884440:T:C', interaction='None')
# plot_qtl(feature_name='PKD1L3', variant_id='chr16:71998478:T:C', interaction='None')
# plot_qtl(feature_name='TOGARAM2', variant_id='chr2:28959915:C:T', interaction='None')
# plot_qtl(feature_name='PPM1K-DT', variant_id='chr4:88285020:T:C', interaction='None')
# plot_qtl(feature_name='LERFS', variant_id='chr9:62890385:C:G', interaction='None')
# plot_qtl(feature_name='PLD5', variant_id='chr1:242148188:A:G', interaction='None')
# plot_qtl(feature_name='PM20D1', variant_id='chr1:205813590:A:G', interaction='None')
# plot_qtl(feature_name='FAM198B-AS1', variant_id='chr4:158171481:A:G', interaction='None')     


# # age-interaction QTL
# qtls <- fread('QTL-output/signif-interaction-qtls.tsv')

# plot_qtl(feature_name='SMIM12', variant_id='chr1:34681413:C:T', interaction='Age')
# plot_qtl(feature_name='GABRB2', variant_id='chr5:161277005:T:C', interaction='Age')
# plot_qtl(feature_name='GABRB2', variant_id='chr5:161279237:C:G', interaction='Age')
# plot_qtl(feature_name='PTPRM', variant_id='chr18:7469023:C:G', interaction='Age')
# plot_qtl(feature_name='SLC7A11', variant_id='chr4:138217741:C:A', interaction='Age')
# plot_qtl(feature_name='PLD5', variant_id='chr1:242147130:T:G', interaction='Age')
# plot_qtl(feature_name='SLC18A1', variant_id='chr8:20175318:C:A', interaction='Age')

		
# plot_qtl(feature_name='SMIM12',  variant_id='chr1:34681413:C:T', interaction='Age')
# plot_qtl(feature_name='SMIM12',  variant_id='chr1:34681972:T:C', interaction='Age')
# plot_qtl(feature_name='AL158064.1', variant_id='chr13:80024792:G:T', interaction='Age')
# plot_qtl(feature_name='AL137781.1', variant_id='chr13:80024792:G:T', interaction='Age')
# plot_qtl(feature_name='SMIM12',  variant_id='chr1:34681413:C:T', interaction='Age')
# plot_qtl(feature_name='SMIM12',  variant_id='chr1:34681972:T:C', interaction='Age')
# plot_qtl(feature_name='AL158064.1', variant_id='chr13:80024792:G:T', interaction='Age')
# plot_qtl(feature_name='AL137781.1', variant_id='chr13:80024792:G:T', interaction='Age')
	
# plot_qtl(feature_name='PRDM5',variant_id='chr4:120622647:G:C', interaction='Age')
# plot_qtl(feature_name='SMIM12'   ,variant_id='chr1:34681972:T:C', interaction='Age')
# plot_qtl(feature_name='AL137781.1'  ,variant_id='chr13:80024792:G:T', interaction='Age')
# plot_qtl(feature_name='chr1_6601970_6603350'    ,variant_id='chr1:6589168:T:G', interaction='Age')
# plot_qtl(feature_name='chr3_69385114_69386888'   ,variant_id='chr3:69377088:G:A', interaction='Age')
# plot_qtl(feature_name='AL158064.1'  ,variant_id='chr13:80024792:G:T', interaction='Age')
# plot_qtl(eature_name='chr1_241639703_241640734'  ,variant_id='chr1:241728320:G:A', interaction='Age')
# plot_qtl(feature_name='chr19_58058668_58059926'  ,variant_id='chr19:58131851:G:A', interaction='Age')
# plot_qtl(feature_name='chr6_16712779_16713673'   ,variant_id='chr6:16782405:C:T', interaction='Age')
# plot_qtl(feature_name='chr20_461526_463296'    ,variant_id='chr20:402763:G:A', interaction='Age')
# plot_qtl(feature_name='KDM6B'   ,variant_id='chr17:7819345:G:A', interaction='Age')
# plot_qtl(feature_name='chr9_4489408_4491889'    ,variant_id='chr9:4407513:G:A', interaction='Age')
# plot_qtl(feature_name='KDM6B'   ,variant_id='chr17:7819999:G:C', interaction='Age')
# plot_qtl(feature_name='CPEB4'  ,variant_id='chr5:173826777:T:C', interaction='Age')
# plot_qtl(feature_name='CPEB4'  ,variant_id='chr5:173822320:A:G', interaction='Age')
# plot_qtl(feature_name='chr3_69385114_69386888'   ,variant_id='chr3:69379706:T:C', interaction='Age')
# plot_qtl(feature_name='chr3_69385114_69386888'   ,variant_id='chr3:69374758:T:C', interaction='Age')
# plot_qtl(feature_name='chr3_69385114_69386888'   ,variant_id='chr3:69375613:T:C', interaction='Age')
# plot_qtl(feature_name='chr17_20189006_20190001'  ,variant_id='chr17:20275827:T:C', interaction='Age')
# plot_qtl(feature_name='chr17_20189006_20190001'  ,variant_id='chr17:20177795:C:T', interaction='Age')
# plot_qtl(feature_name='chr17_20189006_20190001'  ,variant_id='chr17:20200473:A:G', interaction='Age')
# plot_qtl(feature_name='chr17_20189006_20190001'  ,variant_id='chr17:20250134:A:G', interaction='Age')
# plot_qtl(feature_name='chr17_20189006_20190001'  ,variant_id='chr17:20130902:G:A', interaction='Age')
# plot_qtl(feature_name='chr17_20189006_20190001'  ,variant_id='chr17:20133924:G:T', interaction='Age')
# plot_qtl(feature_name='chr17_20189006_20190001'  ,variant_id='chr17:20245036:T:C', interaction='Age')
# plot_qtl(feature_name='chr17_20189006_20190001'  ,variant_id='chr17:20271402:T:A', interaction='Age')
# plot_qtl(feature_name='chr19_58182018_58183985'  ,variant_id='chr19:58131851:G:A', interaction='Age')
# plot_qtl(feature_name='CPEB4'  ,variant_id='chr5:173826777:T:C', interaction='Age')
# plot_qtl(feature_name='hr12_121132313_121133477' ,variant_id='chr12:121138184:G:A', interaction='Age')
# plot_qtl(feature_name='chr17_20189006_20190001'  ,variant_id='chr17:20219862:G:A', interaction='Age')
# plot_qtl(feature_name='chr17_20189006_20190001'  ,variant_id='chr17:20265746:A:C', interaction='Age')
# plot_qtl(feature_name='chr6_16712779_16713673'   ,variant_id='chr6:16780117:G:A', interaction='Age')
# plot_qtl(feature_name='CPEB4'  ,variant_id='chr5:173828289:G:A', interaction='Age')
# plot_qtl(feature_name='ERC2'   ,variant_id='chr3:55551908:C:T', interaction='Age')
# plot_qtl(feature_name='ERC2'   ,variant_id='chr3:55553610:C:T', interaction='Age')
# plot_qtl(feature_name='ERC2'   ,variant_id='chr3:55574312:C:A', interaction='Age')
# plot_qtl(feature_name='ERC2'   ,variant_id='chr3:55555150:C:T', interaction='Age')
# plot_qtl(feature_name='ERC2'   ,variant_id='chr3:55551511:C:T', interaction='Age')
# plot_qtl(feature_name='ERC2'   ,variant_id='chr3:55553514:C:T', interaction='Age')
# plot_qtl(feature_name='ERC2'   ,variant_id='chr3:55572059:G:A', interaction='Age')
# plot_qtl(feature_name='SMIM12'   ,variant_id='chr1:34681413:C:T', interaction='Age')
# plot_qtl(feature_name='chr8_85462949_85464398'   ,variant_id='chr8:85494902:C:G', interaction='Age')
# plot_qtl(feature_name='chr19_58058668_58059926'  ,variant_id='chr19:58157968:C:T', interaction='Age')
# plot_qtl(feature_name='chr11_61947474_61947989'  ,variant_id='chr11:62032773:T:C', interaction='Age')
# plot_qtl(feature_name='chr17_32274901_32275735'  ,variant_id='chr17:32216259:A:G', interaction='Age')
# plot_qtl(feature_name='chr17_32274901_32275735'  ,variant_id='chr17:32206287:C:T', interaction='Age')