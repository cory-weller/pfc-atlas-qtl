#!/usr/bin/env Rscript

## Libraries
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(data.table)

## Functions
plot_cov <- function(DT, pvar, .cohort, .cov) {
    dat <- DT[cohort==.cohort]
    dat[PC1 < -40 | abs(PC2) > 30, 'txt' := SampleID]
    dat[Ethnicity == 'Hispanic or Latino', 'txt' := SampleID]
    dat[cohort == 'NABEC', 'txt' := NA]
    pc1 <- round(pvar[2,1]*100, digits=4)
    pc2 <- round(pvar[2,2]*100, digits=4)
    
    g <- ggplot() +
        geom_point(data=dat[Ethnicity != 'Hispanic or Latino'], aes(x=PC1, y=PC2, color=.data[[.cov]])) +
        geom_point(data=dat[Ethnicity == 'Hispanic or Latino'], aes(x=PC1, y=PC2, color=.data[[.cov]])) +
        geom_text_repel(data=dat, aes(label=txt, x=PC1, y=PC2), size =1, segment.size=0.2, color='black', max.overlaps=20, min.segment.length = unit(0, 'lines'), ) +
        theme_few() +
        labs(x=paste0('PC1: ', pc1, '% var'), y=paste0('PC2: ', pc2, '% var'))
        
    if(! is.factor(dat[[.cov]])) {
        g <- g + scale_color_viridis()
    }
    
    # Save SVGs
    outfile_stem <- paste0('genotypes/PCA-plots/', toupper(.cohort), '-', .cov)
    ggsave(g, file=paste0(outfile_stem, '.svg'), width=12, heigh=12, units='cm')
    ggsave(g, file=paste0(outfile_stem, '.png'), width=12, heigh=12, units='cm')
}


encode_var <- function(DT, variable, prefix=NULL) {
    dt <- copy(DT)
    #new_var_name <- paste0(prefix, variable)
    biggest_var <- dt[, .N, by=variable][order(-N)][1, get(variable)]
    biggest_var <- paste0(prefix, biggest_var)
    old_varnames <- unique(dt[,get(variable)])
    new_varnames <- paste0(prefix, old_varnames)

    dt[, (new_varnames) := lapply(old_varnames, function(x) as.numeric(get(variable) == x))]
    dt[, (biggest_var) := NULL]     # Remove most abundant group
    dt[, (variable) := NULL]        # remove original categorical column
    return(dt[])
}

encode_sex <- function(DT, sex_col_name) {
    dt <- copy(DT)
    
    dt[get(sex_col_name) == 'Female', 'fsex' := 0]
    dt[get(sex_col_name) == 'Male', 'fsex' := 1]
    
    dt[, (sex_col_name) := NULL]     # Remove original col
    return(dt[])
}


fillNA <- function(DT, x) {
    # replaces all NA values with x in data.table DT
    for (j in seq_len(ncol(DT)))
        set(DT,which(is.na(DT[[j]])),j,x)
}



# dat <- fread('HBCC_pruned_for_pcs.traw')
# setnames(dat, gsub('^.*_', '', colnames(dat)))
# dat[, c(1,2,3,4,5,6) := NULL]


# fillNA(dat, 0)

# mat <- t(as.matrix(dat))
# res <- prcomp(mat, center = TRUE, scale = TRUE, rank.=20)
# p_var <- summary(res)$importance
# pcs <- as.data.table(res$x, keep.rownames=T)
# setnames(pcs, 'rn', 'SampleID')
# hbcc.pcs <- copy(pcs)


# covs <- fread('../pfc-metadata-357.csv')
# hbcc.covs <- covs[cohort == 'HBCC']


# pcs <- merge(hbcc.covs, hbcc.pcs, by.x='SampleID', by.y='SampleID')
# covs_to_plot <- c('Ethnicity','Race','PMI','Sex','Age','Brain_bank','Ancestry','cohort','batch')
# pcs[, Ethnicity := factor(Ethnicity)]
# pcs[, Race := factor(Race)]
# pcs[, Sex := factor(Sex)]
# pcs[, Brain_bank := factor(Brain_bank)]
# pcs[, Ancestry := factor(Ancestry)]
# pcs[, cohort := factor(cohort)]
# pcs[, batch := factor(batch)]





# ggplot(data=pcs, aes(x=PC1, y=PC2)) + geom_point() + coord_fixed(ratio = 1) + theme_few()


# plot_cov(pcs, p_var, 'HBCC', 'Ethnicity')

# for(i in covs_to_plot) {
#     for(cohort in c('HBCC')) {
#         plot_cov(pcs, cohort, i)
#     }
# }



# Get HBCC PCs
dat <- fread('HBCC_pruned_for_pcs.traw')
setnames(dat, gsub('^.*_', '', colnames(dat)))
dat[, c(1,2,3,4,5,6) := NULL]
fillNA(dat, 0)
mat <- t(as.matrix(dat))
res <- prcomp(mat, center = TRUE, scale = TRUE, rank.=20)
pc1pct <- summary(res)$importance[2,1]
pc2pct <- summary(res)$importance[2,2]
pcs <- as.data.table(res$x, keep.rownames=T)
setnames(pcs, 'rn', 'SampleID')
hbcc.pcs <- copy(pcs)


# Get NABEC PCs
dat <- fread('NABEC_pruned_for_pcs.traw')
setnames(dat, gsub('^.*_', '', colnames(dat)))
dat[, c(1,2,3,4,5,6) := NULL]
fillNA(dat, 0)
mat <- t(as.matrix(dat))
res <- prcomp(mat, center = TRUE, scale = TRUE, rank.=20)
pc1pct <- summary(res)$importance[2,1]
pc2pct <- summary(res)$importance[2,2]
pcs <- as.data.table(res$x, keep.rownames=T)
setnames(pcs, 'rn', 'SampleID')
nabec.pcs <- copy(pcs)


covs <- fread('../pfc-metadata-357.csv')

dat <- merge(covs, rbind(nabec.pcs, hbcc.pcs), by='SampleID', all=T)

to_exclude <- c('HBCC-1462','HBCC-3029','HBCC-1456','HBCC-1105','HBCC-2267','HBCC-1385','HBCC-2429','HBCC-1058','HBCC-2756','HBCC-1331','HBCC-1560','HBCC-1431','HBCC-2781')


# one-hot encode together


# for performing one-hot encoding separately
hbcc <- dat[cohort=='HBCC'][! SampleID %in% to_exclude][,-c('cohort')]
nabec <- dat[cohort=='NABEC'][! SampleID %in% to_exclude][,-c('cohort')]

# one-hot encode batch
hbcc <- encode_var(hbcc, 'batch', prefix='batch_')
nabec <- encode_var(nabec, 'batch', prefix='batch_')
dat <- encode_var(dat, 'batch', prefix='batch_')

# one-hot encode brain bank
hbcc <- encode_var(hbcc, 'Brain_bank', prefix='bbank_')
nabec <- encode_var(nabec, 'Brain_bank', prefix='bbank_')
dat <- encode_var(dat, 'Brain_bank', prefix='bbank_')

# Encode sex as binary
hbcc <- encode_sex(hbcc, 'Sex')
nabec <- encode_sex(nabec, 'Sex')
dat <- encode_sex(dat, 'Sex')

# Exclude columns
hbcc[, c('Ethnicity','Race','Ancestry') := NULL]
nabec[, c('Ethnicity','Race','Ancestry') := NULL]
dat[, c('Ethnicity','Race','Ancestry') := NULL]

# Export covariates
# fwrite(nabec, file='NABEC_tqtl_covariates.csv', quote=F, row.names=F, col.names=T, sep=',')
# fwrite(hbcc, file='HBCC_tqtl_covariates.csv', quote=F, row.names=F, col.names=T, sep=',')

# quit()

dat_covs <- as.data.frame(t(dat[,-c('SampleID')]))
colnames(dat_covs) <- dat$SampleID
rownames(dat_covs) <- colnames(dat[,-c(1)])
dat_covs <- dat_covs[rownames(dat_covs)[c(1,2,24:36)],]


hbcc_covs <- as.data.frame(t(hbcc[,-c('SampleID')]))
colnames(hbcc_covs) <- hbcc$SampleID
rownames(hbcc_covs) <- colnames(hbcc[,-c(1)])



nabec_covs <- as.data.frame(t(nabec[,-c('SampleID')]))
colnames(nabec_covs) <- nabec$SampleID
rownames(nabec_covs) <- colnames(nabec[,-c(1)])


fwrite(dat_covs, file='COMBINED_one-hot_covariates.csv', quote=F, row.names=T, col.names=T, sep=',')
fwrite(nabec_covs, file='NABEC_tqtl_covariates.csv', quote=F, row.names=T, col.names=T, sep=',')
fwrite(hbcc_covs, file='HBCC_tqtl_covariates.csv', quote=F, row.names=T, col.names=T, sep=',')


