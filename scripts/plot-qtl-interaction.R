#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(ggthemes)
library(viridis)
library(ggrepel)
library(foreach)

# Get non-interaction QTLs
interaction_qtl_filename <- 'significant-interaction-qtls.tsv.gz'
noninteraction_qtl_filename <- 'significant-noninteraction-qtls.tsv.gz'
pbulk_method <- 'sum'

# Load and format data
qtls.noninteraction <- fread(noninteraction_qtl_filename)
qtls.noninteraction <- qtls.noninteraction[pseudobulk_method == pbulk_method]
qtls.interaction <- fread(interaction_qtl_filename)
qtls.interaction <- qtls.interaction[pseudobulk_method == pbulk_method]
qtls.full <- list(qtls.noninteraction, qtls.interaction)

qtls.noninteraction <- qtls.noninteraction[, .N, by=list(celltype, cohort, mode, interaction)]
qtls.interaction <- qtls.interaction[, .N, by=list(celltype, cohort, mode, interaction)]
qtls <- rbindlist(list(qtls.noninteraction, qtls.interaction))


qtls[celltype == 'ExN', celltype := 'Excitatory Neuron']
qtls[celltype == 'InN', celltype := 'Inhibitory Neuron']
qtls[celltype == 'MG', celltype := 'Microglia']
qtls[celltype == 'Oligo', celltype := 'Oligodendrocyte']
qtls[celltype == 'OPC', celltype := 'OPC']
qtls[celltype == 'VC', celltype := 'Vascular Cell']
qtls[celltype == 'Astro', celltype := 'Astrocyte']
qtls[, celltype := factor(celltype, levels=c('Excitatory Neuron','Inhibitory Neuron','Astrocyte','Microglia','Oligodendrocyte','OPC','Vascular Cell'))]

qtls[, cohort_total := sum(N), by=list(cohort, mode, interaction)]
qtls[, Proportion := N/cohort_total]


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
g1 <- plot_nqtls(qtls, plt.colors)
ggsave(g1, file='QTL-plots/N-QTLs.png', width=20, height=20, units='cm')
ggsave(g1, file='QTL-plots/N-QTLs.svg', width=20, height=20, units='cm')

plot_propqtls <- function(DT, clrs) {
    ggplot(DT, aes(x=cohort, y=Proportion, fill=celltype)) +
    geom_bar(stat='identity') +
    theme_few(14) +
    labs(fill='Cell Type', x='Cohort', y='Proportion', title='Significant QTLs (BH correction)') +
    facet_grid(mode~interaction, scales='free_y', labeller=labeller(.rows=label_both,.cols=label_both)) +
    guides(fill = guide_legend()) +
    scale_fill_manual(values=clrs)
}

g2 <- plot_propqtls(qtls, plt.colors)
ggsave(g2, file='QTL-plots/Proportion-QTLs.png', width=20, height=20, units='cm')
ggsave(g2, file='QTL-plots/Proportion-QTLs.svg', width=20, height=20, units='cm')


# # Top 100 gene hits for Genotype x Age
# top_hits <- qtls[Interaction==TRUE][, .N, by=phenotype_id][order(-N)][1:100]$phenotype_id
# hit_counts <- qtls[phenotype_id %chin% top_hits & Interaction == TRUE][, .N, by=list(cohort, phenotype_id)]
# hit_counts[, phenotype_id  := factor(phenotype_id, levels=top_hits)]
# g3 <- ggplot(hit_counts,
#         aes(x=phenotype_id, y=N, fill=cohort, group=phenotype_id)) + geom_bar(stat='identity') +
#         theme_few() +
#         theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
#         labs(x='Gene', y='N hits', fill='Cohort', title='Top 100 Genes by number of  Interaction (Genotype x Age) iQTL hits across both cohorts')

# # Top 100 gene hits for Genotype alone
# top_hits <- qtls[GenoAssoc==TRUE][, .N, by=phenotype_id][order(-N)][1:100]$phenotype_id
# hit_counts <- qtls[phenotype_id %chin% top_hits & GenoAssoc == TRUE][, .N, by=list(cohort, phenotype_id)]
# hit_counts[, phenotype_id  := factor(phenotype_id, levels=top_hits)]
# g4 <- ggplot(hit_counts,
#         aes(x=phenotype_id, y=N, fill=cohort, group=phenotype_id)) + geom_bar(stat='identity') +
#         theme_few() +
#         theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
#         labs(x='Gene', y='N hits', fill='Cohort', title='Top 100 Genes by number of eQTL hits across both cohorts')

# g4.5 <- ggplot(hit_counts[phenotype_id != 'ARL17B'],
#         aes(x=phenotype_id, y=N, fill=cohort, group=phenotype_id)) + geom_bar(stat='identity') +
#         theme_few() +
#         theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
#         labs(x='Gene', y='N hits', fill='Cohort', title='Top 2-100 Genes by number of eQTL hits across both cohorts\nARL17B Excluded')

# # Top 100 gene hits for Age alone
# top_hits <- qtls[AgeAssoc==TRUE][, .N, by=phenotype_id][order(-N)][1:100]$phenotype_id
# hit_counts <- qtls[phenotype_id %chin% top_hits & AgeAssoc == TRUE][, .N, by=list(cohort, phenotype_id)]
# hit_counts[, phenotype_id  := factor(phenotype_id, levels=top_hits)]
# g5 <- ggplot(hit_counts,
#         aes(x=phenotype_id, y=N, fill=cohort, group=phenotype_id)) + geom_bar(stat='identity') +
#         theme_few() +
#         theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
#         labs(x='Gene', y='N hits', fill='Cohort', title='Top 100 Age-associated Genes (by number of Age-associated hits) across both cohorts')

# ggsave(g3, file='QTL-plots/iQTL-barplots.png', width=45, height=12, units='cm')
# ggsave(g3, file='QTL-plots/iQTL-barplots.svg', width=45, height=12, units='cm')

# ggsave(g4, file='QTL-plots/eQTL-barplots.png', width=45, height=12, units='cm')
# ggsave(g4, file='QTL-plots/eQTL-barplots.svg', width=45, height=12, units='cm')
# ggsave(g4.5, file='QTL-plots/eQTL-barplots-no-ARL17B.png', width=45, height=12, units='cm')
# ggsave(g4.5, file='QTL-plots/eQTL-barplots-no-ARL17B.svg', width=45, height=12, units='cm')


# ggsave(g5, file='QTL-plots/Age-QTL-barplots.png', width=45, height=12, units='cm')
# ggsave(g5, file='QTL-plots/Age-QTL-barplots.svg', width=45, height=12, units='cm')


# Plot examples of Age-only Associated genes


# Plot examples of Genotype-only Associated genes (pure eQTL)

# Plot example of iQTL













# plot_celltype_qtl_ranks <- function(DT, .cohort, .celltype, signif_type) {
#     dat <- copy(DT[which(qtls[[signif_type]])][celltype==.celltype][cohort==.cohort])
#     top_hits <- dat[, .N, by=phenotype_id][order(-N)][1:100]$phenotype_id
#     top_hits <- as.character(na.omit(top_hits))
#     hit_counts <- dat[phenotype_id %chin% top_hits, .N, by=phenotype_id][order(-N)]
#     hit_counts[, phenotype_id  := factor(phenotype_id, levels=top_hits)]
#     g <- ggplot(hit_counts,
#         aes(x=phenotype_id, y=N, group=phenotype_id)) + geom_bar(stat='identity') +
#         theme_few() +
#         theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
#         labs(x='Gene', y='N hits', title='Top 100 Genes by number of associated genes across both cohorts')
#     return(g)
# }

if(FALSE) {
# d1 <- plot_celltype_qtl_ranks(qtls, 'HBCC', 'Astrocyte', 'AgeAssoc')
AD_GWAS_fn <- 'data/AD_Bellenguez_GrCh38_full.csv.gz'
AD_GWAS <- fread(AD_GWAS_fn)
AD_GWAS <- AD_GWAS[p < 5e-8]    # 5637 significant hits

setnames(AD_GWAS, c('A1','A2','p','chr','bp'), c('REF','ALT','P','CHROM','POS'))
AD_GWAS <- AD_GWAS[, .SD, .SDcols=c('CHROM','POS','REF','ALT','P')]
AD_GWAS[, CHROM := paste0('chr',CHROM)]
setkey(AD_GWAS, CHROM, POS)

qtls[, c('CHROM','POS','REF','ALT') := tstrsplit(variant_id, split=':')]
qtls[, POS := as.numeric(POS)]
setkey(qtls, CHROM, POS)

dat.merge <- merge(qtls, AD_GWAS, all.x=T)
dat.merge[, c('REF.x','REF.y','ALT.x','ALT.y') := NULL]
dat.merge[!is.na(P), AD_GWAS_SIGNIF := TRUE]
dat.merge[is.na(P), AD_GWAS_SIGNIF := FALSE]
dat.merge[, P := NULL]


fwrite(dat.merge[AD_GWAS_SIGNIF==TRUE], file='AD_GWAS_x_QTL.tsv', sep='\t')


dat.merge <- dat.merge[!is.na(POS)]
dat.merge[, id := paste(CHROM, POS, sep='_')]

dat.merge[GenoAssoc==TRUE][!duplicated(id)] # 44059 variants
dat.merge[GenoAssoc==TRUE & AD_GWAS_SIGNIF==TRUE][!duplicated(id)] # 1473 variants
dat.merge[GenoAssoc==TRUE & AD_GWAS_SIGNIF==TRUE][!duplicated(id)][, .N, by=phenotype_id][order(-N)]


dat.merge[Interaction==TRUE][!duplicated(id)] # 58650 variants
dat.merge[Interaction==TRUE & AD_GWAS_SIGNIF==TRUE][!duplicated(id)] # 17 variants

dat.merge[Interaction==TRUE & AD_GWAS_SIGNIF==TRUE][!duplicated(id)][, .N, by=phenotype_id][order(-N)]

}
# STMN4
# APP
# MAPT
# KANSL1
# PVR
# MS4A7

get_loci <- function(.DT, .gene) {
    o <- .DT[phenotype_id==.gene][GenoAssoc==TRUE]
    return(o$variant_id)
}

# Read in dosage genotype table
genos <- fread('data/alt-allele-dosage.tsv', select=c('ID','Sample','Dosage'))
# Subset loci that we care about in QTL table
qtls <- qtls.full[[1]] # RNA
genos <- genos[ID %chin% qtls$variant_id]
genos[Sample %like% '^HBCC', Cohort := 'HBCC']
genos[Sample %like% '^[USK]', Cohort := 'NABEC']

gc()

get_genotypes <- function(.DT, .loci) {
    return(.DT[ID %chin% .loci])
}

melt_expression <- function(DT) {
    info_cols <- c('#chr','start','end','phenotype_id')
    sample_cols <- colnames(DT)[2:length(colnames(DT))]
    out <- melt(DT, measure.vars=sample_cols, variable.name='Sample', value.name='Expression')
    return(out)
}

get_expression <- function(.genes, .modes=c('rna'), .cohorts=c('HBCC','NABEC'), .celltypes=c('ExN','InN','Astro','OPC','Oligo','VC','MG')) {
    foreach(celltype=.celltypes, .combine='rbind') %do% {
        foreach(cohort=.cohorts, .combine='rbind') %do% {
            foreach(mode=.modes, .combine='rbind') %do% {
                cat(celltype, cohort, mode)
                dt.tmp <- fread(paste0('data/', mode, '/', cohort, '-', celltype, '-counts.bed'), header=TRUE)
                o <- dt.tmp[phenotype_id %chin% .genes]
                o[, '#chr' := NULL]
                o[, 'start' := NULL]
                o[, 'end' := NULL]
                o <- unique(o)
                o <- melt_expression(o)
                o[, 'cohort' := cohort]
                o[, 'mode' := mode]
                o[, 'celltype' := celltype]
                return(o[])
            }
        }
    }
}

genes_of_interest <- c('STMN4','APP','MAPT','KANSL1','PVR','MS4A7','DGKB','PLD5','RPLP0','FAS')

DT.expression <- get_expression(genes_of_interest)
setnames(DT.expression, 'phenotype_id', 'Gene')

DT.genos <- foreach(gene = genes_of_interest, .combine='rbind') %do% {
    goi_genotypes <- get_genotypes(genos, get_loci(qtls, gene))
    goi_genotypes[, 'Gene' := gene]
    return(goi_genotypes[])
}

DT.genos2 <- foreach(celltype=c('ExN','InN','Astro','OPC','Oligo','VC','MG'), .combine='rbind') %do% {
    dat <- copy(DT.genos)
    dat[, 'celltype' := celltype]
    return(dat[])
}

setkey(DT.expression, Sample, Gene, celltype)
setkey(DT.genos2, Sample, Gene, celltype)

dat <- merge(DT.genos2, DT.expression)[!is.na(Dosage)]
dat[, Dosage := factor(Dosage, levels=c(0,1,2))]

plot_qtl_gene <- function(DT, .gene, .variant) {
    dat <- copy(DT[Gene == .gene & ID == .variant])
    g <- ggplot(dat, aes(x=Dosage, y=Expression)) +
        facet_grid(celltype ~ cohort) +
        geom_jitter(shape=21, alpha=0.5) +
        geom_boxplot(outliers=FALSE, fill='white', alpha=0.6) +
        labs(title=.gene, x=paste0(.variant, ' alt allele dosage')) +
        theme_few()
    return(g)
}

# GOOD
plot_qtl_gene(dat, 'STMN4', 'chr8:27247210:G:A')

# Function to read in specific variant using TABIX; output expression boxplot




# TODO: Overlap with African-American AD GWAS sig hits

plot_qtl_gene(dat, 'SNX19', 'chr11:130894514:T:C')

plot_qtl_gene(dat, 'APP', 'chr21:25152142:A:G')

plot_qtl_gene(dat, 'APP', 'chr21:25152142:A:G')
plot_qtl_gene(dat, 'APP', 'chr21:25152502:C:T')
plot_qtl_gene(dat, 'APP', 'chr21:25788268:A:G')
plot_qtl_gene(dat, 'APP', 'chr21:26068891:G:T')
plot_qtl_gene(dat, 'APP', 'chr21:25148434:G:A')
plot_qtl_gene(dat, 'APP', 'chr21:25152142:A:G')
plot_qtl_gene(dat, 'APP', 'chr21:25152502:C:T')
plot_qtl_gene(dat, 'APP', 'chr21:25788268:A:G')
plot_qtl_gene(dat, 'APP', 'chr21:25800873:AAAAC:A')
plot_qtl_gene(dat, 'APP', 'chr21:26129491:T:C')
plot_qtl_gene(dat, 'APP', 'chr21:26146434:C:T')
g.app <- plot_qtl_gene(dat, 'APP', 'chr21:26161943:T:C')
ggsave(g.app, file='QTL-plots/APP.png', height=20, width=10, units='cm')

plot_qtl_gene(dat, 'APP', 'chr21:26661974:G:T')
plot_qtl_gene(dat, 'APP', 'chr21:26666452:CCTAT:C')
plot_qtl_gene(dat, 'APP', 'chr21:26681660:C:G')
plot_qtl_gene(dat, 'APP', 'chr21:26687517:A:T')
plot_qtl_gene(dat, 'APP', 'chr21:26689680:A:G') # Maybe?
plot_qtl_gene(dat, 'APP', 'chr21:26692058:G:A')
plot_qtl_gene(dat, 'APP', 'chr21:26692678:A:G') # Maybe?
plot_qtl_gene(dat, 'APP', 'chr21:26692797:C:T')
plot_qtl_gene(dat, 'APP', 'chr21:26695561:C:T')
plot_qtl_gene(dat, 'APP', 'chr21:26695897:T:C')
plot_qtl_gene(dat, 'APP', 'chr21:26696572:G:A')
plot_qtl_gene(dat, 'APP', 'chr21:26698723:A:T')
plot_qtl_gene(dat, 'APP', 'chr21:26704856:A:G')
plot_qtl_gene(dat, 'APP', 'chr21:26707602:G:A')
plot_qtl_gene(dat, 'APP', 'chr21:26711203:T:C')
plot_qtl_gene(dat, 'APP', 'chr21:26711418:G:A')
plot_qtl_gene(dat, 'APP', 'chr21:26822209:C:T')
plot_qtl_gene(dat, 'APP', 'chr21:26829556:A:C')
plot_qtl_gene(dat, 'APP', 'chr21:26829769:T:A')
plot_qtl_gene(dat, 'APP', 'chr21:26830789:G:T')
plot_qtl_gene(dat, 'APP', 'chr21:26830811:G:T')
plot_qtl_gene(dat, 'APP', 'chr21:26831332:T:C')
plot_qtl_gene(dat, 'APP', 'chr21:26833321:C:T')
plot_qtl_gene(dat, 'APP', 'chr21:26833470:G:A')
plot_qtl_gene(dat, 'APP', 'chr21:26682828:C:T')
plot_qtl_gene(dat, 'APP', 'chr21:26690468:C:T')


g.pvr <- plot_qtl_gene(dat, 'PVR', 'chr19:44892362:A:G') # Yes
ggsave(g.pvr, file='QTL-plots/PVR.png', height=20, width=10, units='cm')

g.pvr2 <- plot_qtl_gene(dat, 'PVR', 'chr19:45222190:G:C')
ggsave(g.pvr2, file='QTL-plots/PVR-2.png', height=20, width=10, units='cm')



plot_qtl_gene(dat, 'PVR', 'chr19:45221485:G:A') # Yes
g.stmn4_1 <- plot_qtl_gene(dat, 'STMN4', 'chr8:26630730:G:T') # Maybe
g.stmn4_2 <- plot_qtl_gene(dat, 'STMN4', 'chr8:27226026:G:A') # Maybe
g.stmn4 <- plot_qtl_gene(dat, 'STMN4', 'chr8:27572213:G:A') # Maybe
ggsave(g.stmn4, file='QTL-plots/STMN4.png', height=20, width=10, units='cm')


qtls[phenotype_id == 'STMN4'][Interaction==TRUE]

dat.merge[GenoAssoc==TRUE, .SD, .SDcols=c('CHROM','POS')]


AD_GWAS_signif <- nrow(AD_GWAS)

dat.merge[Interaction==TRUE & AD_GWAS_SIGNIF==TRUE , .N, by=list(cohort, celltype)][order(-N)]
dat.merge[GenoAssoc==TRUE & AD_GWAS_SIGNIF==TRUE , .N, by=list(cohort, celltype)][order(-N)]


# N AD_GWAS significant iQTL

# N AD_GWAS significant eQTL

dat.merge[Interaction==TRUE & AD_GWAS_SIGNIF==TRUE, .N, by=list( celltype, cohort, AD_GWAS_SIGNIF)][order(-N)]
dat.merge[AgeAssoc==TRUE & AD_GWAS_SIGNIF==TRUE, .N, by=list( celltype, cohort, AD_GWAS_SIGNIF)][order(-N)]
dat.merge[GenoAssoc==TRUE & AD_GWAS_SIGNIF==TRUE, .N, by=list( celltype, cohort, AD_GWAS_SIGNIF)][order(-N)]

rs_id_unique <- unique(rs_ids$SNP)


# Using dbsnp vcf
rs.vcf <- fread('rsids.txt', header=F)
rs.vcf <- rs.vcf[V1 %in% 1:22]
rs.vcf <- rs.vcf[V3 %chin% rs_id_unique]

setnames(rs.vcf, c('CHROM','POS','ID','REF','ALT'))

# Get eQTLs to plot


library(biomaRt)
library(data.table)



snp_mart <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")


# Query BioMart
results <- getBM(
  attributes = c("refsnp_id", "chr_name", "chrom_start"),
  filters = "snp_filter",
  values = rs_ids$SNP,
  mart = snp_mart
)

plt.colors <- data.table(
    'Oligo' = '#945247ff',
    'ExN' = '#ff7500ff',
    'InN' = '#00a162ff'
    'Astro' = '#0078b8ff',
    'MG' = '#ea0016ff',
    'OPC' = '#b735ffff',
    'VR' = '#f26ec5ff'    
    )
