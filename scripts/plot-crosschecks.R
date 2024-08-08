#!/usr/bin/env Rscript

library(data.table)
library(foreach)
library(ggplot2)
library(ggthemes)
library(viridis)

dat <- fread('fingerprints/all-crosscheck.tsv')

dat[RIGHT_SAMPLE %like% 'HBCC', 'cohort' := 'HBCC']
dat[! RIGHT_SAMPLE %like% 'HBCC', 'cohort' := 'NABEC']
dat[, c('LOD_SCORE_NORMAL_TUMOR','DATA_TYPE','LOD_SCORE_TUMOR_NORMAL','LEFT_FILE','RIGHT_FILE','LEFT_GROUP_VALUE','RIGHT_GROUP_VALUE') := NULL]

samples <- data.table(id=unique(c(dat$LEFT_SAMPLE, dat$RIGHT_SAMPLE)))
samples[id %like% 'UMARY|KEN', c('BrainBank','ID_N') := tstrsplit(id, split='-')]
samples[id %like% 'SH', c('BrainBank','ID_N1','ID_N2') := tstrsplit(id, split='-')]
samples[id %like% 'SH', 'ID_N' := paste0(ID_N1, ID_N2)]
samples[, c('ID_N1','ID_N2') := NULL]
samples[id %like% '^HBCC', BrainBank := 'HBCC']
samples[BrainBank == 'HBCC', 'ID_N' := tstrsplit(id, split='_')[2]]
samples[, ID_N := as.numeric(ID_N)]
setkey(samples, BrainBank, ID_N)

cols_to_remove <- c('_RUN_BARCODE','_LANE','_MOLECULAR_BARCODE_SEQUENCE','_LIBRARY')
dat[, paste0('LEFT', cols_to_remove) := NULL]
dat[, paste0('RIGHT', cols_to_remove) := NULL]

dat[, LEFT_SAMPLE := factor(LEFT_SAMPLE, levels=samples$id)]
dat[, RIGHT_SAMPLE := factor(RIGHT_SAMPLE, levels=samples$id)]


dat[LEFT_SAMPLE==RIGHT_SAMPLE, self := TRUE]
dat[LEFT_SAMPLE!=RIGHT_SAMPLE, self := FALSE]

dat <- merge(dat, dat[LOD_SCORE > 0, .N, by=LEFT_SAMPLE], by='LEFT_SAMPLE', all=T)
setnames(dat, 'N', 'N_LEFT')
dat[is.na(N_LEFT), N_LEFT := 0]
dat <- merge(dat, dat[LOD_SCORE > 0, .N, by=RIGHT_SAMPLE], by='RIGHT_SAMPLE', all=T)
setnames(dat, 'N', 'N_RIGHT')
dat[is.na(N_RIGHT), N_RIGHT := 0]

# Left is genotype
# HBCC_823   HBCC_1058  HBCC_1134  HBCC_1228  HBCC_1331  HBCC_1385 
# HBCC_1431  HBCC_1560  HBCC_1832  HBCC_2429  HBCC_2756  HBCC_2781 
# HBCC_2998  KEN-845    KEN-1156   SH-06-05   SH-96-35   UMARY-1064
# 339 samples with genotype data
# 326 samples with RNA data

# Right is RNA
# KEN-845    KEN-1156   UMARY-1818 UMARY-4789

clearcut <- as.character(dat[N_LEFT==1 & N_RIGHT==1 & LEFT_SAMPLE==RIGHT_SAMPLE & LOD_SCORE > 0]$LEFT_SAMPLE)

dat[LEFT_SAMPLE %in% clearcut & RIGHT_SAMPLE %in% clearcut, 'clearcut' := TRUE]
dat[is.na(clearcut), clearcut := FALSE]

# 339 samples
# 339 with genotypes
# 326 with rnaseq
# 317 total samples with genotypes and snRNAseq
# 310 reciprocal best match, clear-cut good samples
ambig_samples <- unique(c(dat[clearcut==FALSE & LOD_SCORE > 0]$LEFT_SAMPLE, dat[clearcut==FALSE & LOD_SCORE > 0]$RIGHT_SAMPLE))


plotSubset <- function(DT, fn) {
    n_x <- length(unique(DT[, LEFT_SAMPLE]))
    n_y <- length(unique(DT[, RIGHT_SAMPLE]))

    g1 <- ggplot(DT, aes(x=LEFT_SAMPLE, y=RIGHT_SAMPLE, fill=LOD_SCORE)) +
        geom_tile() +
        theme_few() +
        scale_x_discrete(guide = guide_axis(angle = 90)) +
        scale_fill_viridis() +
        labs(x='Sample Genotype', y='Sample snRNA')

    ggsave(g1, file=fn, width=10+n_x/3, height=10+n_y/3, units='cm', limitsize=FALSE)

}

plotSubset(dat, fn='plots/fingerprint-all.png')
plotSubset(dat[clearcut==TRUE], fn='plots/fingerprint-clearcut.png')
plotSubset(dat[LEFT_SAMPLE %in% ambig_samples & RIGHT_SAMPLE %in% ambig_samples], fn='plots/fingerprint-ambiguous.png')


dat[LEFT_SAMPLE %chin% ambig_samples | RIGHT_SAMPLE %chin% ambig_samples]

plotSubset(dat[LEFT_SAMPLE %in% ambig_samples & RIGHT_SAMPLE %in% ambig_samples])




# LEFT_SAMPLE matches 1 thing, and it's itself:
left_selfmatch <- dat[LEFT_SAMPLE %chin% dat[LOD_SCORE > 0, .N, by=list(LEFT_SAMPLE)][N==1, LEFT_SAMPLE]][LEFT_SAMPLE==RIGHT_SAMPLE, LEFT_SAMPLE]
right_selfmatch <- dat[RIGHT_SAMPLE %chin% dat[LOD_SCORE > 0, .N, by=list(RIGHT_SAMPLE)][N==1, RIGHT_SAMPLE]][LEFT_SAMPLE==RIGHT_SAMPLE, RIGHT_SAMPLE]
reciprocal_bestmatch <- intersect(left_selfmatch, right_selfmatch)

dat[LEFT_SAMPLE %chin% reciprocal_bestmatch & RIGHT_SAMPLE %chin% reciprocal_bestmatch, status := 'SELFMATCH']
dat[! (LEFT_SAMPLE %chin% reciprocal_bestmatch & RIGHT_SAMPLE %chin% reciprocal_bestmatch), status := 'AMBIGUOUS']


plotSubset <- function(DT) {
    n_x <- length(unique(DT[, LEFT_SAMPLE]))
    n_y <- length(unique(DT[, RIGHT_SAMPLE]))

    ggplot(DT, aes(x=LEFT_SAMPLE, y=RIGHT_SAMPLE, fill=LOD_SCORE)) +
        geom_tile() +
        theme_few() +
        scale_x_discrete(guide = guide_axis(angle = 90)) +
        scale_fill_viridis() +
        labs(x='Sample Genotype', y='Sample snRNA')
}

plotSubset(dat[status=='SELFMATCH'])

fwrite(dat, file=paste0('~/pfc-atlas-qtl/fingerprints/all-crosscheck.tsv'), quote=F, row.names=F, col.names=T, sep='\t')



ggsave(g, file=paste0('~/pfc-atlas-qtl/fingerprints/crosscheck.png'), width=10+n_x/3, height=n_y/3, units='cm', limitsize=FALSE)

