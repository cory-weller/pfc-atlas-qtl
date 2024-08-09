#!/usr/bin/env Rscript

library(data.table)
library(foreach)
library(ggplot2)
library(ggthemes)
library(viridis)

dat <- fread('fingerprints/crosscheck/all-pairs.txt')

dat[RIGHT_SAMPLE %like% 'HBCC', 'cohort' := 'HBCC']
dat[! RIGHT_SAMPLE %like% 'HBCC', 'cohort' := 'NABEC']
dat[, c('LOD_SCORE_NORMAL_TUMOR','DATA_TYPE','LOD_SCORE_TUMOR_NORMAL','LEFT_FILE','RIGHT_FILE','LEFT_GROUP_VALUE','RIGHT_GROUP_VALUE') := NULL]
all_samples <- fread('sample_bam.txt', header=F)$V1

dat <- dat[LEFT_SAMPLE %in% all_samples & RIGHT_SAMPLE %in% all_samples]

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

318 'clear cut' samples

clearcut.samples <- as.character(unique(dat[clearcut==TRUE, LEFT_SAMPLE]))
remaining.samples <- 
# remaining not-clear-cut include

# RIGHT_SAMPLE LEFT_SAMPLE              RESULT LOD_SCORE 
#       <fctr>      <fctr>              <char>     <num> 
#    HBCC_1832   HBCC_1832 UNEXPECTED_MISMATCH -3830.879 
#      KEN-845     KEN-845 UNEXPECTED_MISMATCH -3739.422 
#     KEN-1156    KEN-1156 UNEXPECTED_MISMATCH -3915.812 
#     SH-06-05    SH-06-05 UNEXPECTED_MISMATCH -3872.430 
#     SH-96-35    SH-96-35 UNEXPECTED_MISMATCH -3737.910 
#   UMARY-1544  UMARY-1544 UNEXPECTED_MISMATCH -4079.143 
#   UMARY-1818  UMARY-1818 UNEXPECTED_MISMATCH -3223.763 
#   UMARY-1845  UMARY-1845 UNEXPECTED_MISMATCH -4042.400 
#   UMARY-4789  UMARY-4789 UNEXPECTED_MISMATCH -3642.717 


#341 total samples in single cell set currently
all_samples

# 334 'tested' samples
tested.samples <- as.character(dat[, .N, by=LEFT_SAMPLE]$LEFT_SAMPLE)

# 7 untested samples that didn't have WGS data, presumably we keep
untested.samples <- base::setdiff(all_samples, tested.samples)

# "SH-94-35"   "UMARY-1158" "HBCC_1357"  "HBCC_2304"  "HBCC_1099" "HBCC_2003"  "HBCC_972"  

# of the 334 'tested' samples we can split into two groups:
# 318 'clear-cut good samples'
clearcut.samples

# and 16 'other' ambiguous samples
ambiguous.samples <- base::setdiff(tested.samples, clearcut.samples)
# HBCC_943 which has two RNAseq samples--itself, and mislabeled RNAseq HBCC_1832
# HBCC_1536 which is sample-swapped with HBCC_3081
# KEN-845 which doesn't match any genotype data, even itself
# KEN-1156 which doesn't match any genotype data, even itself
# SH-06-05 which is acutally duplicae SH-06-25
#  SH-96-35 which is actually SH-96-39
# UMARY-1544 which is actually duplicate for UMARY-544
# UMARY-1818 which doesn't match any genotype data, even itself
# UMARY-1845 which is actually duplicate UMARY-1843
# UMARY-4789 which doesn't match any genotype, even itself

# Giving us seven 'resolved' samples to add back
resolved.samples <- c('HBCC_943','HBCC_1536','HBCC_3081','SH-06-25','SH-96-39','UMARY-544','UMARY-1843')
# [1] "HBCC_943"   "HBCC_1536"  "HBCC_3081"  "SH-06-25"   "SH-96-39"  
# [6] "UMARY-544"  "UMARY-1843"

# 332 UMAP samples
UMAP.samples <- c(clearcut.samples, resolved.samples, untested.samples)
writeLines(sort(UMAP.samples), con='UMAP-samples.txt')

# 325 eQTL samples
eQTL.samples <- c(clearcut.samples, resolved.samples)

## RECOMMENDATIONS:
# - SWAP LABELS between HBCC_1536 and HBCC_3081
# - REMOVE current HBCC_943 and relabel HBCC_1832 as new replacement for HBCC_943
# - REMOVE KEN-845 (snRNA for KEN-845 has no match to any genotype, not even itself)
# - REMOVE KEN-1156 (snRNA for KEN-1156 has no match to any genotype, not even itself)
# - REMOVE SH-06-05 (is actually duplicate SH-06-25)
# - REMOVE SH-96-35 (is actually duplicate SH-96-39)
# - REMOVE UMARY-1544 (is actually duplicate UMARY-544)
# - REMOVE UMARY-1818 (snRNA for UMARY-1818 has no match to any genotype, not even itself)
# - REMOVE UMARY-1845 (is actually duplicate UMARY-1843)
# - REMOVE UMARY-4789 (snRNA for UMARY-4789 has no match to any genotype, not even itself)

# Leaving 332 UMAP samples (318 clear-cut matches of RNA x WGS; 7 untested samples with no WGS; and samples resolved by fingerprinting)
# Leaving 325 eQTL samples (above, less 7 samples without WGS)