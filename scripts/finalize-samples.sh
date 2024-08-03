#!/usr/bin/env bash

COHORT=${1}

bfile="data/genotypes/${COHORT}-forQTL"

cat $bfile.fam | \
    cut -d ' ' -f 1 | \
    sed "s/${COHORT}_//g" | \
    sed "s/$/-ARC/g" > ${COHORT}_samples_in_bfile.txt

wc -l ${COHORT}_samples_in_bfile.txt
# 135 samples in HBCC bfile
# 196 samples in NABEC bfile


module load R/4.3

Rscript - ${COHORT} <<EOF
library(data.table)
args <- commandArgs(trailingOnly=TRUE)
cohort <- args[1]

mycommand <- paste0('grep -f ', cohort, '_samples_in_bfile.txt ', cohort, '_bam_file_paths.txt')
dat <- fread(cmd=mycommand, header=F)
dat[, c('p1','p2','p3') := tstrsplit(V1, split='/')[6:8]]
dat[, 'p2' := NULL]
dat[, list(p1,.N), by=p3]
dat[, batch := tstrsplit(p1, split='batch')[2]]
dat[, batch := as.numeric(batch)]

dat.chosen <- dat[, .SD[batch==max(batch)], by=p3]
fwrite(dat.chosen[, .SD, .SDcols=c('p3','V1')], file=paste0(cohort, '.samples.txt'), quote=F, row.names=F, col.names=F, sep='\t')

EOF