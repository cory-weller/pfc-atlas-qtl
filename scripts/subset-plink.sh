#!/usr/bin/env bash

COHORT=$1
MODE=$2
CELLTYPE=$3
PROJDIR=$4

module load plink/1.9


plink \
  --bfile genotypes \
  --make-bed \
  --output-chr chrM \
  --keep-allele-order \
  --keep samples.txt \
  --indiv-sort file samples.txt \
  --out genotypes-forqtl
