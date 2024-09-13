#!/usr/bin/env bash

module load plink/1.9

parallel -j 1 plink \
  --bfile data/genotypes/{1}-forQTL \
  --make-bed \
  --keep-allele-order \
  --output-chr chrM \
  --geno 0.05 \
  --hwe 0.000001 \
  --maf 0.01 \
  --keep tensorqtl-subsets/{1}-{2}-{3}.txt \
  --indiv-sort file tensorqtl-subsets/{1}-{2}-{3}.txt \
  --out tensorqtl-subsets/{1}-{2}-{3}-plink \
  ::: HBCC NABEC ::: rna atac ::: ExN InN MG VC OPC Oligo Astro