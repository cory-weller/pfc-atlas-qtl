#!/usr/bin/env bash

awk '$7 > 0.2 {print $1,$2}' genotypes/HBCC_polarized.pruned_pc.txt > genotypes/HBCC_to_remove.txt

module load plink/1.9
plink --bfile genotypes/HBCC_polarized \
      --keep-allele-order \
      --remove genotypes/HBCC_to_remove.txt \
      --make-bed \
      --out genotypes/HBCC_polarized_nooutliers