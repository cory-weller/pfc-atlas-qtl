#!/usr/bin/env bash
#SBATCH --partition quick
#SBATCH --time 3:30:00
#SBATCH --mem 12G

dirs=$(find ${PWD} -mindepth 1 -maxdepth 1 -type d -name '*atac')

for dir in ${dirs}; do
    cd $dir
    python3 /data/CARD_singlecell/users/wellerca/pfc-atlas-qtl/QTL/get-pval-correction.py
done