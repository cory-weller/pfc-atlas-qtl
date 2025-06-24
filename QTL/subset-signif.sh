#!/usr/bin/env bash
#SBATCH --partition quick
#SBATCH --time 3:30:00
#SBATCH --mem 12G

dir=${1}
cd $dir

ml R/4.3

threshold=$(Rscript - <<EOF
library(data.table)
dat <- fread('pval-correction-mapping.txt')

dat <- dat[p_Bonferroni < 0.05]

p_threshold <- max(dat[['p_nominal']])

cat(p_threshold)

EOF
)

head -n 1 chr1.csv > bonferroni-significant-qtls.csv

for file in chr*.csv; do
    awk -F',' -v T=$threshold '$7 < T' ${file} >> bonferroni-significant-qtls.csv
done