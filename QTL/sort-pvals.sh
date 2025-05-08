#!/usr/bin/env bash
#SBATCH --partition quick
#SBATCH --mem 40G


file=${1}

echo $file
awk -F',' 'NR>1 {print $7}' $file | sort -g > ${file}.sort

# parallel -j 1 sbatch sort-pvals.sh {} ::: *-ALL.csv
