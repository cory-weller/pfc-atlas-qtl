#!/usr/bin/env bash
#SBATCH --partition quick
#SBATCH --time 3:30:00
#SBATCH --mem 12G

dir=${1}

cd $dir

for file in chr*.csv; do
    awk -F',' 'NR>1 {print $2}' ${file} | sort -u >> variants-tested.txt
done
