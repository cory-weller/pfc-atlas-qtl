#!/usr/bin/env bash
#SBATCH --partition quick
#SBATCH --time 3:30:00
#SBATCH --mem 12G


dir=${1}
cd $dir

for file in *.csv; do
    if [[ ! -f ${file}.sort ]]; then
        echo $file
        awk -F',' 'NR>1 {print $7}' $file | sort -g > ${file}.sort
    fi
done

sort -g -m --batch-size=22 *.csv.sort > chrALL.sorted.txt

# parallel -j 1 echo sbatch sort-pvals.sh {} ::: *-atac
