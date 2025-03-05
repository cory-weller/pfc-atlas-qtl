#!/usr/bin/env bash
#SBATCH --time 3:59:00
#SBATCH --mem 200G
#SBATCH --nodes 1
#SBATCH --partition quick,norm
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4

module unload
module load R/4.3

export QTLMETHOD='sum'
Rscript scripts/parse-parquet-files.R