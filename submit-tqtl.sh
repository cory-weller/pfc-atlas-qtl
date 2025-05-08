#!/usr/bin/env bash
#SBATCH --mem 24g
#SBATCH --cpus-per-task 2
#SBATCH --time=3:59:00
#SBATCH --partition quick,norm

# If GPU needed in future...
####SBATCH --gres gpu:v100x:1,lscratch:50

module load tensorqtl

windowsize=1000000

# parallel -j 1 echo {}  ::: Astro ExN InN VC Oligo OPC MG ::: rna  atac ::: HBCC NABEC | sort -k3,3 -k2,2r > array-params.txt
N=${SLURM_ARRAY_TASK_ID}
paramsfile='array-params.txt'
celltype=$(sed -n ${N}p ${paramsfile}  | awk '{print $1}')
mode=$(sed -n ${N}p ${paramsfile} | awk '{print $2}')
cohort=$(sed -n ${N}p ${paramsfile} | awk '{print $3}')

echo $celltype $mode $cohort

python-tqtl /data/CARD_singlecell/users/wellerca/pfc-atlas-qtl/run-tensorqtl.py \
    --traw genotypes/${cohort}_autosomes.traw \
    --counts /data/CARD_singlecell/PFC_atlas/data/celltypes/${celltype}/pseudobulk_rna.csv \
    --counts-index SampleID \
    --outdir QTL/${cohort}-${celltype}-${mode} \
    --window ${windowsize} \
    --covariates genotypes/${cohort}_tqtl_covariates.csv \
    --covariates-index SampleID \
    --mode ${mode}

# HBCC RNA
# sbatch --array=1-7 submit-tqtl.sh

# HBCC ATAC
# sbatch --array=8-14 submit-tqtl.sh

# NABEC RNA
# sbatch --array=15-21 submit-tqtl.sh

# NABEC ATAC
# sbatch --array=22-28 submit-tqtl.sh