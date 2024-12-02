#!/usr/bin/env bash
#SBATCH --mem 24g
#SBATCH --cpus-per-task 2
#SBATCH --time=3:59:00
#SBATCH --partition quick,norm

# If GPU needed in future...
####SBATCH --gres gpu:v100x:1,lscratch:50

module load tensorqtl

windowsize=1000000

N=${SLURM_ARRAY_TASK_ID}

paramsfile='data/array-params.txt'

celltype=$(sed -n ${N}p ${paramsfile}  | awk '{print $1}')
mode=$(sed -n ${N}p ${paramsfile} | awk '{print $2}')
cohort=$(sed -n ${N}p ${paramsfile} | awk '{print $3}')
chr=$(sed -n ${N}p ${paramsfile} | awk '{print $4}')


for method in mean sum; do
    # run Nominal, no interaction
    python-tqtl /data/CARD_singlecell/users/wellerca/pfc-atlas-qtl/run-tensorqtl.py \
        --mode ${mode} \
        --celltype ${celltype} \
        --covariates PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 Sex \
        --cohort ${cohort} \
        --interaction None \
        --qtlmethod nominal \
        --outdir QTL-output \
        --pseudobulkmethod ${method} \
        --window ${windowsize} \
        --chr chr${chr}

    # run Nominal, interaction with Age
    python-tqtl /data/CARD_singlecell/users/wellerca/pfc-atlas-qtl/run-tensorqtl.py \
        --mode ${mode} \
        --celltype ${celltype} \
        --covariates PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 Sex \
        --cohort ${cohort} \
        --interaction Age \
        --qtlmethod nominal \
        --outdir QTL-output \
        --pseudobulkmethod ${method} \
        --window ${windowsize} \
        --chr chr${chr}

done

