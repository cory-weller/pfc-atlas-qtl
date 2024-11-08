#!/usr/bin/env bash
#SBATCH --mem 24g
#SBATCH --cpus-per-task 2
#SBATCH --time=2:00:00
#SBATCH --partition quick

# If GPU needed in future...
####SBATCH --gres gpu:v100x:1,lscratch:50

module load tensorqtl

windowsize=100000

celltype=${1}
mode=${2}
cohort=${3}

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
        --window ${windowsize}

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
        --window ${windowsize}
done

