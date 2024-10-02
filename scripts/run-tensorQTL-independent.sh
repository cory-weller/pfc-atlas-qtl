#!/usr/bin/env bash
#SBATCH --mem 24g
#SBATCH --gres gpu:v100x:1,lscratch:50
#SBATCH --cpus-per-task 56
#SBATCH --time=4:00:00
#SBATCH --partition gpu

module load singularity/4
module load CUDA/11.4

cd /data/CARD_singlecell/users/wellerca/pfc-atlas-qtl

N=${SLURM_ARRAY_TASK_ID}


PARAMS='data/array-params.tsv'

CELLTYPE=$(sed -n ${N}p ${PARAMS} | cut -f 1)
MODE=$(sed -n ${N}p ${PARAMS} | cut -f 2)
COHORT=$(sed -n ${N}p ${PARAMS} | cut -f 3)

PREFIX="${COHORT}_${CELLTYPE}_${MODE}"
PLINK="tensorqtl-subsets/${COHORT}-${MODE}-${CELLTYPE}-plink"
COUNTS="tensorqtl-subsets/${COHORT}-${MODE}-${CELLTYPE}-counts.bed"
COVS="tensorqtl-subsets/${COHORT}-${MODE}-${CELLTYPE}-covariates.txt"
INTERACTION="tensorqtl-subsets/${COHORT}-${MODE}-${CELLTYPE}-interaction.tsv"

tqtl_img='/usr/local/apps/tensorqtl/1.0.9/libexec/TensorQTL-1.0.9_from_docker.sif'

OUTDIR='tensorqtl_permutations'
mkdir -p ${OUTDIR}

singularity exec -B ${PWD} --no-home --nv ${tqtl_img} tensorqtl \
    ${PLINK} ${COUNTS} ${PREFIX} \
    --mode cis \
    --seed 2024 \
    --output_dir ${OUTDIR} \
    --covariates ${COVS}

singularity exec -B ${PWD} --no-home --nv ${tqtl_img} tensorqtl \
    ${PLINK} ${COUNTS} ${PREFIX} \
    --mode cis_independent \
    --seed 2024 \
    --output_dir ${OUTDIR} \
    --covariates ${COVS} \
    --cis_output ${OUTDIR}/${PREFIX}.cis_qtl.txt.gz

