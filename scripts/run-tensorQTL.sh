#!/usr/bin/env bash
#SBATCH --mem 24g
#SBATCH --gres gpu:v100x:1,lscratch:50
#SBATCH --cpus-per-task 56
#SBATCH --time=4:00:00
#SBATCH --partition gpu
#SBATCH --output=logs/%j.out


BATCHNAME=${1}

module load singularity/4
module load CUDA/11.4

cd /data/CARD_singlecell/users/wellerca/pfc-atlas-qtl
PROJDIR=$(realpath ${PWD})
mkdir -p QTL-output/${BATCHNAME}

PARAMS="${PROJDIR}/data/array-params.tsv"

N=${SLURM_ARRAY_TASK_ID}


CELLTYPE=$(sed -n ${N}p ${PARAMS} | cut -f 1)
MODE=$(sed -n ${N}p ${PARAMS} | cut -f 2)
COHORT=$(sed -n ${N}p ${PARAMS} | cut -f 3)
EXCLUDE=$(sed -n ${N}p ${PARAMS} | cut -f 4)

TMPDIR="/lscratch/${SLURM_JOB_ID}"

mkdir -p ${TMPDIR} && cd ${TMPDIR}

module load R/4.3

MYCMD=(Rscript ${PROJDIR}/scripts/intersect-files.R \
    --cohort ${COHORT} \
    --mode ${MODE} \
    --celltype ${CELLTYPE} \
    --projdir ${PROJDIR})

if [[ "${EXCLUDE}" != '' ]]; then
    MYCMD+=( "--exclude" "${EXCLUDE}" )
fi
    
echo "Running command:"
echo "${MYCMD[@]}"

${MYCMD[@]}

## Subset plink files
# Link original genotype files
if [[ ${COHORT} == 'HBCC' ]]; then
    ln -s ${PROJDIR}/genotypes/HBCC_polarized_nooutliers.bed genotypes.bed
    ln -s ${PROJDIR}/genotypes/HBCC_polarized_nooutliers.bim genotypes.bim
    ln -s ${PROJDIR}/genotypes/HBCC_polarized_nooutliers.fam genotypes.fam
fi

if [[ ${COHORT} == 'NABEC' ]]; then
    ln -s ${PROJDIR}/genotypes/NABEC_polarized.bed genotypes.bed
    ln -s ${PROJDIR}/genotypes/NABEC_polarized.bim genotypes.bim
    ln -s ${PROJDIR}/genotypes/NABEC_polarized.fam genotypes.fam
fi

bash ${PROJDIR}/scripts/subset-plink.sh 






PREFIX="${COHORT}_${CELLTYPE}_${MODE}"
PLINK=genotypes-forqtl
COUNTS=pseudobulk-counts.bed
COVS=covariates.txt
INTERACTION="interaction.tsv"
OUTPUT=${COHORT}-${MODE}-${CELLTYPE}
tqtl_img='/usr/local/apps/tensorqtl/1.0.9/libexec/TensorQTL-1.0.9_from_docker.sif'

mkdir -p ${PREFIX}

singularity exec -B ${PWD} --no-home --nv ${tqtl_img} tensorqtl \
    ${PLINK} ${COUNTS} ${PREFIX} \
    --mode cis_nominal \
    --seed 2024 \
    --output_dir ${PREFIX} \
    --covariates ${COVS} \
    --interaction ${INTERACTION}

cp -r ${PREFIX} ${PROJDIR}/QTL-output/${BATCHNAME}

