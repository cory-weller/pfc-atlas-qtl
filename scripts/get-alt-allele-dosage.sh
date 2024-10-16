#!/usr/bin/env bash

PROJDIR='/data/CARD_singlecell/users/wellerca/pfc-atlas-qtl' && \
cd ${PROJDIR}


HBCC='genotypes/HBCC_polarized_nooutliers'
HBCC=$(realpath $HBCC)
NABEC='genotypes/NABEC_polarized'
NABEC=$(realpath $NABEC)

# work in lscratch tmp directory
TMPDIR="/lscratch/${SLURM_JOB_ID}"
mkdir -p ${TMPDIR} && cd ${TMPDIR}

# Get HBCC .traw
module load plink/1.9
plink --bfile ${HBCC} \
    --recode A-transpose \
    --keep-allele-order \
    --output-chr chrM \
    --out HBCC-alt-dosage

# Get NABEC .traw
plink --bfile ${NABEC} \
    --recode A-transpose \
    --keep-allele-order \
    --output-chr chrM \
    --out NABEC-alt-dosage

sed -i -E '1 s/_HBCC_[0-9]+//g' HBCC-alt-dosage.traw
sed -i -E '1 s/_[A-Z]+\-[0-9]+//g' NABEC-alt-dosage.traw

mv HBCC-alt-dosage.traw ${PROJDIR}/genotypes/
mv NABEC-alt-dosage.traw ${PROJDIR}/genotypes/
