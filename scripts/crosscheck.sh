#!/usr/bin/env bash
#SBATCH --time 1:30:00
#SBATCH --mem 60g
#SBATCH --partition quick
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4

module load GATK

RNAFILE=${1}
DNAFILE=${2}
MAPFILE=${3}
OUTPUTFILE=${4}

if [[ -z ${RNAFILE} ]]; then echo '4 arguments required: RNAFILE DNAFILE MAPFILE OUTPUTFILE'; exit 1; fi
if [[ -z ${DNAFILE} ]]; then echo '4 arguments required: RNAFILE DNAFILE MAPFILE OUTPUTFILE'; exit 1; fi
if [[ -z ${MAPFILE} ]]; then echo '4 arguments required: RNAFILE DNAFILE MAPFILE OUTPUTFILE'; exit 1; fi
if [[ -z ${OUTPUTFILE} ]]; then echo '4 arguments required: RNAFILE DNAFILE MAPFILE OUTPUTFILE'; exit 1; fi

gatk CrosscheckFingerprints \
    --INPUT ${DNAFILE} \
    --SECOND_INPUT ${RNAFILE} \
    --HAPLOTYPE_MAP ${MAPFILE} \
    --CROSSCHECK_BY SAMPLE \
    --CROSSCHECK_MODE CHECK_ALL_OTHERS \
    --OUTPUT ${OUTPUTFILE}

