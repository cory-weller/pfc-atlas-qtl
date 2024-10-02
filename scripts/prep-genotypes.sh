#!/usr/bin/env bash
# Generates plink files and genotype PCs for meffil analysis

# Check plink
if ! command -v plink &> /dev/null; then
    echo "INFO: plink command not found, looking for plink module"
    if ! command -v module &> /dev/null; then
        echo "ERROR: module command not found. Did you mean to run this on an HPC?"
        exit 1
    else
        module load plink/1.9
    fi
fi

# Check king
if ! command -v king &> /dev/null; then
    echo "INFO: king command not found, looking for king module"
    if ! command -v module &> /dev/null; then
        echo "ERROR: module command not found. Did you mean to run this on an HPC?"
        exit 1
    else
        module load king/2.2.7
    fi
fi

bedfile=${1}
bfile=${bedfile%.bed}


# Generates plink.prune.in for pruned loci for GRM/PCA
prune_fn="${bfile}.prune.in"
if [[ ! -f "${prune_fn}" ]]; then
    plink \
        --bfile ${bfile} \
        --indep-pairwise 1000 10 0.02 && \
    mv plink.prune.in ${prune_fn} && \
    rm plink.prune.out plink.log plink.nosex
else
    echo "List of LD-Pruned snps already exists,  skipping"
fi

pruned_bfile="${bfile}.pruned"
if [[ ! -f "${pruned_bfile}.bed" ]]; then
    plink \
        --bfile ${bfile} \
        --make-bed \
        --extract ${prune_fn} \
        --keep-allele-order \
        --out ${pruned_bfile}
else
    echo "LD-Pruned plink file already exists, skipping"
fi

# Generate PCs for population structure
if [[ ! -f "${pruned_bfile}_pc.txt" ]]; then
    king -b ${pruned_bfile}.bed \
        --pca \
        --prefix ${pruned_bfile}_
else
    echo "Genetic PCs from king already exist, skipping"
fi

echo "Done!"

exit 0