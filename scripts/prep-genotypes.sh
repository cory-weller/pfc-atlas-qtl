#!/usr/bin/env bash

module load plink/1.9

###############################################################################
# Section 1
# Rename HBCC IDs and filter both cohorts to only autosomal biallelic SNPs 
###############################################################################
plink --bfile HBCC_oldID \
    --biallelic-only \
    --snps-only \
    --not-chr X,Y,chrX,chrY,23,chr23,M,MT,chrM,chrMT \
    --make-bed \
    --update-ids HBCC_id_update.txt \
    --remove <(echo -e "HBCC-1832 HBCC-1832\nHBCC-887 HBCC-887") \
    --out HBCC_biallelic_autosomal

plink --bfile NABEC \
    --biallelic-only \
    --snps-only \
    --not-chr X,Y,chrX,chrY,23,chr23,M,MT,chrM,chrMT \
    --make-bed \
    --out NABEC_biallelic_autosomal





###############################################################################
# Section 2
# merge cohorts into a single plink file to get shared REF/ALT allele order
###############################################################################
# First-pass merge 
# This prints (in log file) list of SNPs that are triallelic when combined into one dataset
plink -bfile NABEC_biallelic_autosomal \
    --bmerge HBCC_biallelic_autosomal \
    --make-bed \
    --out triallelic_test

# Extract triallelic sites into text file
grep 'Warning: Variants' triallelic_test.log | \
    cut -d ' ' -f 3,5 | tr -d "'" | sed 's/ /\n/g' > triallelic.txt && \
    find . -name 'triallelic_test.*' -delete

# Remove triallelic SNPs from NABEC
plink -bfile NABEC_biallelic_autosomal \
    --exclude triallelic.txt \
    --make-bed \
    --out NABEC_biallelic_autosomal2 && \
    find . -name 'NABEC_biallelic_autosomal.*' -delete

# Remove triallelic SNPs from HBCC
plink -bfile HBCC_biallelic_autosomal \
    --exclude triallelic.txt \
    --make-bed \
    --out HBCC_biallelic_autosomal2 && \
    find . -name 'HBCC_biallelic_autosomal.*' -delete

# Redo merge without triallelic SNPs
plink -bfile HBCC_biallelic_autosomal2 \
    --bmerge NABEC_biallelic_autosomal2 \
    --merge-equal-pos \
    --make-bed \
    --out HBCC_plus_NABEC_biallelic_autosomal


###############################################################################
# Section 3
# pulls out cohort-specific files while retaining shared REF/ALT allele order
###############################################################################


# Generate HBCC-specific file with proper filters and retain same allele order as NABEC
plink --bfile HBCC_plus_NABEC_biallelic_autosomal \
    --keep HBCC_biallelic_autosomal2.fam \
    --maf 0.05 \
    --hwe 1e-6 \
    --geno 0.05 \
    --make-bed \
    --keep-allele-order \
    --out HBCC_polarized

# Generate NABEC-specific file with proper filters and retain same allele order as HBCC
plink --bfile HBCC_plus_NABEC_biallelic_autosomal \
    --keep NABEC_biallelic_autosomal2.fam \
    --maf 0.05 \
    --hwe 1e-6 \
    --geno 0.05 \
    --make-bed \
    --keep-allele-order \
    --out NABEC_polarized

find . -name '*_biallelic_autosomal*' -delete


###############################################################################
# Section 4
# LD pruning and calculating principal components
###############################################################################

calc_pruned_pcs() {
    
    bedfile=${1}
    bfile=${bedfile%.bed}
    prune_fn="${bfile}.prune.in"
    pruned_bfile="${bfile}.pruned"

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



    # Generates plink.prune.in for pruned loci for GRM/PCA
    if [[ ! -f "${prune_fn}" ]]; then
        plink \
            --bfile ${bfile} \
            --indep-pairwise 1000 10 0.02 && \
        mv plink.prune.in ${prune_fn} && \
        rm plink.prune.out plink.log plink.nosex
    else
        echo "List of LD-Pruned snps already exists,  skipping"
    fi

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

}

calc_pruned_pcs HBCC_polarized
calc_pruned_pcs NABEC_polarized

for cohort in HBCC NABEC; do
    plink --bfile ${cohort}_polarized.pruned \
        --keep-allele-order \
        --recode A-transpose \
        --output-chr chrM \
        --out ${cohort}_pruned_for_pcs
done



echo "Done LD-pruning and calculating PCs"



###############################################################################
# Section 5
# Generate transposed raw genotype dosage files and compress to gzip
###############################################################################
for cohort in HBCC NABEC; do
    for i in $(seq 1 22); do
        chr="chr$i"
        plink --bfile ${cohort}_polarized \
            --keep-allele-order \
            --recode A-transpose \
            --output-chr chrM \
            --chr ${chr} \
            --out ${cohort}_${chr}
        
        sed -r '1s/_[^\t]+//g' ${cohort}${chr}.traw | pigz -c > ${cohort}_${chr}.traw.gz && rm ${cohort_chr}.traw
        pigz ${cohort_chr}.traw
    done
done

# Single genome-wide file:
for cohort in HBCC NABEC; do
    plink --bfile ${cohort}_polarized \
        --keep-allele-order \
        --recode A-transpose \
        --output-chr chrM \
        --chr chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22 \
        --out ${cohort}_autosomes
    
    # Fix header to only include IID (instead of FID_IID)
    sed -r '1s/_[^\t]+//g' ${cohort}_autosomes.traw | pigz -c > ${cohort}_autosomes.traw.gz && rm ${cohort}_autosomes.traw
done


find . -name '*.nosex' -delete
find . -name '*.log' -delete


echo "Done generating allele dosage files"

