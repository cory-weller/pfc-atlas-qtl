# pfc-atlas-qtl

# Covariates
In addition to `sex`, we will include principal components accounting for genetic ancestry
(`genetics_PCs`) and chromatin accessibility (`ATAC_PCs`).

## Genetics PCs
Genotypes from whole genome sequencing, LD-pruned, and principal components calculated.

## ATAC PCs
Pseudobulked ATAC features per `donor` x `cell_type` combination, principal components calculated.


# Formatting


Raph’s original notebook link
https://github.com/FOUNDINPD/foundin_qtl


https://github.com/FOUNDINPD/foundin_qtl/blob/main/analyses/cis_tensorqtl.ipynb


/data/CARDPB/projects/daidak2/singlecell_QTL/notebooks/cis_qtl_tensorqtl_SV_covariates_2024_paper.ipynb

##  Prepare NABEC genotypes
```bash
cp /data/CARDPB/projects/daidak2/singlecell_QTL/genotypes/NABEC/SNV_illumina/MERGED_MAF_GENO005_plink19_NABEC_scRNA.{bed,bim,fam} data/genotypes
for file in $(find data/genotypes/ -name 'MERGED_MAF*'); do
    newname=$(echo $file | sed 's/MERGED_.*scRNA/NABEC/g')
    mv $file $newname
done

mv data/genotypes/NABEC.bed data/genotypes/NABEC-forQTL.bed
mv data/genotypes/NABEC.bim data/genotypes/NABEC-forQTL.bim
mv data/genotypes/NABEC.fam data/genotypes/NABEC-forQTL.fam
```
## Prepare HBCC genotypes
```bash
# Copy files
cp /data/CARDPB/projects/daidak2/singlecell_QTL/genotypes/HBCC/hbcc_DV_gvcf.deepvariant_SC_MAF_GENO_005_HWE_0001_updateid_plink19.{bed,bim,fam} data/genotypes
for file in $(find data/genotypes/ -name 'hbcc_*'); do
    newname=$(echo $file | sed 's/hbcc_.*plink19/HBCC_oldID/g')
    mv $file $newname
done

# Rename IDs
awk '{print $1,$1,$2,$2}' /data/CARDPB/projects/daidak2/singlecell_QTL/sample_info/HBCC_rename.tsv > data/HBCC_id_update.txt
plink --bfile data/genotypes/HBCC_oldID \
    --make-bed \
    --keep-allele-order \
    --update-ids data/HBCC_id_update.txt \
    --out data/genotypes/HBCC
```

## Copy covariates
```bash
cp /data/CARD_singlecell/brain_atlas_wnn/input/PFC_covariates.csv data/covariates.csv
```

## LD Prune genotypes
```bash
bash scripts/prep-genotypes.sh data/genotypes/HBCC
bash scripts/prep-genotypes.sh data/genotypes/NABEC
```


## Plot covarites along Principal Component axes
```bash
Rscript scripts/plot-covariates.R
```


<details>
    <summary>HBCC covariate plots</summary

![](plots/HBCC-Age.png)
![](plots/HBCC-Sex.png)
![](plots/HBCC-Ancestry.png)
![](plots/HBCC-Homogenization.png)
![](plots/HBCC-LibraryPrep.png)
![](plots/HBCC-Sequencing.png)

</details>


<details>
    <summary>NABEC covariate plots</summary

![](plots/NABEC-Age.png)
![](plots/NABEC-Sex.png)
![](plots/NABEC-Ancestry.png)
![](plots/NABEC-Homogenization.png)
![](plots/NABEC-LibraryPrep.png)
![](plots/NABEC-Sequencing.png)

</details

Based on separation along principal component #1 for cohort HBCC, these 8 samples were excluded:
```
         FID    PC1
      <char>  <num>
1: HBCC_1058 0.2968
2: HBCC_1331 0.3027
3: HBCC_1385 0.3055
4: HBCC_1431 0.2655
5: HBCC_1560 0.3043
6: HBCC_2429 0.3139
7: HBCC_2756 0.3283
8: HBCC_2781 0.2680
```

Generate [`data/HBCC-remove.txt`](data/HBCC-remove.txt) for removing with `plink`:

```bash
echo "
awk '$7 > 0.2 {print $1,$2}' data/genotypes/HBCC.pruned_pc.txt > data/HBCC-remove.txt


module load plink/1.9.0-beta4.4 
plink --bfile data/genotypes/HBCC \
      --remove data/HBCC-remove.txt \
      --make-bed \
      --out data/genotypes/HBCC-forQTL
```



## Prepare final tensorQTL covariates files

`prep-covariates.R`


# Format phenotypes as BED format

`prep-QTL-bedfiles.R`



# Generate final files with appropriate intersections of samples
```bash
# Subset everything except plink files
```bash
Rscript scripts/intersect-files.sh
```

# Subset plink files
```bash
bash scripts/subset-plink.sh
```
```
# Run tensorQTL
```bash
sbatch --array=2-24%8 scripts/run-tensorQTL.sh
```

# Combine results

```bash
Rscript scripts/rbind-qtl-results.R

# cis-caQTL.tsv
# cis-eQTL.tsv
# cis-QTL-combined.tsv
```

