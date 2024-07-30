# pfc-atlas-qtl



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
```

# Rename IDs
```bash
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
```bash
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

To remove these files, generated [`data/HBCC-remove.txt`](data/HBCC-remove.txt) then excluded samples with `plink`:

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
See [`prep-covariates.R`](scripts/prep-covariates.R)
```bash
Rscript scripts/prep-covariates.R
```

## Format phenotypes as BED format
See [`prep-QTL-bedfiles.R`](scripts/prep-QTL-bedfiles.R)
```bash
Rscript scripts/prep-QTL-bedfiles.R
```


## Generate final files with appropriate intersections of samples
See [`intersect-files.R`](scripts/intersect-files.R) which is executed by [`intersect-files.sh`](scripts/intersect-files.sh)
```bash
# Subset everything except plink files
Rscript scripts/intersect-files.sh
```

## Subset plink files
This script generates a separate `plink` fileset for every combination of `cohort` `celltype` and `mode` (24 total sets).
See [`subset-plink.sh`](scripts/subset-plink.sh)
```bash
bash scripts/subset-plink.sh
```


## Run tensorQTL
See [`run-tensorQTL.sh`](scripts/run-tensorQTL.sh)
```bash
sbatch --array=1-24%8 scripts/run-tensorQTL.sh
```

# Combine results
See [`rbind-qtl-results.R`](scripts/rbind-qtl-results.R) which generates three files:
- `cis-caQTL.tsv`
- `cis-eQTL.tsv`
- `cis-QTL-combined.tsv`
```bash
Rscript scripts/rbind-qtl-results.R
```

