# pfc-atlas-qtl

# Files

`ATAC` caQTL blacklist files were prepared using 1) RefSeq TSS positions,
and 2) from the Boyle Lab [ENCODE blacklist](https://github.com/Boyle-Lab/Blacklist/) to yield
[`TSS-blacklist.bed`](data/TSS-blacklist.bed), [`boyle-blacklist.bed`](data/boyle-blacklist.bed),
and [`boyle-plus-TSS-blacklist.bed`](data/boyle-plus-TSS-blacklist.bed)

# Fingerprinting

See [GATK documentation](https://gatk.broadinstitute.org/hc/en-us/articles/360037594711-CrosscheckFingerprints-Picard)

## 0. Prepare reference data
The tool needs a haplotype map file (VCF style) with identical chromosome order as your other inputs.
I manually reordered the map with a series of `grep` commands, as well as prepare reference and
sequence dictionary by running [prepare-refdata.sh](scripts/prepare-refdata.sh).

## 1. Get RNA fingerprints from `BAM`
Generate a two-column file that contains `SAMPLEID` and `FULL_BAMFILE_PATH` (no header), in my case [`sample_bam.txt`](sample_bam.txt).

Then submit job array with one job per RNA `bam`, each running [get-rna-fingerprints.sh](scripts/get-rna-fingerprints.sh)
```bash
mkdir -p fingerprints/rna
samplefile='sample_bam.txt'
njobs=$(wc -l $samplefile)
sbatch --array=1-${njobs}%50 scripts/get-rna-fingerprints.sh ${samplefile}
```

## 2. Combine RNA fingerprints into one file
Once all jobs complete and separate fingrprints files exist (one per sample), merge all into a single file `rna-merged.vcf.gz`
```bash
cd fingerprints/rna
module load samtools
ls *.vcf.gz > files.txt
bcftools merge --file-list files.txt -o rna-merged.vcf
bgzip rna-merged.vcf 
tabix -p vcf rna-merged.vcf.gz
```

## 3. Get DNA fingerprints from plink `bed/bim/fam`
[`merge-nabec-hbcc-genotypes.sh`](scripts/merge-nabec-hbcc-genotypes.sh) generates fingerprints file `dna-merged.vcf.gz` 

The script [`merge-nabec-hbcc-genotypes.sh`](scripts/merge-nabec-hbcc-genotypes.sh) is currently _very_ specific to this dataset. 
```bash
bash scripts/merge-nabec-hbcc-genotypes.sh
```

## 4. Calculate fingerprints on whole set
```bash
rna='fingerprints/rna/rna-merged.vcf.gz'
dna='fingerprints/dna/dna-merged.vcf.gz'
map='hg38_chr.reorder.map'
output='fingerprints/crosscheck/all-pairs.txt'
sbatch scripts/crosscheck.sh ${rna} ${dna} ${map} ${output}
```

## Plot results
The script [`plot-crosschecks.R`](scripts/plot-crosschecks.R) is _very_ specific to the samples ran in this project, so modify as needed.
```bash
module load R/4.3
Rscript scripts/plot-crosschecks.R
```


# QTL Analysis
Generate table of samples along with batch and `bam` location
```bash
module load R/4.3 && \
Rscript scripts/finalize-samples.R
```

Pseudobulked counts in `/data/CARD_singlecell/brain_atlas_wnn/output/rna/`

## Prepare genotypes
See [`genotypes`](genotypes) directory's `README` file.


## Plot covarites along Principal Component axes
```bash
module load R/4.3 && \
Rscript scripts/plot-covariates.R
```

<details>
    <summary>HBCC covariate plots</summary

![](PCA-plots/HBCC-Age.png)
![](PCA-plots/HBCC-Sex.png)
![](PCA-plots/HBCC-Ancestry.png)
![](PCA-plots/HBCC-Homogenization.png)
![](PCA-plots/HBCC-LibraryPrep.png)
![](PCA-plots/HBCC-Sequencing.png)

</details>


<details>
    <summary>NABEC covariate plots</summary

![](PCA-plots/NABEC-Age.png)
![](PCA-plots/NABEC-Sex.png)
![](PCA-plots/NABEC-Ancestry.png)
![](PCA-plots/NABEC-Homogenization.png)
![](PCA-plots/NABEC-LibraryPrep.png)
![](PCA-plots/NABEC-Sequencing.png)

</details>

Based on separation along principal component #1 for cohort HBCC, these 8 samples were excluded
by running [`remove-hbcc-outliers.sh`](scripts/remove-hbcc-outliers.sh) to generate `genotypes/HBCC_polarized_nooutliers.{bed,bim,fam}`.


| FID       | PC1    |
| --------- | ------ |
| HBCC_1058 | 0.2968 |
| HBCC_1331 | 0.3027 |
| HBCC_1385 | 0.3055 |
| HBCC_1431 | 0.2655 |
| HBCC_1560 | 0.3043 |
| HBCC_2429 | 0.3139 |
| HBCC_2756 | 0.3283 |
| HBCC_2781 | 0.2680 |


## Prepare final tensorQTL covariates files
See [`prep-covariates.R`](scripts/prep-covariates.R). The script generates a data frame for subsetting with particular cell type QTL runs.
```bash
module load R/4.3 && \
Rscript scripts/prep-covariates.R
```

## Format expression counts as BED format
See [`prep-QTL-bedfiles.R`](scripts/prep-QTL-bedfiles.R). Generates a file `{mode}-{cohort}-{celltype}-counts.bed` within the
directory `QTL-pseudobulk-counts`
```bash
module load R/4.3 && \
Rscript scripts/prep-QTL-bedfiles.R
```

## Manually add feature blacklists
The file `data/array-params.tsv` can include four columns. By default it will be generated with the
first three. The fourth column specifying a feature blacklist can be added manually. I used `awk` 
to add `data/TSS-blacklist.bed` for the ATAC runs that used the blacklist.

```bash
awk '{OFS="\t"; if ($0~/atac/) {print $1,$2,$3,"data/TSS-blacklist.bed"} else {print $0,""}}'  data/array-params.tsv \
> .tmpparams && \
mv .tmpparams data/array-params.tsv
```


| Column    | Value  |
| --------- | ------ |
| 1 | Celltype, `{Astro,ExN,InN,MG,OPC,Oligo,VC}` |
| 2 | Mode of data, . `{rna, atac}` |
| 3 | Cohort, `{HBCC,NABEC}` |
| 4 | Blacklist `bed` file with (at minimum) headers `{chr,start,stop}` |

## Run TensorQTL

See [`run-tensorQTL.sh`](scripts/run-tensorQTL.sh). It takes a single argument, a `$BATCHNAME` that
will create a folder for outputs, i.e. `QTL-output/$BATCHNAME`. 
The script must be submitted as a job array. Each job corresponds to a row from `data/array-params.tsv`.

Briefly, the script does the following:
1. Imports the `N`th row (of job array, using `$SLURM_ARRAY_TASK_ID`) from `data/array-params.tsv`
2. Creates and uses a temporary `lscratch` working directory
3. Executes [`intersect-files.R`](scripts/intersect-files.R) to generate run-specific subset of pseudobulk counts, covariates, and interaction terms
4. Executes [`subset-plink.sh`](scripts/subset-plink.sh) to generate run-specific subset of genotypes
5. Executes `tensorqtl` job using singularity container
6. Copies output from `lscratch` to `$BATCHNAME` within this project directory


```bash
sbatch --array=1-28%4 scripts/run-tensorQTL-nominal.sh scvicounts-nominal-oct2
sbatch --array=1-28%4 scripts/run-tensorQTL-interaction.sh scvicounts-interaction-oct2
```

# Analyze QTL results
Retrieved list of rsIDs with [`get-rsids.sh`](scripts/get-rsids.sh)
```bash
bash scripts/get-rsids.sh
```

# Combine results
See [`rbind-qtl-results.R`](scripts/rbind-qtl-results.R) which generates three files:
- `cis-caQTL.tsv`
- `cis-eQTL.tsv`
- `cis-QTL-combined.tsv`
```bash
Rscript scripts/rbind-qtl-results.R
```



# eQTL Analysis
To assess expression across genotype, we need a file with alternate allele dosage per variant.
[`get-alt-allele-dosage.sh`](scripts/get-alt-allele-dosage.sh) extracts such info from `plink` files
to generate `.traw` files (transposed raw plink genotypes) where columns=samples, rows=variants.

Get list of variants `variant-ids.txt` associated with ATAC or RNA data with `scripts/analyze-eqtls.R`
