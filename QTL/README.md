# README

Generating BH-adjusted p-values for ATAC required some on-disk sorting:

1. generate a separate p-value file for every chromosome, as done for rna as well.
2. sort nominal p-values with `sort-pvals.sh` to generate a single file per chromosome.
3. use `get-pval-correction.py` to generate a file that has three columns:
    - `p_nominal` `p_bh` `p_Bonferroni`
4. subset significant rows with `subset-signif.sh`