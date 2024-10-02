#!/usr/bin/env bash

wget -O data/00-common_all.vcf.gz https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-common_all.vcf.gz

zgrep -v -F '#' data/00-common_all.vcf.gz | awk '{print $1,$2,$3,$4,$5}' | gzip -c > rsids.txt.gz