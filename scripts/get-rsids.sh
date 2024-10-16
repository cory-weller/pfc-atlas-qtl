#!/usr/bin/env bash

vcf='data/00-common_all.vcf.gz'
url='https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-common_all.vcf.gz'
rsids='data/rsids.txt.gz'

if [[ ! -f $vcf ]]; then
    wget -O $vcf $url
fi

if [[ ! -f $rsids ]]; then
    zgrep -v -F '#' $vcf | awk '{print $1,$2,$3,$4,$5}' | gzip -c > $rsids
fi
