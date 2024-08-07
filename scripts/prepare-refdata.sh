#!/usr/bin/env bash

module load GATK
module load samtools

REF='/fdb/cellranger-arc/refdata-cellranger-arc-GRCh38-2020-A/fasta/genome.fa'
mkdir -p "/data/${USER}/.hg38"
if [[ ! -f "/data/${USER}/.hg38/hg38.fa" ]]; then
    ln -s ${REF} "/data/${USER}/.hg38/hg38.fa"
fi

if [[ ! -f "/data/${USER}/.hg38/hg38.dict" ]]; then
    gatk CreateSequenceDictionary R="/data/${USER}/.hg38/hg38.fa"
fi

if [[ ! -f "/data/${USER}/.hg38/hg38.fa.fai" ]]; then
    samtools faidx "/data/${USER}/.hg38/hg38.fa"
fi

if [[ ! -f "/data/${USER}/.hg38/hg38_chr.map" ]]; then
    wget -O "/data/${USER}/.hg38/hg38_chr.map" https://github.com/naumanjaved/fingerprint_maps/raw/master/map_files/hg38_chr.map
fi

# Build new map
function reorderMap() {
    local fn=${1}
    local out_fn=${fn/map/reorder.map}
    echo "$fn being reordered to $out_fn"
    grep '^@SQ' ${fn}  > header.txt
    grep 'chr1\s' header.txt > ${out_fn}
    grep 'chr10\s' header.txt >> ${out_fn}
    grep 'chr11\s' header.txt >> ${out_fn}
    grep 'chr12\s' header.txt >> ${out_fn}
    grep 'chr13\s' header.txt >> ${out_fn}
    grep 'chr14\s' header.txt >> ${out_fn}
    grep 'chr15\s' header.txt >> ${out_fn}
    grep 'chr16\s' header.txt >> ${out_fn}
    grep 'chr17\s' header.txt >> ${out_fn}
    grep 'chr18\s' header.txt >> ${out_fn}
    grep 'chr19\s' header.txt >> ${out_fn}
    grep 'chr2\s' header.txt >> ${out_fn}
    grep 'chr20\s' header.txt >> ${out_fn}
    grep 'chr21\s' header.txt >> ${out_fn}
    grep 'chr22\s' header.txt >> ${out_fn}
    grep 'chr3\s' header.txt >> ${out_fn}
    grep 'chr4\s' header.txt >> ${out_fn}
    grep 'chr5\s' header.txt >> ${out_fn}
    grep 'chr6\s' header.txt >> ${out_fn}
    grep 'chr7\s' header.txt >> ${out_fn}
    grep 'chr8\s' header.txt >> ${out_fn}
    grep 'chr9\s' header.txt >> ${out_fn}
    grep '#' ${fn} >> ${out_fn}
    grep '^chr1\s' ${fn} >> ${out_fn}
    grep '^chr10\s' ${fn} >> ${out_fn}
    grep '^chr11\s' ${fn} >> ${out_fn}
    grep '^chr12\s' ${fn} >> ${out_fn}
    grep '^chr13\s' ${fn} >> ${out_fn}
    grep '^chr14\s' ${fn} >> ${out_fn}
    grep '^chr15\s' ${fn} >> ${out_fn}
    grep '^chr16\s' ${fn} >> ${out_fn}
    grep '^chr17\s' ${fn} >> ${out_fn}
    grep '^chr18\s' ${fn} >> ${out_fn}
    grep '^chr19\s' ${fn} >> ${out_fn}
    grep '^chr2\s' ${fn} >> ${out_fn}
    grep '^chr20\s' ${fn} >> ${out_fn}
    grep '^chr21\s' ${fn} >> ${out_fn}
    grep '^chr22\s' ${fn} >> ${out_fn}
    grep '^chr3\s' ${fn} >> ${out_fn}
    grep '^chr4\s' ${fn} >> ${out_fn}
    grep '^chr5\s' ${fn} >> ${out_fn}
    grep '^chr6\s' ${fn} >> ${out_fn}
    grep '^chr7\s' ${fn} >> ${out_fn}
    grep '^chr8\s' ${fn} >> ${out_fn}
    grep '^chr9\s' ${fn} >> ${out_fn}
    rm header.txt
}

if [[ ! -f "/data/${USER}/.hg38/hg38_chr.reorder.map" ]]; then
    reorderMap "/data/${USER}/.hg38/hg38_chr.map"
fi
