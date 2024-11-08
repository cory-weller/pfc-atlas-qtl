#!/usr/bin/env Rscript

library(data.table)

rna_features_file <- 'data/rna/RNA-features.tsv'
if(!file.exists(rna_features_file)) {
    # Import cellranger transcriptome
    dat <- fread('/fdb/cellranger/refdata-gex-GRCh38-2020-A/genes/genes.gtf')

    # Take only gene entries
    dat <- dat[V3 == 'gene']

    # Ignore lncRNAs for genes with duplicate entries

    # Pull out `gene_name` from column 9
    dat[, gn := gsub('^.*gene_name', 'SYM:', V9)]
    dat[, gn := gsub(';.*$', '', gn)]
    dat[, gn := gsub('SYM: ', '', gn)]

    # Remove annoying quote marks in string
    dat[, gn := noquote(gn)]
    dat[, gn := gsub('"', '', gn)]
    dat[, gn := as.character(gn)]

    # Subset and rename columns
    dat <- dat[, .SD, .SDcols=c('##date: 2019-09-05', 'V4','V5','gn')]
    dat[, V5 := V4+1]
    setnames(dat, c('chr','start','end','phenotype_id'))
    dat <- dat[chr %in% c(paste0('chr', 1:22))]
    setkey(dat, chr, start, end)

    dat <- dat[!duplicated(phenotype_id)]
    setnames(dat, 'chr','#chr')
    fwrite(dat, file=rna_features_file, quote=F, row.names=F, col.names=T, sep='\t')
} else {
    dat <- fread(rna_features_file)
}


# generate RNA bed files for tensorQTL
# HBCC-samples
# NABEC-samples
modes <- c('rna','atac')
celltypes <- c('ExN','InN','Astro','MG','Oligo','OPC','VC')

mode <- 'rna'
for(celltype in celltypes) {
    countsfile <- paste0('/data/CARD_singlecell/brain_atlas_subtype/output/', mode, '/', celltype, '/cpm_log_pseudobulk_model_counts.csv')
    counts <- fread(countsfile)
    setnames(counts, 'V1', 'symbol')
    samplenames <- grep('-ARC$', colnames(counts), value=TRUE)
    # add 'HBCC_' to samples that begin with number
    dt.tmp <- data.table(orig=samplenames)
    dt.tmp[orig %like% '^[0-9]', v2 := paste0('HBCC_', orig)]
    dt.tmp[orig %like% '^[0-9]', cohort := 'HBCC']
    
    dt.tmp[orig %like% '^[SUK]', v2 := orig]
    dt.tmp[orig %like% '^[SUK]', cohort := 'NABEC']
    dt.tmp[, new := gsub('-ARC$', '', v2)]
    nabec_samples <- dt.tmp[cohort=='NABEC', new]
    hbcc_samples <- dt.tmp[cohort=='HBCC', new]

    setnames(counts, dt.tmp$orig, dt.tmp$new)
    counts <- merge(counts, dat, by.x='symbol', by.y='phenotype_id')
    setcolorder(counts, c('#chr','start','end','symbol',dt.tmp$new))
    setnames(counts, 'symbol','phenotype_id')
    setkey(counts, '#chr', 'start', 'end')
    fwrite(counts[, .SD, .SDcols=c('#chr','start','end','phenotype_id',hbcc_samples)], file=paste0('QTL-pseudobulk-counts/', mode, '-HBCC-', celltype, '-counts.bed'), quote=F, row.names=F, col.names=T, sep='\t')
    fwrite(counts[, .SD, .SDcols=c('#chr','start','end','phenotype_id',nabec_samples)], file=paste0('QTL-pseudobulk-counts/', mode, '-NABEC-', celltype, '-counts.bed'), quote=F, row.names=F, col.names=T, sep='\t')
}



mode <- 'atac'
for(celltype in celltypes) {
    countsfile <- paste0('/data/CARD_singlecell/brain_atlas_subtype/output/', mode, '/', celltype, '/cpm_log_pseudobulk_model_counts.csv')
    counts <- fread(countsfile)
    counts[, c('chr','start','end') := tstrsplit(peaks, split='[-:]')]
    counts[, start := as.numeric(start)]
    counts[, end := as.numeric(end)]
    counts[, start := ceiling((start+end)/2)]
    counts[, end := start+1]
    samplenames <- grep('-ARC$', colnames(counts), value=TRUE)
    # add 'HBCC_' to samples that begin with number
    dt.tmp <- data.table(orig=samplenames)
    dt.tmp[orig %like% '^[0-9]', v2 := paste0('HBCC_', orig)]
    dt.tmp[orig %like% '^[0-9]', cohort := 'HBCC']
    
    dt.tmp[orig %like% '^[SUK]', v2 := orig]
    dt.tmp[orig %like% '^[SUK]', cohort := 'NABEC']
    dt.tmp[, new := gsub('-ARC$', '', v2)]
    nabec_samples <- dt.tmp[cohort=='NABEC', new]
    hbcc_samples <- dt.tmp[cohort=='HBCC', new]

    setnames(counts, dt.tmp$orig, dt.tmp$new)
    setcolorder(counts, c('chr','start','end','peaks',dt.tmp$new))
    setnames(counts, 'chr','#chr')
    setnames(counts, 'peaks','phenotype_id')
    setkey(counts, '#chr', 'start', 'end')
    fwrite(counts[, .SD, .SDcols=c('#chr','start','end','phenotype_id',hbcc_samples)], file=paste0('QTL-pseudobulk-counts/', mode, '-HBCC-', celltype, '-counts.bed'), quote=F, row.names=F, col.names=T, sep='\t')
    fwrite(counts[, .SD, .SDcols=c('#chr','start','end','phenotype_id',nabec_samples)], file=paste0('QTL-pseudobulk-counts/', mode, '-NABEC-', celltype, '-counts.bed'), quote=F, row.names=F, col.names=T, sep='\t')
}

# Create array-params.tsv
combos <- CJ('celltype'=celltypes, 'mode'=c('atac','rna'), cohort=c('NABEC','HBCC'))
fwrite(combos, file='data/array-params.tsv', quote=F, row.names=F, col.names=F, sep='\t')