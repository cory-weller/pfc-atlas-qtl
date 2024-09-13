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
    setnames(dat, c('chr','TSS_start','TSS_start_plus_1','symbol'))
    dat <- dat[chr %in% c(paste0('chr', 1:22), 'chrX')]
    setkey(dat, chr, TSS_start, TSS_start_plus_1)
    # Remove duplicate gene entries
    #       chr     start       end         symbol
    #     <char>     <int>     <int>         <char>
    #  1:   chrX 102599512 102714671 ARMCX5-GPRASP2
    #  2:   chrX 102712495 102753530 ARMCX5-GPRASP2
    #  3:   chr3  50350695  50358460       CYB561D2
    #  4:   chr3  50365334  50368197       CYB561D2
    #  5:  chr22  24556007  24629005           GGT1
    #  6:  chr22  24594811  24629005           GGT1
    #  7:  chr15  28701954  28738384        GOLGA8M
    #  8:  chr15  28719377  28738431        GOLGA8M
    #  9:  chr10  14838160  14847018         HSPA14
    # 10:  chr10  14838306  14871741         HSPA14
    # 11:   chr2 241970683 241977276      LINC01238
    # 12:   chr2 242087351 242088457      LINC01238
    # 13:   chr9 105993310 106740875      LINC01505
    # 14:   chr9 106664754 106702662      LINC01505
    # 15:   chr5 139273752 139331671          MATR3
    # 16:   chr5 139293674 139331677          MATR3
    # 17:   chr1 235328570 235448952           TBCE
    # 18:   chr1 235367360 235452443           TBCE
    # 19:   chrX 103918896 103966712        TMSB15B
    # 20:   chrX 104063871 104076212        TMSB15B
    dat <- dat[!duplicated(symbol)]
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
    countsfile <- paste0('/data/CARD_singlecell/brain_atlas_subtype/output/', mode, '/', celltype, '/cpm_log_pseudobulk_sum_counts.csv')
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
    counts <- merge(counts, dat, by='symbol')
    setcolorder(counts, c('chr','TSS_start','TSS_start_plus_1','symbol',dt.tmp$new))
    setnames(counts, 'chr','#chr')
    setnames(counts, 'TSS_start','start')
    setnames(counts, 'TSS_start_plus_1','end')
    setnames(counts, 'symbol','phenotype_id')
    setkey(counts, '#chr', 'start', 'end')
    fwrite(counts[, .SD, .SDcols=c('#chr','start','end','phenotype_id',hbcc_samples)], file=paste0('data/', mode, '/HBCC-', celltype, '-counts.bed'), quote=F, row.names=F, col.names=T, sep='\t')
    fwrite(counts[, .SD, .SDcols=c('#chr','start','end','phenotype_id',nabec_samples)], file=paste0('data/', mode, '/NABEC-', celltype, '-counts.bed'), quote=F, row.names=F, col.names=T, sep='\t')
}



mode <- 'atac'
for(celltype in celltypes) {
    countsfile <- paste0('/data/CARD_singlecell/brain_atlas_subtype/output/', mode, '/', celltype, '/cpm_log_pseudobulk_sum_counts.csv')
    counts <- fread(countsfile)
    counts[, c('chr','start','end') := tstrsplit(regions, split='[-:]')]
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
    setcolorder(counts, c('chr','start','end','regions',dt.tmp$new))
    setnames(counts, 'chr','#chr')
    setnames(counts, 'regions','phenotype_id')
    setkey(counts, '#chr', 'start', 'end')
    fwrite(counts[, .SD, .SDcols=c('#chr','start','end','phenotype_id',hbcc_samples)], file=paste0('data/', mode, '/HBCC-', celltype, '-counts.bed'), quote=F, row.names=F, col.names=T, sep='\t')
    fwrite(counts[, .SD, .SDcols=c('#chr','start','end','phenotype_id',nabec_samples)], file=paste0('data/', mode, '/NABEC-', celltype, '-counts.bed'), quote=F, row.names=F, col.names=T, sep='\t')
}


combos <- CJ('celltype'=celltypes, 'mode'=c('atac','rna'), cohort=c('NABEC','HBCC'))
fwrite(combos, file='data/array-params.tsv', quote=F, row.names=F, col.names=F, sep='\t')