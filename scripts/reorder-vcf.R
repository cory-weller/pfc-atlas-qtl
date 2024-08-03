    #!/usr/bin/env Rscript
    
    library(data.table)
    args <- commandArgs(trailingOnly=TRUE)
    sample <- args[1]

    # Import genotypes vcf, fingerprints file, and set keys
    vcf.genotypes <- fread(paste0(sample,'.vcf'))
    vcf.fingerprints <- fread(paste0(sample, '.10Xrna.vcf'), skip='#CHROM')
    setnames(vcf.genotypes, '#CHROM', 'CHROM')
    setnames(vcf.fingerprints, '#CHROM', 'CHROM')
    vcf_cols <- c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT')
    setnames(vcf.fingerprints, c(vcf_cols, 'GT.fingerprint'))
    setnames(vcf.genotypes, c(vcf_cols, 'GT.genotype'))

    setkey(vcf.genotypes, 'CHROM','POS')
    setkey(vcf.fingerprints, 'CHROM','POS')

    # Merge to find overlapping loci
    dat.merge <- merge(vcf.genotypes, vcf.fingerprints)
    setnames(dat.merge, 'REF.y','REF')
    setnames(dat.merge, 'ALT.y','ALT')

    # Fix loci where ref/alt allele has been swapped
    dat.merge[,GT_pre := GT.genotype]
    dat.merge[REF != REF.x & GT_pre == '0/0', GT_post := '1/1']
    dat.merge[REF != REF.x & GT_pre == '1/1', GT_post := '0/0']
    dat.merge[is.na(GT_post), GT_post := GT_pre]
    dat.merge[, ID.x := paste(CHROM, POS, REF, ALT, sep=':')]
    setnames(dat.merge, 'ID.x', 'ID')
    setnames(dat.merge, 'QUAL.x', 'QUAL')
    setnames(dat.merge, 'FILTER.x', 'FILTER')
    setnames(dat.merge, 'INFO.x', 'INFO')
    setnames(dat.merge, 'FORMAT.x', 'FORMAT')

    # Structure output VCF
    vcf.out <- dat.merge[, .SD, .SDcols=c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','GT_post')]
    setnames(vcf.out, 'GT_post', sample)
    setnames(vcf.out, 'CHROM', '#CHROM')

    # Output
    fwrite(vcf.out, file=paste0(sample, '.genotype-fingerprints-noheader.vcf'), quote=F, row.names=F, col.names=T, sep='\t')