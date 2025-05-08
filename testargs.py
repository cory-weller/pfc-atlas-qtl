args = [
    '--traw','genotypes/HBCC_autosomes.traw',
    '--counts','/data/CARD_singlecell/PFC_atlas/data/celltypes/Oligo/pseudobulk_rna.csv',
    '--counts-index','SampleID',
    '--outdir', 'QTL-test/HBCC-Oligo-ATAC',
    '--window', '1000000',
    '--covariates', 'genotypes/HBCC_tqtl_covariates.csv',
    '--covariates-index', 'SampleID',
    '--mode', 'rna'
]
