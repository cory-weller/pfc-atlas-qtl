args = [
    '--traw','genotypes/HBCC_autosomes.traw',
    '--counts','/data/CARD_singlecell/PFC_atlas/data/celltypes/Astro/pseudobulk_consensus_atac.csv',
    '--counts-index','sample_id',
    '--outdir', 'QTL-test/HBCC-Astro-atac',
    '--window', '1000000',
    '--covariates', 'genotypes/HBCC_tqtl_covariates.csv',
    '--covariates-index', 'SampleID',
    '--mode', 'atac'
]
