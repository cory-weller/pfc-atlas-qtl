args = [
    "--mode","rna",
    "--celltype","Astro",
    "--covariates", "PC1","PC2","PC3","PC4","PC5","Sex","Age",
    "--modality-covariates", "PC1","PC2","PC3","PC4","PC5",
    "--cohort", "NABEC",
    "--qtlmethod", "nominal",
    "--interaction", "None",
    "--outdir", "QTL-output",
    "--window", '1000000',
    "--pseudobulkmethod","sum",
    "--chr","chr1"
]
