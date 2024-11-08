args = [
    "--mode","atac",
    "--celltype","Astro",
    "--covariates", "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","Sex","Age",
    "--cohort", "NABEC",
    "--qtlmethod", "nominal",
    "--interaction", "None",
    "--outdir", "QTL-output",
    "--window", '1000000',
    "--pseudobulkmethod","sum"
]