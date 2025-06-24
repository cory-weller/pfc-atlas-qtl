#!/usr/bin/env bash
#SBATCH --mem 150G
#SBATCH --cpus-per-task 2
#SBATCH --time=2:00:00
#SBATCH --partition quick,norm

ml R/4.3

Rscript scripts/save-signif-parquet.R ${1}


exit


# sbatch parse-parquet.sh QTL/HBCC-Astro-atac
# sbatch parse-parquet.sh QTL/HBCC-Astro-rna
# sbatch parse-parquet.sh QTL/HBCC-ExN-atac
# sbatch parse-parquet.sh QTL/HBCC-ExN-rna
# sbatch parse-parquet.sh QTL/HBCC-InN-atac
# sbatch parse-parquet.sh QTL/HBCC-InN-rna
# sbatch parse-parquet.sh QTL/HBCC-MG-atac
# sbatch parse-parquet.sh QTL/HBCC-MG-rna
# sbatch parse-parquet.sh QTL/HBCC-OPC-atac
# sbatch parse-parquet.sh QTL/HBCC-OPC-rna
# sbatch parse-parquet.sh QTL/HBCC-Oligo-atac
# sbatch parse-parquet.sh QTL/HBCC-Oligo-rna
# sbatch parse-parquet.sh QTL/HBCC-VC-atac
# sbatch parse-parquet.sh QTL/HBCC-VC-rna
# sbatch parse-parquet.sh QTL/NABEC-Astro-atac
# sbatch parse-parquet.sh QTL/NABEC-Astro-rna
# sbatch parse-parquet.sh QTL/NABEC-ExN-atac
# sbatch parse-parquet.sh QTL/NABEC-ExN-rna
# sbatch parse-parquet.sh QTL/NABEC-InN-atac
# sbatch parse-parquet.sh QTL/NABEC-InN-rna
# sbatch parse-parquet.sh QTL/NABEC-MG-atac
# sbatch parse-parquet.sh QTL/NABEC-MG-rna
# sbatch parse-parquet.sh QTL/NABEC-OPC-atac
# sbatch parse-parquet.sh QTL/NABEC-OPC-rna
# sbatch parse-parquet.sh QTL/NABEC-Oligo-atac
# sbatch parse-parquet.sh QTL/NABEC-Oligo-rna
# sbatch parse-parquet.sh QTL/NABEC-VC-atac
# sbatch parse-parquet.sh QTL/NABEC-VC-rna