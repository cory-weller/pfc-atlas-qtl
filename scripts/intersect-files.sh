#!/usr/bin/env bash

module load R/4.3

parallel -j 1 Rscript scripts/intersect-files.R {} ::: HBCC NABEC ::: rna atac ::: ExN InN MG VC OPC Oligo