library(Seurat)
liberary(tidyverse)
library(dplyr)
library(data.table)


##### Read in seurat with genes #####
seurat <- readRDS("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Distribution_tests/seurat_integrated_all_times_clustered_1pct_expressing.rds")

dir.create("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/variance_partitioning_all_cells/", recursive = TRUE)

fwrite(data.table(Gene = rownames(seurat)), "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/variance_partitioning_all_cells/seurat_integrated_Sydney_1pct_expressing_genes.tsv", sep = "\t")
