library(Seurat)
library(tivdyverse)
library(data.table)


dir.create("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Kilpinen_eQTLs/")



seurat <- readRDS("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/preQC/time-integrated_filtered_seurat_1pct_expressing.rds")


### Make DF for modeling ###
df_hier_unscale <- data.frame("Passage" = as.factor(gsub("Village_", "", seurat@meta.data$Pool)), "Line" = seurat@meta.data$AtLeastHalfSinglet_Individual_Assignment)
df_hier_unscale$Barcode <- colnames(seurat)

fwrite(df_hier_unscale, paste0("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/multi-passage/cell_meta.tsv"), sep = "\t")

