library(Seurat)
library(tivdyverse)
library(data.table)


dir.create("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Kilpinen_eQTLs/")



seurat <- readRDS("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/variance_partition_post_review/seurat_integrated_noncryo_1pct_expressing.rds")


### Make DF for modeling ###
df_hier_unscale <- data.frame("Village" = as.factor(ifelse(seurat@meta.data$Time == "Baseline", 0, 1)), "Line" = seurat@meta.data$Final_Assignment, "Replicate" = as.factor(gsub("[A-Z][a-z]+", "", seurat@meta.data$MULTI_ID)), "Site" = seurat$Location)
df_hier_unscale$Barcode <- rownames(df_hier_unscale)

fwrite(df_hier_unscale, paste0("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Kilpinen_eQTLs/cell_meta.tsv"), sep = "\t")

