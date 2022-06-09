library(Seurat)
liberary(tidyverse)
library(dplyr)
library(data.table)


##### Read in seurat with genes #####
seurat <- readRDS("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Distribution_tests/seurat_integrated_all_times_clustered_1pct_expressing.rds")

seurat@meta.data$Location <- gsub("_Baseline", "", seurat@meta.data$Location) %>% gsub("_Village.+", "", .) %>% gsub("Thawed", "Cryopreserved",.)
seurat@meta.data$Time <- gsub("Thawed Village Day 0", "Baseline", seurat@meta.data$Time) %>% gsub("Thawed Village Day 7", "Village", .) %>% gsub(" Day 4", "", .)

seurat <- subset(seurat, subset = Location != "Sydney_Cryopreserved")

seurat <- subset(seurat, features = rownames(seurat)[which(rowSums(seurat[["SCT"]]@counts > 0)/ncol(seurat[["SCT"]]@counts) >= 0.01)])

seurat <- SCTransform(seurat, verbose = TRUE, vars.to.regress = c("scores.G1", "scores.S", "scores.G2M", "percent.mt", "percent.rb"), return.only.var.genes = FALSE)

saveRDS(seurat, "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/variance_partition_post_review/seurat_integrated_noncryo_1pct_expressing.rds")

fwrite(data.table(Gene = rownames(seurat)), "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/variance_partition_post_review/seurat_integrated_noncryo_1pct_expressing_genes.tsv", sep = "\t")