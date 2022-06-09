library(tidyverse)
library(haven)
library(ggplot2)
library(lme4)
library(glmmTMB)
library(Seurat)
library(data.table)




dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
outdir <- paste0(dir,"output/Variance/RNAvelocity/post_review/data/")
dir.create(outdir, recursive = TRUE)



##### Read in Data #####
### seurat object ###
seurat <- readRDS(paste0(dir,"output/Distribution_tests/seurat_integrated_all_times_clustered_1pct_expressing.rds"))

### velocity metadata ###
velo_meta <- fread("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/scVelo/velocyto/scvelo_combat_corrected_all2/metadata.csv", sep = ",")
rownames(velo_meta) <- velo_meta$V1
velo_meta$V1 <- NULL

velo_meta_sub <- velo_meta[,c("n_unspliced_counts", "latent_time")]
rownames(velo_meta_sub) <- rownames(velo_meta)

seurat <- AddMetaData(seurat, velo_meta_sub)


seurat@meta.data$Location <- gsub("_Baseline", "", seurat@meta.data$Location) %>% gsub("_Village.+", "", .) %>% gsub("Thawed", "Cryopreserved",.)

seurat_noNA <- subset(seurat, subset = latent_time >= 0)


seurat_noNA@meta.data$Location <- gsub("_Baseline", "", seurat_noNA@meta.data$Location) %>% gsub("_Village.+", "", .) %>% gsub("Thawed", "Cryopreserved",.)
seurat_noNA@meta.data$Time <- gsub("Thawed Village Day 0", "Baseline", seurat_noNA@meta.data$Time) %>% gsub("Thawed Village Day 7", "Village", .) %>% gsub(" Day 4", "", .)
seurat_noNA@meta.data$Cryopreserved <-ifelse(seurat_noNA@meta.data$Location == "Sydney_Cryopreserved", "Cryopreserved", "Fresh")
seurat_noNA@meta.data$Location <- gsub("_Cryopreserved", "", seurat_noNA@meta.data$Location)


saveRDS(seurat_noNA, paste0(outdir, "seurat_integrated_all_times_clustered_1pct_expressing_pseudotime.rds"))