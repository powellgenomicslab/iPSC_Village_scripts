library(haven)
library(ggplot2)
library(lme4)
library(glmmTMB)
library(Seurat)
library(data.table)
library(tidyverse)



dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
outdir <- paste0(dir,"output/Variance/RNAvelocity/nb_partitioning_rep_separate/data/")
dir.create(outdir, recursive = TRUE)


##### Read in Data #####
### seurat object ###
seurat <- readRDS(paste0(dir,"output/Distribution_tests/seurat_integrated_all_times_clustered_1pct_expressing.rds"))

### velocity metadata ###
velo_meta <- fread("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/scVelo/velocyto/scvelo_combat_corrected_all2/metadata.csv", sep = ",")
rownames(velo_meta) <- velo_meta$V1
velo_meta$V1 <- NULL

velo_meta_sub <- velo_meta[,c("Location", "n_unspliced_counts", "latent_time")]
rownames(velo_meta_sub) <- rownames(velo_meta)

seurat <- AddMetaData(seurat, velo_meta_sub)


seurat@meta.data$MULTI_ID <- ifelse(seurat@meta.data$Time %in% c("Thawed Village Day 0", "Thawed Village Day 7"), gsub("Sydney", "Sydney_cryopreserved", seurat@meta.data$MULTI_ID), seurat@meta.data$MULTI_ID) 

##### Subset for cells with latent time and #####
seurat_noNA <- subset(seurat, subset = latent_time >= 0)


seurat_noNA@meta.data$Quintile  <- with(seurat_noNA@meta.data, factor(
                            		findInterval(latent_time, c(-Inf,
                               		quantile(latent_time, probs=c(0.2, 0.4, 0.6, 0.8)), Inf)), 
                            		labels=c("Q1","Q2","Q3","Q4","Q5")))

seurat_noNA@meta.data$Rep_Quintile <- paste0(seurat_noNA@meta.data$MULTI_ID, "_", seurat_noNA@meta.data$Quintile)


seurat_sub <- list()

for (rep in unique(seurat_noNA@meta.data$Rep_Quintile)){
	seurat_sub[[rep]] <- subset(seurat_noNA, subset = Rep_Quintile == rep)
}



seurat_sub <- lapply(seurat_sub, function(y){
	subset(y, features = rownames(y)[which(rowSums(y[["SCT"]]@counts > 0)/ncol(y[["SCT"]]@counts) >= 0.01)])
})



seurat_sub_SCT <- lapply(seurat_sub, function(x){
	SCTransform(x, verbose = TRUE, vars.to.regress = c("scores.G1", "scores.S", "scores.G2M", "percent.mt", "percent.rb"), return.only.var.genes = FALSE)
})


### Update the cryopreserved Time names
seurat_sub_SCT <- lapply(seurat_sub_SCT, function(x){
	x@meta.data$Time <- gsub("Thawed Village Day 0", "Baseline", x@meta.data$Time) %>% gsub("Thawed Village Day 7", "Village", .)
	return(x)
})

df <- data.frame(matrix(nrow = length(seurat_sub_SCT), ncol = 2))
colnames(df) <- c("rows", "columns")
rownames(df) <- names(seurat_sub_SCT)

for (name in names(seurat_sub_SCT)){
	df[name,] <- dim(seurat_sub_SCT[[name]])
} 

max(df$rows) ## max rows = 15106


lapply(names(seurat_sub_SCT), function(rep){
	saveRDS(seurat_sub_SCT[[rep]], paste0(outdir,rep,"_seurat_1pct.rds"))
})









