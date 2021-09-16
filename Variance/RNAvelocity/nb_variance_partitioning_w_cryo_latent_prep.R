library(tidyverse)
library(haven)
library(ggplot2)
library(lme4)
library(glmmTMB)
library(Seurat)
library(data.table)




dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
outdir <- paste0(dir,"output/Variance/RNAvelocity/nb_partitioning_w_latent/data/")
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


seurat_sub <- list()
for (location in unique(seurat_noNA@meta.data$Location)){
	seurat_sub[[location]] <- subset(seurat_noNA, subset = Location == location)
}


seurat_sub <- lapply(seurat_sub, function(x){
	subset(x, features = rownames(x)[which(rowSums(x[["SCT"]]@counts > 0)/ncol(x[["SCT"]]@counts) >= 0.01)])
})


seurat_sub_SCT <- lapply(seurat_sub, function(x){
	SCTransform(x, verbose = TRUE, vars.to.regress = c("scores.G1", "scores.S", "scores.G2M", "percent.mt", "percent.rb"), return.only.var.genes = FALSE)
})


df <- data.frame(matrix(nrow = length(seurat_sub_SCT), ncol = 2))
colnames(df) <- c("rows", "columns")
rownames(df) <- names(seurat_sub_SCT)

for (name in names(seurat_sub_SCT)){
	df[name,] <- dim(seurat_sub_SCT[[name]])
} 

max(df$rows) ## max rows = 15129


### Update the cryopreserved Time names
seurat_sub_SCT <- lapply(seurat_sub_SCT, function(x){
	x@meta.data$Time <- gsub("Thawed Village Day 0", "Baseline", x@meta.data$Time) %>% gsub("Thawed Village Day 7", "Village", .)
	return(x)
})



lapply(names(seurat_sub_SCT), function(location){
	saveRDS(seurat_sub_SCT[[location]], paste0(outdir, location, "_SCT_seurat_1pct.rds"))
})


