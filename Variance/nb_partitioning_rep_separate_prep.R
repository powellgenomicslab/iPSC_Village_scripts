library(haven)
library(ggplot2)
library(lme4)
library(glmmTMB)
library(Seurat)



dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
outdir <- paste0(dir,"output/nb_partitioning_rep_separate/data/")
dir.create(outdir, recursive = TRUE)



seurat <- readRDS(paste0(dir,"output/Distribution_tests/seurat_integrated_all_times_clustered_1pct_expressing.rds"))


seurat@meta.data$MULTI_ID <- ifelse(seurat@meta.data$Time %in% c("Thawed Village Day 0", "Thawed Village Day 7"), gsub("Sydney", "Sydney_cryopreserved", seurat@meta.data$MULTI_ID), seurat@meta.data$MULTI_ID) 


seurat_sub <- list()

for (rep in unique(seurat@meta.data$MULTI_ID)){
	seurat_sub[[rep]] <- subset(seurat, subset = MULTI_ID == rep)
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


lapply(seurat_sub_SCT, dim) ## max rows = 15128


lapply(names(seurat_sub_SCT), function(rep){
	saveRDS(seurat_sub_SCT[[rep]], paste0(outdir,rep,"_seurat_1pct.rds"))
})









