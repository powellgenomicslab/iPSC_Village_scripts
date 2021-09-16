#library(lme4)
library(tidyverse)
library(haven)
library(ggplot2)
library(lme4)
library(glmmTMB)
library(Seurat)



# ##### Bring in variables #####
# ### Bring in arguments
# args <- commandArgs(trailingOnly = TRUE)
# outdir <- paste0(args[1])
# number <- as.numeric(args[2])




dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
outdir <- paste0(dir,"output/nb_partitioning_cryo/data/")
dir.create(outdir, recursive = TRUE)



##### Read in seurat with genes #####
seurat <- readRDS(paste0(dir,"output/Distribution_tests/seurat_integrated_all_times_clustered_1pct_expressing.rds"))

seurat@meta.data$Location <- gsub("_Baseline", "", seurat@meta.data$Location) %>% gsub("_Village.+", "", .) %>% gsub("Thawed", "Cryopreserved",.)


seurat_sub <- list()
for (location in unique(seurat@meta.data$Location)){
	seurat_sub[[location]] <- subset(seurat, subset = Location == location)
}


seurat_sub <- lapply(seurat_sub, function(x){
	subset(x, features = rownames(x)[which(rowSums(x[["SCT"]]@counts > 0)/ncol(x[["SCT"]]@counts) >= 0.01)])
})


seurat_sub_SCT <- lapply(seurat_sub, function(x){
	SCTransform(x, verbose = TRUE, vars.to.regress = c("scores.G1", "scores.S", "scores.G2M", "percent.mt", "percent.rb"), return.only.var.genes = FALSE)
})


lapply(seurat_sub_SCT, dim) ## max rows = 15128

### Update the cryopreserved Time names
seurat_sub_SCT <- lapply(seurat_sub_SCT, function(x){
	x@meta.data$Time <- gsub("Thawed Village Day 0", "Baseline", x@meta.data$Time) %>% gsub("Thawed Village Day 7", "Village", .)
	return(x)
})



lapply(names(seurat_sub_SCT), function(location){
	saveRDS(seurat_sub_SCT[[location]], paste0(outdir, location, "_SCT_seurat_1pct.rds"))
})


