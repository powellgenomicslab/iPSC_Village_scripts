#library(lme4)
library(tidyverse)
library(haven)
library(ggplot2)
library(lme4)
library(glmmTMB)
library(Seurat)



dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
datadir <- paste0(dir,"output/nb_partitioning_cryo/data/")
outdir <- paste0(dir,"output/Variance/Clusters/nb_variance_partitioning_w_cryo/data/")
dir.create(outdir, recursive = TRUE)



##### Read in seurat with genes #####
seurat_files <- list.files(datadir, pattern = "1pct.rds")

seurat_list <- lapply(seurat_files, function(x){
	readRDS(paste0(datadir,x))
})

### Check that Location and Time have indeed been updated
lapply(seurat_list, function(x){
	print(unique(x@meta.data$Location)) # <- gsub("_Baseline", "", seurat@meta.data$Location) %>% gsub("_Village.+", "", .) %>% gsub("Thawed", "Cryopreserved",.)
})

lapply(seurat_list, function(x){
	print(unique(x@meta.data$Time)) ## <- gsub("Thawed Village Day 0", "Baseline", x@meta.data$Time) %>% gsub("Thawed Village Day 7", "Village", .)
})

### Fix Time to all be Baseline and Village (won't have mattered for how coded the negative binomial assessment before but would be good for consistency)
seurat_list <- lapply(seurat_list, function(x){
	x@meta.data$Time <- gsub("Village Day 4", "Village", x@meta.data$Time)
	return(x)
})


### Subset the Clusters ###
seurat_list_cluster <- lapply(seurat_list, function(x){
	tmp <- list()
	for (cluster in levels(Idents(x))){
		tmp[[cluster]] <- subset(x, idents = cluster)
	}
	return(tmp)
})


seurat_sub <- lapply(seurat_list_cluster, function(x){
	lapply(x, function(y){
		subset(y, features = rownames(y)[which(rowSums(y[["SCT"]]@counts > 0)/ncol(y[["SCT"]]@counts) >= 0.01)])
	})
})
names(seurat_sub) <- gsub("_SCT_seurat_1pct.rds", "", seurat_files)

lapply(seurat_sub, function(x) lapply(x, dim)) ## max rows = 15087


lapply(names(seurat_sub), function(location){
	lapply(names(seurat_sub[[location]]), function(cluster){
		saveRDS(seurat_sub[[location]][[cluster]], paste0(outdir, location, "_cluster", cluster,"_SCT_seurat_1pct.rds"))
	})
})


