library(haven)
library(ggplot2)
library(lme4)
library(glmmTMB)
library(Seurat)



dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
datadir <- paste0(dir,"output/nb_partitioning_rep_separate/data/")
outdir <- paste0(dir,"output/Variance/Clusters/nb_partitioning_rep_separate/data/")
dir.create(outdir, recursive = TRUE)



##### Read in seurat with genes #####
seurat_files <- list.files(datadir, pattern = "1pct.rds")

seurat_list <- lapply(seurat_files, function(x){
	readRDS(paste0(datadir,x))
})



### Check that Location and Time has indeed been updated
lapply(seurat_list, function(x){
	print(unique(x@meta.data$Time)) 
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
names(seurat_sub) <- gsub("_seurat_1pct.rds", "", seurat_files)

lapply(seurat_sub, function(x) lapply(x, dim)) ## max rows = 15077


lapply(names(seurat_sub), function(rep){
	lapply(names(seurat_sub[[rep]]), function(cluster){
		saveRDS(seurat_sub[[rep]][[cluster]], paste0(outdir, rep, "_cluster", cluster,"_SCT_seurat_1pct.rds"))
	})
})









