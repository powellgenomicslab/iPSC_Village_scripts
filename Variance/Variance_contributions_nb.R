library(tidyverse)
library(hier.part.negbinom)
library(Seurat)


##### Bring in variables #####
### Bring in arguments
args <- commandArgs(trailingOnly = TRUE)
outdir <- paste0(args[1])
number <- as.numeric(args[2])




dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"




##### Read in seurat with genes #####
seurat <- readRDS(paste0(dir,"output/Distribution_tests/seurat_integrated_all_times_clustered_1pct_expressing.rds"))
seurat_sub <- subset(seurat, subset = Time != "Thawed Village Day 0")
seurat_sub <- subset(seurat_sub, subset = Time != "Thawed Village Day 7")
gene <- rownames(seurat[["SCT"]]@counts)[number]

seurat_sub@meta.data$Location <- gsub("_.+", "", seurat_sub@meta.data$Location)

seurat_sub_list  <- list()

for (location in unique(seurat_sub@meta.data$Location)){
	seurat_sub_list[[location]] <- subset(seurat_sub, subset = Location == location)
}


hier_models <- list()
for(x in names(seurat_sub_list)){
	df <- data.frame("Expression" = data.frame(seurat_sub_list[[x]][["SCT"]]@counts[gene,]), "Line" = seurat_sub_list[[x]]@meta.data$Final_Assignment, "Village" = ifelse(seurat_sub_list[[x]]@meta.data$Time == "Baseline", 0, 1), "Replicate" = gsub("[A-Z][a-z]+", "", seurat_sub_list[[x]]@meta.data$MULTI_ID), "Phase" = seurat_sub_list[[x]]@meta.data$phases,  "MT" = seurat_sub_list[[x]]@meta.data$percent.mt, "RB" = seurat_sub_list[[x]]@meta.data$percent.rb)
	# df <- data.frame("Expression" = data.frame(seurat_sub_list[[x]][["SCT"]]@counts[gene,]), "Line" = seurat_sub_list[[x]]@meta.data$Final_Assignment, "Village" = ifelse(seurat_sub_list[[x]]@meta.data$Time == "Baseline", 0, 1), "Replicate" = gsub("[A-Z][a-z]+", "", seurat_sub_list[[x]]@meta.data$MULTI_ID), "Phase" = seurat_sub_list[[x]]@meta.data$phases,  "MT" = seurat_sub_list[[x]]@meta.data$percent.mt, "RB" = seurat_sub_list[[x]]@meta.data$percent.rb)
	colnames(df)[1] <- "Expression"
	hier_models[[x]][[gene]] <- hier.part.negbinom::hier.part(y = df$Expression, xcan = df[,2:ncol(df)], family = "negbinom",gof = "logLik")
}


saveRDS(hier_models, paste0(outdir,gene,"_nb_hier_models.rds"))