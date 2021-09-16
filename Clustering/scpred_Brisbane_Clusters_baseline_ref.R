
library(Seurat)
library(tidyverse)
library(ggplot2)
library(scPred)

dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
outdir <- paste0(dir,"output/scPred_Cluster_Assignments_Baseline_ref/")
dir.create(outdir)


##### Read in Seurat Objects #####
seurat_brisbane_baseline <- readRDS(paste0(dir,"output/Brisbane_Clustering_Testing/Make_Final_Brisbane_Baseline_Seurat/Brisbane_Seurat_SCT_CellCycle_cov_Indiv_integration_w_meta.rds"))

file_list <- list.files(paste0(dir,"output/Seurat_w_Hash_freeze_thaw/CellCycle"), pattern = ".rds")
file_list <- file_list[!grepl("Brisbane_Baseline", file_list)]
file_list <- file_list[!grepl("seurat_filtered_cell_cycle.rds", file_list)]
seurat_list <- lapply(file_list, function(x){
	readRDS(paste0(dir,"output/Seurat_w_Hash_freeze_thaw/CellCycle/",x))
})
names(seurat_list) <- gsub("_seurat.rds", "", file_list) %>% gsub(".rds","",.)
names <- gsub("_seurat.rds", "", file_list) %>% gsub(".rds","",.)


## Add location to the seurat metadata ##
temp <- list()
seurat_list <- lapply(names(seurat_list), function(x){
	temp[[x]] <- list()
	print(x)
	seurat_list[[x]]$Location <- gsub("\\d","", seurat_list[[x]]$Site_rep)
	temp[[x]] <- lapply(unique(seurat_list[[x]]$Final_Assignment), function(y){
		subset(seurat_list[[x]], subset = Final_Assignment == y)
	})
	names(temp[[x]]) <- unique(seurat_list[[x]]$Final_Assignment)
	return(temp[[x]])
})

names(seurat_list) <- gsub("_seurat.rds", "", file_list) %>% gsub(".rds","",.)



##### Train the model #####
reference <- getFeatureSpace(seurat_brisbane_baseline, "Clusters_Res_0.09")
reference <- trainModel(reference)

get_probabilities(reference) %>% head()

get_scpred(reference)

plot <- plot_probabilities(reference)

ggsave(plot, filename = paste0(outdir, "scPred_Performance_caret.png"))

#### To retrain #####
reference <- trainModel(reference, model = "mda")


get_probabilities(reference) %>% head()

get_scpred(reference)

plot_mda <- plot_probabilities(reference)

ggsave(plot_mda, filename = paste0(outdir, "scPred_Performance_mda.png"))


### Check prediction for Brisbane Village Day 4 ###
query_list <- lapply(seurat_list, function(x){
	temp_list <- lapply(x, function(y) {
		SCTransform(y, verbose = TRUE, vars.to.regress = c("scores.G1", "scores.S", "scores.G2M"))
	})

	seurat_features <- SelectIntegrationFeatures(object.list = temp_list, nfeatures = 3000)
	temp_list <- PrepSCTIntegration(object.list = temp_list, anchor.features = seurat_features)
	seurat_anchors <- FindIntegrationAnchors(object.list = temp_list, normalization.method = "SCT", anchor.features = seurat_features)

	integrated <- IntegrateData(anchorset = seurat_anchors, normalization.method = "SCT")
	return(integrated)
})

query_list <- lapply(query_list, function(x){
	scPredict(x, reference)
})
names(query_list) <- names
 

PCA_plot_list <- lapply(query_list, function(x){
	DimPlot(x, group.by = "scpred_prediction", reduction = "scpred")
})

lapply(names(PCA_plot_list), function(x){
	ggsave(PCA_plot_list[[x]], filename = paste0(outdir, x,"_predicted_PCA.png"))
})

query_list <- lapply(query_list, function(x){
	RunUMAP(x, reduction = "scpred", dims = 1:30)
})

UMAP_plot_list <- lapply(query_list, function(x){
	DimPlot(x, group.by = "scpred_prediction", label = TRUE, repel = TRUE)
})

lapply(names(UMAP_plot_list), function(x){
	ggsave(UMAP_plot_list[[x]], filename = paste0(outdir, x, "_predicted_UMAP.png"))
})

UMAP_plot_Assignment_list <- lapply(query_list, function(x){
	DimPlot(x, group.by = "scpred_prediction", label = TRUE, repel = TRUE, split.by = "Final_Assignment")
})

lapply(names(UMAP_plot_Assignment_list), function(x){
	ggsave(UMAP_plot_Assignment_list[[x]], filename = paste0(outdir, x, "_predicted_UMAP_individual.png"), width = 21, height = 7)
})


query_list <- lapply(query_list, function(x){
	x %>%
	RunPCA() %>%
	RunUMAP(., reduction = "pca", dims = 1:30, reduction.name = "pca_umap")
}) 

UMAP_plot_Location_list <- lapply(query_list, function(x){
	DimPlot(x, group.by = "scpred_prediction", label = TRUE, repel = TRUE, split.by = "Final_Assignment", reduction = "pca_umap")
})

lapply(names(UMAP_plot_Location_list), function(x){
	ggsave(UMAP_plot_Location_list[[x]], filename = paste0(outdir, x, "_predicted_PCA_UMAP.png"), width = 21, height = 7)
})


##### Save the seurat objecst #####
saveRDS(reference, paste0(outdir, "Brisbane_Baseline_post_scpred.rds"))
# reference <- readRDS(paste0(outdir, "Brisbane_Baseline_post_scpred.rds"))

lapply(names(query_list), function(x){
	saveRDS(query_list[[x]], paste0(outdir, x, "_post_scpred.rds"))
})


query_list <- list()
query_list <- lapply(names, function(x){
	readRDS(paste0(outdir, x, "_post_scpred.rds"))
})

##### Combine all data together to investigate similarities and differences
combined_assigned <- merge(reference, query_list)


##### Look at table of number of each cell type in each dataset #####
lapply(query_list, function(x){
	prop.table(table(x$scpred_prediction))
})


##### Check to see if any other covariates are confounding the PCA #####
pca_cov_outdir <- paste0(outdir,"scpred_reduc_covariates/")
dir.create(pca_cov_outdir)

scpred_plot_cov_list <- lapply(query_list, function(x){
temp <- list()
	for (cov in c("phases", "Final_Assignment","Pool","Time" )){
		temp[[cov]] <- DimPlot(x, group.by = cov, reduction = "scpred")
	}
	for (cov in c( "percent.rb", "percent.mt", "nCount_RNA", "nFeature_RNA"))
	{
		temp[[cov]] <- FeaturePlot(x, reduction = "scpred", features = cov)
	}
	return(temp)
})

lapply(names(PCA_plot_cov_list), function(x){
	lapply(names(PCA_plot_cov_list[[x]]), function(cov){
		ggsave(PCA_plot_cov_list[[x]][[cov]], filename = paste0(pca_cov_outdir, x,"_",cov,"_predicted_PCA.png"))
	})
})


pca_cov_outdir <- paste0(outdir,"scpred_reduc_umap_covariates/")
dir.create(pca_cov_outdir)

scpred_umap_cov_list <- lapply(query_list, function(x){
temp <- list()
	for (cov in c("phases", "Final_Assignment","Pool","Time" )){
		temp[[cov]] <- DimPlot(x, group.by = cov, reduction = "umap")
	}
	for (cov in c( "percent.rb", "percent.mt", "nCount_RNA", "nFeature_RNA"))
	{
		temp[[cov]] <- FeaturePlot(x, reduction = "umap", features = cov)
	}
	return(temp)
})

lapply(names(scpred_umap_cov_list), function(x){
	lapply(names(scpred_umap_cov_list[[x]]), function(cov){
		ggsave(scpred_umap_cov_list[[x]][[cov]], filename = paste0(pca_cov_outdir, x,"_",cov,"_predicted_UMAP.png"))
	})
})




pca_cov_outdir <- paste0(outdir,"PCA_covariates/")
dir.create(pca_cov_outdir)

PCA_plot_cov_list <- lapply(query_list, function(x){
temp <- list()
	for (cov in c("phases", "Final_Assignment","Pool","Time", "scpred_prediction" )){
		temp[[cov]] <- DimPlot(x, group.by = cov, reduction = "pca")
	}
	for (cov in c( "percent.rb", "percent.mt", "nCount_RNA", "nFeature_RNA"))
	{
		temp[[cov]] <- FeaturePlot(x, reduction = "pca", features = cov)
	}
	return(temp)
})

lapply(names(PCA_plot_cov_list), function(x){
	lapply(names(PCA_plot_cov_list[[x]]), function(cov){
		ggsave(PCA_plot_cov_list[[x]][[cov]], filename = paste0(pca_cov_outdir, x,"_",cov,"_predicted_PCA.png"))
	})
})
