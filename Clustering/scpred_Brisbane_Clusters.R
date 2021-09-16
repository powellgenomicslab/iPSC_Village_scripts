
library(Seurat)
library(tidyverse)
library(ggplot2)
library(scPred)

dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
outdir <- paste0(dir,"output/scPred_Cluster_Assignments/")
dir.create(outdir)


##### Read in Seurat Objects #####
seurat_brisbane <- readRDS(paste0(dir,"output/Brisbane_Clustering_Testing/Make_Final_Brisbane_Seurat/Brisbane_Seurat_SCT_CellCycle_cov_Indiv_integration_w_meta.rds"))
seurat <- readRDS(paste0(dir,"output/Seurat_w_Hash_freeze_thaw/CellCycle/seurat_filtered_cell_cycle.rds"))

## Subset seurat to include ##
seurat$Location <- gsub("\\d","", seurat$Site_rep)
seurat <- subset(seurat, subset = Location != "Brisbane")


##### Train the model #####
seurat_brisbane_baseline <- subset(seurat_brisbane, subset = Location_Time == "Brisbane_Baseline")
seurat_brisbane_4day <- subset(seurat_brisbane, subset = Location_Time == "Brisbane_Village_Day_4")

seurat_brisbane_baseline <- seurat_brisbane_baseline %>%
	SCTransform(., verbose = TRUE, vars.to.regress = c("scores.G1", "scores.S", "scores.G2M")) %>%
	RunPCA() %>%
	RunUMAP(dims = 1:30)

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
query <- SCTransform(seurat_brisbane_4day, verbose = TRUE, vars.to.regress = c("scores.G1", "scores.S", "scores.G2M"))

query <- scPredict(query, reference)


df <- data.frame("Cell_type" = unique(query$scpred_prediction), "n" = NA, "Sens" = NA, "Spec" = NA)
for (number in 0:4){
	df$n[which(df$Cell_type == number)] <- nrow(query@meta.data[which(query$scpred_prediction == number),])
	df$Sens[which(df$Cell_type == number)] <- nrow(query@meta.data[which(query$scpred_prediction == number & query$Clusters_Res_0.09 == number),])/nrow(query@meta.data[which(query$Clusters_Res_0.09 == number),])
	df$Spec[which(df$Cell_type == number)] <- nrow(query@meta.data[which(query$scpred_prediction != number & query$Clusters_Res_0.09 != number),])/nrow(query@meta.data[which(query$Clusters_Res_0.09 != number),])
}


PCA_plot <- DimPlot(query, group.by = "scpred_prediction", reduction = "scpred")
ggsave(PCA_plot, filename = paste0(outdir, "predicted_PCA.png"))

query <- RunUMAP(query, reduction = "scpred", dims = 1:30)
UMAP_plot <- DimPlot(query, group.by = "scpred_prediction", label = TRUE, repel = TRUE)
ggsave(UMAP_plot, filename = paste0(outdir, "predicted_UMAP.png"))

UMAP_plot_Assignment <- DimPlot(query, group.by = "scpred_prediction", label = TRUE, repel = TRUE, split.by = "Final_Assignment")

ggsave(UMAP_plot_Assignment, filename = paste0(outdir, "predicted_UMAP_individual.png"), width = 21, height = 7)



##### Try with separate times  but integrated Assay #####
DefaultAssay(seurat_brisbane_baseline) <- "integrated"
DefaultAssay(seurat_brisbane_4day) <- "integrated"

seurat_brisbane_baseline <- seurat_brisbane_baseline %>%
	RunPCA() %>%
	RunUMAP(dims = 1:30)

reference_integrated <- getFeatureSpace(seurat_brisbane_baseline, "Clusters_Res_0.09")
reference_integrated <- trainModel(reference_integrated)

get_probabilities(reference_integrated) %>% head()

get_scpred(reference_integrated)

plot_integrated <- plot_probabilities(reference_integrated)

ggsave(plot_integrated, filename = paste0(outdir, "scPred_Performance_caret_integrated.png"))



### Normalize query data ###
query_integrated <- scPredict(seurat_brisbane_4day, reference)


#### Generate specificity, sensitivity ####
df <- data.frame("Cell_type" = unique(query_integrated$scpred_prediction), "n" = NA, "Sens" = NA, "Spec" = NA)
for (number in 0:4){
	df$n[which(df$Cell_type == number)] <- nrow(query_integrated@meta.data[which(query_integrated$scpred_prediction == number),])
	df$Sens[which(df$Cell_type == number)] <- nrow(query_integrated@meta.data[which(query_integrated$scpred_prediction == number & query_integrated$Clusters_Res_0.09 == number),])/nrow(query_integrated@meta.data[which(query_integrated$Clusters_Res_0.09 == number),])
	df$Spec[which(df$Cell_type == number)] <- nrow(query_integrated@meta.data[which(query_integrated$scpred_prediction != number & query_integrated$Clusters_Res_0.09 != number),])/nrow(query_integrated@meta.data[which(query_integrated$Clusters_Res_0.09 != number),])
}


PCA_plot_integrated <- DimPlot(query_integrated, group.by = "scpred_prediction", reduction = "scpred")
ggsave(PCA_plot_integrated, filename = paste0(outdir, "predicted_PCA_integrated.png"))

query <- RunUMAP(query_integrated, reduction = "scpred", dims = 1:30)
UMAP_plot_integrated <- DimPlot(query, group.by = "scpred_prediction", label = TRUE, repel = TRUE)
ggsave(UMAP_plot_integrated, filename = paste0(outdir, "predicted_UMAP_integrated.png"))

UMAP_plot_Location_integrated <- DimPlot(query, group.by = "MULTI_ID", label = TRUE, repel = TRUE)
ggsave(UMAP_plot_Location_integrated, filename = paste0(outdir, "predicted_UMAP_location_integrated.png"))

query <- RunUMAP(query_integrated, reduction = "pca", dims = 1:30)
UMAP_plot <- DimPlot(query, group.by = "scpred_prediction", label = TRUE, repel = TRUE)
ggsave(UMAP_plot, filename = paste0(outdir, "predicted_UMAP_pca_integrated.png"))

UMAP_plot_Location <- DimPlot(query_integrated, group.by = "MULTI_ID", label = TRUE, repel = TRUE)
ggsave(UMAP_plot_Location, filename = paste0(outdir, "predicted_UMAP_location_pca_integrated.png"))

UMAP_plot_Location_Time <- DimPlot(query_integrated, group.by = "Location_Time", label = TRUE, repel = TRUE)
ggsave(UMAP_plot_Location_Time, filename = paste0(outdir, "predicted_UMAP_location_time_pca_integrated.png"))


### Try predicting in SCT (4 days) from integrated ref (baseline) ###
query_integrated2sct <- scPredict(query, reference_integrated)

#### Generate specificity, sensitivity ####
df <- data.frame("Cell_type" = unique(query_integrated2sct$scpred_prediction), "n" = NA, "Sens" = NA, "Spec" = NA)
for (number in 0:4){
	df$n[which(df$Cell_type == number)] <- nrow(query_integrated2sct@meta.data[which(query_integrated2sct$scpred_prediction == number),])
	df$Sens[which(df$Cell_type == number)] <- nrow(query_integrated2sct@meta.data[which(query_integrated2sct$scpred_prediction == number & query_integrated2sct$Clusters_Res_0.09 == number),])/nrow(query_integrated2sct@meta.data[which(query_integrated2sct$Clusters_Res_0.09 == number),])
	df$Spec[which(df$Cell_type == number)] <- nrow(query_integrated2sct@meta.data[which(query_integrated2sct$scpred_prediction != number & query_integrated2sct$Clusters_Res_0.09 != number),])/nrow(query_integrated2sct@meta.data[which(query_integrated2sct$Clusters_Res_0.09 != number),])
}



###  Try predictin in integrated - newly integrated with just cell lines
seurat_list <- list()
for (line in unique(seurat_brisbane_4day@meta.data$Final_Assignment)){
	seurat_list[[line]] <- subset(seurat_brisbane_4day, subset = Final_Assignment == line)
}
message("The seurat objects for integration:")
print(seurat_list)

seurat_list <- lapply(seurat_list, function(x) {
	SCTransform(x, verbose = TRUE, vars.to.regress = c("scores.G1", "scores.S", "scores.G2M"))
})

seurat_features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)
seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = seurat_features)
seurat_list <- lapply(seurat_list, function(x){
	x[["integrated"]] <- NULL
	return(x)
})

seurat_anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT", anchor.features = seurat_features)


seurat_brisbane_4day_integrated <- IntegrateData(anchorset = seurat_anchors, normalization.method = "SCT")


### Normalize query data ###
query_integrated_separately <- scPredict(seurat_brisbane_4day, reference_integrated)

#### Generate specificity, sensitivity ####
df <- data.frame("Cell_type" = unique(query_integrated_separately$scpred_prediction), "n" = NA, "Sens" = NA, "Spec" = NA)
for (number in 0:4){
	df$n[which(df$Cell_type == number)] <- nrow(query_integrated_separately@meta.data[which(query_integrated_separately$scpred_prediction == number),])
	df$Sens[which(df$Cell_type == number)] <- nrow(query_integrated_separately@meta.data[which(query_integrated_separately$scpred_prediction == number & query_integrated_separately$Clusters_Res_0.09 == number),])/nrow(query_integrated_separately@meta.data[which(query_integrated_separately$Clusters_Res_0.09 == number),])
	df$Spec[which(df$Cell_type == number)] <- nrow(query_integrated_separately@meta.data[which(query_integrated_separately$scpred_prediction != number & query_integrated_separately$Clusters_Res_0.09 != number),])/nrow(query_integrated_separately@meta.data[which(query_integrated_separately$Clusters_Res_0.09 != number),])
}














##### Try predicting for each location & time separately #####
seurat$Locations <- gsub("\\d","", seurat$Site_rep)


#### Train the model ####
DefaultAssay(seurat_brisbane) <- "SCT"

seurat_brisbane <- seurat_brisbane %>%
	RunPCA() %>%
	RunUMAP(dims = 1:30)

reference <- getFeatureSpace(seurat_brisbane, "Clusters_Res_0.09")
reference <- trainModel(reference, model = "mda", number = 30)

get_probabilities(reference) %>% head()

get_scpred(reference)

plot_mda <- plot_probabilities(reference)

ggsave(plot_mda, filename = paste0(outdir, "scPred_Performance_mda.png"))


### Separate seurat object by Location_Time ###
seurat_list <- lapply(unique(seurat$Location_Time), function(x){
	subset(seurat, subset = Location_Time == x)
})
names(seurat_list) <- unique(seurat$Location_Time)

seurat_list <- lapply(seurat_list, function(x){
	SCTransform(x, verbose = TRUE, vars.to.regress = c("scores.G1", "scores.S", "scores.G2M"))
})

query_list <- lapply(seurat_list, function(x){
	scPredict(x, reference)
})

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
	ggsave(UMAP_plot_list[[x]], filename = paste0(outdir, x, "_predicted_UMAP_SCT.png"))
})

query_list <- lapply(query_list, function(x){
	x %>%
	RunPCA() %>%
	RunUMAP(., reduction = "pca", dims = 1:30)
}) 

UMAP_plot_Location_list <- lapply(query_list, function(x){
	DimPlot(x, group.by = "scpred_prediction", label = TRUE, repel = TRUE, split.by = "Final_Assignment")
})

lapply(names(UMAP_plot_Location_list), function(x){
	ggsave(UMAP_plot_Location_list[[x]], filename = paste0(outdir, x, "_predicted_UMAP_SCT_normal.png"), width = 21, height = 7)
})



##### Try predicting for each location & time separately #####
seurat_sydney$Locations <- gsub("\\d","", seurat$Site_rep)


#### Train the model ####
DefaultAssay(seurat_brisbane) <- "integrated"

# seurat_brisbane <- seurat_brisbane %>%
# 	SCTransform() %>%
# 	RunPCA() %>%
# 	RunUMAP(dims = 1:30)

reference <- getFeatureSpace(seurat_brisbane, "Clusters_Res_0.09")
reference <- trainModel(reference, model = "mda", number = 30)

get_probabilities(reference) %>% head()

get_scpred(reference)

plot_mda <- plot_probabilities(reference)

ggsave(plot_mda, filename = paste0(outdir, "scPred_Performance_mda.png"))


### Separate seurat object by Location_Time With Integrat + SCT ###
seurat_list <- lapply(unique(seurat$Location_Time), function(x){
	subset(seurat, subset = Location_Time == x)
})
names(seurat_list) <- unique(seurat$Location_Time)

seurat_list <- lapply(seurat_list, function(x){
	x %>%
	SCTransform()
})

seurat_list <- lapply(seurat_list, function(x){

})
    seurat_features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)
    seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = seurat_features)
    seurat_anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT", 
        anchor.features = seurat_features)
    seurat_integrated <- IntegrateData(anchorset = seurat_anchors, normalization.method = "SCT")


query_list <- lapply(seurat_list, function(x){
	scPredict(x, reference)
})

PCA_plot_list <- lapply(query_list, function(x){
	DimPlot(x, group.by = "scpred_prediction", reduction = "scpred")
})

lapply(names(PCA_plot_list), function(x){
	ggsave(PCA_plot_list[[x]], filename = paste0(outdir, x,"_predicted_PCA.png"))
})

query_list <- lapply(query_list, function(x){
	x %>%
	RunPCA() %>%
	RunUMAP(., reduction = "scpred", dims = 1:30)
}) 


UMAP_plot_list <- lapply(query_list, function(x){
	DimPlot(x, group.by = "scpred_prediction_integrated", label = TRUE, repel = TRUE)
})

ggsave(UMAP_plot3, filename = paste0(outdir, "predicted_UMAP_integrated.png"))

UMAP_plot_Location_Time3 <- DimPlot(query3, group.by = "Location_Time", label = TRUE, repel = TRUE)
ggsave(UMAP_plot_Location_Time3, filename = paste0(outdir, "predicted_UMAP_location_time_integrated.png"))











##### Using reference with both timepoints to predict in each  #####
file_list <- list.files(paste0(dir,"output/Seurat_w_Hash_freeze_thaw/CellCycle"), pattern = ".rds")
file_list <- file_list[!grepl("Brisbane_", file_list)]
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
reference <- getFeatureSpace(seurat_brisbane, "Clusters_Res_0.09")
reference <- trainModel(reference)

get_probabilities(reference) %>% head()

get_scpred(reference)

plot <- plot_probabilities(reference)

ggsave(plot, filename = paste0(outdir, "scPred_Performance_integrated_base_v4_caret.png"))

#### To retrain #####
reference <- trainModel(reference, model = "mda")


get_probabilities(reference) %>% head()

get_scpred(reference)

plot_mda <- plot_probabilities(reference)

ggsave(plot_mda, filename = paste0(outdir, "scPred_Performance_integrated_base_v4_mda.png"))


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
 
outdir_combined_int_ref <- paste0(outdir, "Brisbane_Combinede_Integrated/")
dir.create(outdir_combined_int_ref)


PCA_plot_list <- lapply(query_list, function(x){
	DimPlot(x, group.by = "scpred_prediction", reduction = "scpred")
})

lapply(names(PCA_plot_list), function(x){
	ggsave(PCA_plot_list[[x]], filename = paste0(outdir_combined_int_ref, x,"_predicted_Brisbane_int_ref_PCA.png"))
})

query_list <- lapply(query_list, function(x){
	RunUMAP(x, reduction = "scpred", dims = 1:30)
})

UMAP_plot_list <- lapply(query_list, function(x){
	DimPlot(x, group.by = "scpred_prediction", label = TRUE, repel = TRUE)
})

lapply(names(UMAP_plot_list), function(x){
	ggsave(UMAP_plot_list[[x]], filename = paste0(outdir_combined_int_ref, x, "_predicted_Brisbane_int_ref_UMAP.png"))
})

UMAP_plot_Assignment_list <- lapply(query_list, function(x){
	DimPlot(x, group.by = "scpred_prediction", label = TRUE, repel = TRUE, split.by = "Final_Assignment")
})

lapply(names(UMAP_plot_Assignment_list), function(x){
	ggsave(UMAP_plot_Assignment_list[[x]], filename = paste0(outdir_combined_int_ref, x, "_predicted_Brisbane_int_ref_UMAP_individual.png"), width = 21, height = 7)
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
	ggsave(UMAP_plot_Location_list[[x]], filename = paste0(outdir_combined_int_ref, x, "_predicted_Brisbane_int_ref_PCA_UMAP.png"), width = 21, height = 7)
})


##### Save the seurat objecst #####
saveRDS(reference, paste0(outdir_combined_int_ref, "Brisbane_Baseline_post_scpred.rds"))
# reference <- readRDS(paste0(outdir_combined_int_ref, "Brisbane_Baseline_post_scpred.rds"))

lapply(names(query_list), function(x){
	saveRDS(query_list[[x]], paste0(outdir_combined_int_ref, x, "_post_scpred.rds"))
})


# query_list <- list()
# query_list <- lapply(names, function(x){
# 	readRDS(paste0(outdir_combined_int_ref, x, "_post_scpred.rds"))
# })

##### Combine all data together to investigate similarities and differences
combined_assigned <- merge(reference, query_list)


##### Look at table of number of each cell type in each dataset #####
lapply(query_list, function(x){
	knitr::kable(prop.table(table(x$scpred_prediction)), "pipe")
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
