##### Predict by each cell line separately #####
library(Seurat)
library(tidyverse)
library(ggplot2)
library(scPred)
library(plyr)

dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
outdir <- paste0(dir,"output/scPred_Cluster_Assignments_Baseline_ref_Res0.13_5cluster_all_sep_SCT/")
dir.create(outdir)



##### Read in Seurat Objects #####
seurat_brisbane_baseline <- readRDS(paste0(dir,"output/Brisbane_Clustering_Testing/Make_Brisbane_Baseline_Seurat_covs_res_0.13/Brisbane_Seurat_SCT_CellCycle_cov_Indiv_integration_w_meta.rds"))


file_list <- list.files(paste0(dir,"output/Seurat_w_Hash_freeze_thaw/CellCycle"), pattern = ".rds")
file_list <- file_list[!grepl("Brisbane_Baseline", file_list)]
file_list <- file_list[!grepl("seurat_filtered_cell_cycle.rds", file_list)]
seurat_list <- lapply(file_list, function(x){
	readRDS(paste0(dir,"output/Seurat_w_Hash_freeze_thaw/CellCycle/",x))
})
names(seurat_list) <- gsub("_seurat.rds", "", file_list) %>% gsub(".rds","",.)
names <- gsub("_seurat.rds", "", file_list) %>% gsub(".rds","",.)



##### Subset brisbane baseline by cell line and compute SCT #####
seurat_brisbane_baseline_list <- lapply(unique(seurat_brisbane_baseline@meta.data$Final_Assignment), function(line){
	subset(seurat_brisbane_baseline, subset = Final_Assignment == line)
})
names(seurat_brisbane_baseline_list) <- unique(seurat_brisbane_baseline@meta.data$Final_Assignment)


seurat_brisbane_baseline_list <- lapply(seurat_brisbane_baseline_list, function(x){
	print(x)
	SCTransform(x, verbose = TRUE, vars.to.regress = c("scores.G1", "scores.S", "scores.G2M", "percent.mt"))
})

seurat_brisbane_baseline_list <- lapply(seurat_brisbane_baseline_list, function(x){
	RunPCA(x, features = VariableFeatures(object = x)) 
})





## Add location to the seurat metadata ##
seurat_list_line <- list()
i <- 1
for (x in names(seurat_list)){
	print(x)
	seurat_list[[x]]$Location <- gsub("\\d","", seurat_list[[x]]$Site_rep)
	for (line in unique(seurat_list[[x]]$Final_Assignment)){
		if (i == 1){
			seurat_list_line[[line]] <- list()
		}
		seurat_list_line[[line]][[x]] <- subset(seurat_list[[x]], subset = Final_Assignment == line)
	}
	i <- i + 1
}



##### Train the model #####
reference <- lapply(seurat_brisbane_baseline_list, function(x){
	getFeatureSpace(x, "Clusters_Res_0.13")
})

reference <- lapply(reference, function(x){
	trainModel(x)
})

lapply(reference, function(x){
	get_scpred(x)
})

plot_list <- lapply(reference, function(x){
	plot_probabilities(x)
})

lapply(names(plot_list), function(x){
	ggsave(plot_list[[x]], filename = paste0(outdir, x, "_scPred_Performance_caret.png"))
})


#### To retrain #####
reference <- lapply(seurat_brisbane_baseline_list, function(x){
	trainModel(x, model = "mda")
})

lapply(reference, function(x){
	get_scpred(x)
})

plot_list <- lapply(reference, function(x){
	plot_probabilities(x)
})

lapply(names(plot_list), function(x){
	ggsave(plot_list[[x]], filename = paste0(outdir, x, "_scPred_Performance_mda.png"))
})




### Check prediction for Brisbane Village Day 4 ###
query_list <- lapply(seurat_list_line, function(x){
	lapply(x, function(y){
		SCTransform(y, verbose = TRUE, vars.to.regress = c("scores.G1", "scores.S", "scores.G2M", "percent.mt"))
	})
})

# seurat_features <- SelectIntegrationFeatures(object.list = temp_list, nfeatures = 6000)
# temp_list <- PrepSCTIntegration(object.list = temp_list, anchor.features = seurat_features)
# seurat_anchors <- FindIntegrationAnchors(object.list = temp_list, normalization.method = "SCT", anchor.features = seurat_features)

# integrated <- IntegrateData(anchorset = seurat_anchors, normalization.method = "SCT")
# return(integrated)

query_list <- lapply(names(query_list), function(x){
	lapply(query_list[[x]], function(y){
		scPredict(y, reference[[x]])
	})
})
names(query_list) <- names(reference)
 

PCA_plot_list <- lapply(query_list, function(x){
	lapply(x, function(y){
		DimPlot(y, group.by = "scpred_prediction", reduction = "scpred")
	})
})


lapply(names(PCA_plot_list), function(x){
	lapply(names(PCA_plot_list[[x]]), function(y){
		ggsave(PCA_plot_list[[x]][[y]], filename = paste0(outdir, x,"_", y, "_predicted_PCA.png"))
	})
})


query_list <- lapply(query_list, function(x){
	lapply(x, function(y){
		RunUMAP(y, reduction = "scpred", dims = 1:30)
	})
})


UMAP_plot_list <- lapply(query_list, function(x){
	lapply(x, function(y){
		DimPlot(y, group.by = "scpred_prediction", label = TRUE, repel = TRUE)
	})
})

lapply(names(UMAP_plot_list), function(x){
	lapply(names(UMAP_plot_list[[x]]), function(y){
		ggsave(UMAP_plot_list[[x]][[y]], filename = paste0(outdir, x,"_", y, "_predicted_UMAP.png"))
	})
})


query_list <- lapply(query_list, function(x){
	lapply(x, function(y){
		y %>%
		RunPCA() %>%
		RunUMAP(., reduction = "pca", dims = 1:30, reduction.name = "pca_umap")
	})
}) 

UMAP_plot_Location_list <- lapply(query_list, function(x){
	lapply(x, function(y){
		DimPlot(y, group.by = "scpred_prediction", label = TRUE, repel = TRUE,  reduction = "pca_umap")
	})
})

lapply(names(UMAP_plot_Location_list), function(x){
	lapply(names(UMAP_plot_Location_list[[x]]), function(y){
		ggsave(UMAP_plot_Location_list[[x]][[y]], filename = paste0(outdir, x,"_", y, "_predicted_PCA_UMAP.png"))
	})
})


##### Save the seurat objecst #####
saveRDS(reference, paste0(outdir, "Brisbane_Baseline_post_scpred.rds"))
# reference <- readRDS(paste0(outdir, "Brisbane_Baseline_post_scpred.rds"))

lapply(names(query_list), function(x){
	saveRDS(query_list[[x]], paste0(outdir, x, "_post_scpred.rds"))
})








##### Check to see if any other covariates are confounding the PCA #####
pca_cov_outdir <- paste0(outdir,"scpred_reduc_covariates/")
dir.create(pca_cov_outdir)

scpred_plot_cov_list <- lapply(query_list, function(x){
	lapply(x, function(y){
		temp <- list()
			for (cov in c("phases", "Final_Assignment","Pool","Time" )){
				temp[[cov]] <- DimPlot(y, group.by = cov, reduction = "scpred")
			}
			for (cov in c( "percent.rb", "percent.mt", "nCount_RNA", "nFeature_RNA"))
			{
				temp[[cov]] <- FeaturePlot(y, reduction = "scpred", features = cov)
			}
			return(temp)
	})
})

lapply(names(scpred_plot_cov_list), function(x){
	lapply(names(scpred_plot_cov_list[[x]]), function(y){
		lapply(names(scpred_plot_cov_list[[x]][[y]]), function(cov){
			ggsave(scpred_plot_cov_list[[x]][[y]][[cov]], filename = paste0(pca_cov_outdir, x,"_",y,"_",cov,"_predicted_PCA.png"))
		})
	})
})


pca_cov_outdir <- paste0(outdir,"scpred_reduc_umap_covariates/")
dir.create(pca_cov_outdir)

scpred_umap_cov_list <- lapply(query_list, function(x){
	lapply(x, function(y){
		temp <- list()
		for (cov in c("phases", "Final_Assignment","Pool","Time" )){
			temp[[cov]] <- DimPlot(y, group.by = cov, reduction = "umap")
		}
		for (cov in c( "percent.rb", "percent.mt", "nCount_RNA", "nFeature_RNA"))
		{
			temp[[cov]] <- FeaturePlot(y, reduction = "umap", features = cov)
		}
		return(temp)
	})
})

lapply(names(scpred_umap_cov_list), function(x){
	lapply(names(scpred_plot_cov_list[[x]]), function(y){
		lapply(names(scpred_plot_cov_list[[x]][[y]]), function(cov){
			ggsave(scpred_umap_cov_list[[x]][[cov]], filename = paste0(pca_cov_outdir, x,"_",cov,"_predicted_UMAP.png"))
		})
	})
})




pca_cov_outdir <- paste0(outdir,"PCA_covariates/")
dir.create(pca_cov_outdir)

PCA_plot_cov_list <- lapply(query_list, function(x){
	lapply(x, function(y){
		temp <- list()
		for (cov in c("phases", "Final_Assignment","Pool","Time", "scpred_prediction" )){
			temp[[cov]] <- DimPlot(y, group.by = cov, reduction = "pca")
		}
		for (cov in c( "percent.rb", "percent.mt", "nCount_RNA", "nFeature_RNA"))
		{
			temp[[cov]] <- FeaturePlot(y, reduction = "pca", features = cov)
		}
		return(temp)
	})
})

lapply(names(PCA_plot_cov_list), function(x){
	lapply(names(scpred_plot_cov_list[[x]]), function(y){
		lapply(names(scpred_plot_cov_list[[x]][[y]]), function(cov){
			ggsave(PCA_plot_cov_list[[x]][[cov]], filename = paste0(pca_cov_outdir, x,"_",cov,"_predicted_PCA.png"))
		})
	})
})



##### Look at table of number of each cell type in each dataset #####
prop_table_list <- lapply(names(query_list), function(x){
	lapply(names(query_list[[x]]), function(y){
		temp <- as.data.frame(prop.table(table(query_list[[x]][[y]]$scpred_prediction)))
		colnames(temp) <- c("Cluster",paste0(y,"_",x))
		return(temp)
	})
})


prop_table_list2 <- lapply(prop_table_list, function(x){
	join_all(x)
})

prop_table <- join_all(prop_table_list2)

write_delim(prop_table, paste0(outdir,"Cluster_proportions.tsv"))

prop_table_long <- pivot_longer(prop_table, names_to = "Location_Time", values_to = "Proportion", cols = colnames(prop_table)[2:ncol(prop_table)])

prop_table_Brisbane_Baseline <- as.data.frame(prop.table(table(seurat_brisbane_baseline$Clusters_Res_0.13, seurat_brisbane_baseline$Final_Assignment), margin = 2))
colnames(prop_table_Brisbane_Baseline) <- c("Cluster","Location_Time", "Proportion")
prop_table_Brisbane_Baseline$Location_Time <- paste0("Brisbane_Baseline_", prop_table_Brisbane_Baseline$Location_Time)


prop_table_long <- rbind(prop_table_long, prop_table_Brisbane_Baseline)


prop_table_long$Location_Time <- factor(prop_table_long$Location_Time, levels = c("Brisbane_Baseline_FSA0006", "Brisbane_Baseline_MBE1006", "Brisbane_Baseline_TOB0421",
																					"Brisbane_Village_Day_4_FSA0006", "Brisbane_Village_Day_4_MBE1006", "Brisbane_Village_Day_4_TOB0421", 
																					"Melbourne_Baseline_FSA0006", "Melbourne_Baseline_MBE1006", "Melbourne_Baseline_TOB0421", 
																					"Melbourne_Village_Day_4_TOB0421", "Melbourne_Village_Day_4_MBE1006", "Melbourne_Village_Day_4_FSA0006",
																					"Sydney_Baseline_FSA0006", "Sydney_Baseline_MBE1006", "Sydney_Baseline_TOB0421",
																					"Sydney_Village_Day_4_FSA0006", "Sydney_Village_Day_4_MBE1006", "Sydney_Village_Day_4_TOB0421",
																					"Sydney_Thawed_Village_Day_0_FSA0006", "Sydney_Thawed_Village_Day_0_MBE1006", "Sydney_Thawed_Village_Day_0_TOB0421",
																					"Sydney_Thawed_Village_Day_7_FSA0006", "Sydney_Thawed_Village_Day_7_MBE1006", "Sydney_Thawed_Village_Day_7_TOB0421"
																					))


prop_plot <- ggplot(prop_table_long, aes(Location_Time, Proportion, fill = Cluster)) +
	geom_bar(stat = "identity",position = "stack") +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 45, hjust=1))

ggsave(prop_plot, filename = paste0(outdir, "Cluster_proportions.png"))
