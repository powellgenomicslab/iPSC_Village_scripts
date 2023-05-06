##### Read in Arguments #####
print("Reading and assigning input arguments")


dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
outdir <- paste0(dir,"output/All_data_integrated_remove_bad/")
datadir <- paste0(dir,"output/All_data_integrated/")
datadir_multi_cluster <- paste0(datadir,"/Multi_Resolution/")

dir.create(outdir)


library(Seurat)
library(tidyverse)
library(ggplot2)
library(jcolors)
library(cowplot)
library(RColorBrewer)
library(readr)
library(purrr)
library(clustree)
library(reticulate)
library(awtools)
library(Nebulosa)
library(plyr)




print(sessionInfo())


### Set names ###
pools <- dir("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/")
old_name <- c("22_FSA", "29_MBE", "36_TOB00421_i_E8")
new_name <- c("FSA0006","MBE1006","TOB0421")
cell_key <- data.frame(old_name, new_name)


##### Read in Data #####
seurat_integrated <- readRDS(paste0(datadir, "seurat_integrated_all_times.rds"))
clusters <- read_delim(paste0(datadir_multi_cluster, "/Identities_Resolution_0.08.txt"), delim = "\t")


##### Decided to use remove two of the clusters - one with low mt % and one with high Mt % #####
rownames(clusters) <- rownames(seurat_integrated@meta.data)
seurat_integrated <- AddMetaData(seurat_integrated, clusters)

seurat_integrated <- subset(seurat_integrated, subset = Resolution_0.08 < 4) ### removes 584 cells


##### Rerun scaling and 
seurat_integrated <- ScaleData(seurat_integrated, verbose = FALSE, vars.to.regress = c("scores.G1", "scores.S", "scores.G2M", "percent.mt"))
seurat_integrated <- RunPCA(seurat_integrated, verbose = FALSE)
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:30)
seurat_integrated <- FindNeighbors(seurat_integrated, reduction = "pca", assay = "integrated", dims = 1:20)







##### PCAs #####
plots <- DimPlot(seurat_integrated, reduction = "pca", group.by = c("Final_Assignment"))
ggsave(plots, filename = paste0(outdir,"pca.png"),width = 9, height = 9)

pca <- DimPlot(seurat_integrated, reduction = "pca", group.by = "Pool", pt.size = .01, split.by = 'Pool')
ggsave(pca, filename = paste0(outdir, "pca_pool.png"), height = 4, width = 10)

pca <- DimPlot(seurat_integrated, reduction = "pca", group.by = "Time", pt.size = .01, split.by = 'Time')
ggsave(pca, filename = paste0(outdir, "pca_Time.png"), height = 4, width = 10)

pca <- DimPlot(seurat_integrated, reduction = "pca", group.by = "MULTI_ID", pt.size = .01, split.by = 'MULTI_ID')
ggsave(pca, filename = paste0(outdir, "pca_SiteSampleReplicate.png"), height = 4, width = 10)

pca <- DimPlot(seurat_integrated, reduction = "pca", group.by = "Final_Assignment", pt.size = .01, split.by = 'Final_Assignment')
ggsave(pca, filename = paste0(outdir, "pca_CellLine.png"), height = 4, width = 10)

pca <- DimPlot(seurat_integrated, reduction = "pca", group.by = "Time", pt.size = .001, split.by = 'MULTI_ID')
ggsave(pca, filename = paste0(outdir, "pca_SiteSampleReplicate_Time.png"), height = 4, width = 10)

pca <- DimPlot(seurat_integrated, reduction = "pca", group.by = "Time", pt.size = .001, split.by = 'Final_Assignment')
ggsave(pca, filename = paste0(outdir, "pca_CellLine_Time.png"), height = 4, width = 10)

seurat_integrated$Time_Final_Assignment <- paste0(seurat_integrated$Time, seurat_integrated$Final_Assignment)
pca <- DimPlot(seurat_integrated, reduction = "pca", group.by = "Time_Final_Assignment", pt.size = .01, split.by = 'MULTI_ID')
ggsave(pca, filename = paste0(outdir, "pca_CellLine_Time_Site.png"), height = 4, width = 21)

pca <- DimPlot(seurat_integrated,  pt.size = .01, group.by = c("phases"))  
ggsave(pca, filename = paste0(outdir, "pca_CellCycle.png"))

mt_pca <- FeaturePlot(seurat_integrated, features = "percent.mt")
ggsave(mt_pca, filename = paste0(outdir,"mt_percent_pca_seurat.png"))

rb_pca <- FeaturePlot(seurat_integrated, features = "percent.rb")
ggsave(rb_pca, filename = paste0(outdir,"rb_percent_pca_seurat.png"))

libsize_pca <- FeaturePlot(seurat_integrated, features = "nCount_RNA")
ggsave(libsize_pca, filename = paste0(outdir,"libsize_pca_seurat.png"))

Ngene_pca <- FeaturePlot(seurat_integrated, features = "nFeature_RNA")
ggsave(Ngene_pca, filename = paste0(outdir,"Ngene_pca_seurat.png"))

Brisbane_clusters_pca <- DimPlot(seurat_integrated, reduction = "pca", group.by = "Clusters_Res_0.13", cells = rownames(seurat_integrated@meta.data)[which(!is.na(seurat_integrated@meta.data$Clusters_Res_0.13))])
ggsave(Brisbane_clusters_pca, filename = paste0(outdir,"Brisbane_clusters_pca_seurat.png"))



##### UMAPs #####
plots <- DimPlot(seurat_integrated, group.by = c("Day", "Final_Assignment"), combine = FALSE)
plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 4, 
	byrow = TRUE, override.aes = list(size = 2.5))))
ggsave(CombinePlots(plots), filename = paste0(outdir,"umap_integrated.png"),width = 12, height = 9)


plots <- DimPlot(seurat_integrated, group.by = c("Final_Assignment"))
ggsave(plots, filename = paste0(outdir,"umap.png"),width = 20, height = 9)

umap <- DimPlot(seurat_integrated, reduction = "umap", group.by = "Pool", pt.size = .01, split.by = 'Pool')
ggsave(umap, filename = paste0(outdir, "umap_pool.png"), height = 4, width = 10)

umap <- DimPlot(seurat_integrated, reduction = "umap", group.by = "Time", pt.size = .01, split.by = 'Time')
ggsave(umap, filename = paste0(outdir, "umap_Time.png"), height = 4, width = 10)

umap <- DimPlot(seurat_integrated, reduction = "umap", group.by = "MULTI_ID", pt.size = .01, split.by = 'MULTI_ID')
ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate.png"), height = 4, width = 10)

umap <- DimPlot(seurat_integrated, reduction = "umap", group.by = "Final_Assignment", pt.size = .01, split.by = 'Final_Assignment')
ggsave(umap, filename = paste0(outdir, "umap_CellLine.png"), height = 4, width = 10)

umap <- DimPlot(seurat_integrated, reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'MULTI_ID')
ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate_Time.png"), height = 4, width = 10)

umap <- DimPlot(seurat_integrated, reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'Final_Assignment')
ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time.png"), height = 4, width = 10)

seurat_integrated$Time_Final_Assignment <- paste0(seurat_integrated$Time, seurat_integrated$Final_Assignment)
umap <- DimPlot(seurat_integrated, reduction = "umap", group.by = "Time_Final_Assignment", pt.size = .01, split.by = 'MULTI_ID')
ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time_Site.png"), height = 4, width = 21)

umap <- DimPlot(seurat_integrated,  pt.size = .01, group.by = c("phases"))  
ggsave(umap, filename = paste0(outdir, "umap_CellCycle.png"))

mt_umap <- FeaturePlot(seurat_integrated, features = "percent.mt")
ggsave(mt_umap, filename = paste0(outdir,"mt_percent_umap_seurat.png"))

rb_umap <- FeaturePlot(seurat_integrated, features = "percent.rb")
ggsave(rb_umap, filename = paste0(outdir,"rb_percent_umap_seurat.png"))

libsize_umap <- FeaturePlot(seurat_integrated, features = "nCount_RNA")
ggsave(libsize_umap, filename = paste0(outdir,"libsize_umap_seurat.png"))

Ngene_umap <- FeaturePlot(seurat_integrated, features = "nFeature_RNA")
ggsave(Ngene_umap, filename = paste0(outdir,"Ngene_umap_seurat.png"))

Brisbane_clusters_umap <- DimPlot(seurat_integrated, group.by = "Clusters_Res_0.13", cells = rownames(seurat_integrated@meta.data)[which(!is.na(seurat_integrated@meta.data$Clusters_Res_0.13))])
ggsave(Brisbane_clusters_umap, filename = paste0(outdir,"Brisbane_clusters_umap_seurat.png"))








##### Test multiple resolutions for clustering
outdir_multi_cluster <- paste0(outdir,"/Multi_Resolution/")
dir.create(outdir_multi_cluster)

### Find the clusters for this resolution ###
for (res in seq(0,0.5,0.01)){
	print(paste0("Finding Clusters for resolution=",res))
	seurat_integrated <- FindClusters(seurat_integrated, resolution = res)

	### Plot the UMAP for this resolution ###
	UMAP <- DimPlot(seurat_integrated, reduction = "umap")
	ggsave(file = paste0(outdir_multi_cluster,"UMAPresolution",res,".png"), plot = UMAP)

	### Save the identities for this resolution ###
	identities <- as.data.frame(Idents(seurat_integrated))
	colnames(identities) <- c(paste0("Resolution_",res))

	write_delim(identities, paste0(outdir_multi_cluster,"Identities_Resolution_",res,".txt"), delim = "\t")

}



##### Clustree for multiple resolutions for clustering #####
outdir_clustree <- paste0(outdir,"/Clustree/")
dir.create(outdir_clustree)

ResolutionsFileList <- list.files(outdir_multi_cluster, pattern = "Identities_Resolution_")

ResolutionsList <- lapply(ResolutionsFileList, FUN = function(x){
	read_delim(paste0(outdir_multi_cluster, x), delim = "\t")
})

ResolutionsDF <- do.call(cbind, ResolutionsList)

##### CLUSTREE #####
print("Using Clustree to diagram changes with increasing resolutions")
pClustree <- clustree(ResolutionsDF, prefix = "Resolution_")
ggsave(file = paste0(outdir_clustree, "ClusteringClustree.png"), plot = pClustree)






##### Generate figures for DE genes at each resolution #####
outdir_DEGs <- paste0(outdir,"/DEGs/")
dir.create(outdir_DEGs)

resolutions <- c(0.03, 0.05, 0.07, 0.08, 0.11, 0.16, 0.17)


resolutions_assignments_list <- lapply(resolutions, function(res){
	read_delim(paste0(outdir_multi_cluster, "/Identities_Resolution_",res,".txt"), delim = "\t")
})

ResolutionsDF <- do.call(cbind, resolutions_assignments_list)

if(!file.exists(paste0(outdir_DEGs,"DEG_dif_res.rds"))){
	degs_list <- lapply(resolutions, function(x){
		Idents(seurat_integrated) <- ResolutionsDF[,paste0("Resolution_", x)]
		FindAllMarkers(seurat_integrated, assay = "RNA")
	})
	names(degs_list) <- paste0("Resolution_", resolutions)

	saveRDS(degs_list, paste0(outdir_DEGs,"DEG_dif_res.rds"))
} else {
	degs_list <- readRDS(paste0(outdir_DEGs,"DEG_dif_res.rds"))
}

degs_list_res <- lapply(names(degs_list), function(x){
	colnames(degs_list[[x]]) <- c(paste0(colnames(degs_list[[x]])[1:(ncol(degs_list[[x]])-1)],"_",x),colnames(degs_list[[x]])[ncol(degs_list[[x]])])
	return(degs_list[[x]])
})
names(degs_list_res) <- paste0("Resolution_", resolutions)


##### Add gene IDs for easy identification downstream #####
GeneConversion1 <- read_delim("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/DRENEA_1/outs/filtered_feature_bc_matrix/features.tsv.gz", col_names = F, delim = "\t")
GeneConversion2 <- read_delim("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/Village_B_1_week/outs/filtered_feature_bc_matrix/features.tsv.gz", col_names = F, delim = "\t")

GeneConversion <- unique(rbind(GeneConversion1, GeneConversion2))
GeneConversion <- GeneConversion[!duplicated(GeneConversion$X1),]
GeneConversion <- left_join(data.frame(rownames(seurat_integrated[["RNA"]])), GeneConversion, by= c("rownames.seurat_integrated...RNA...." = "X1"))
GeneConversion$X3 <- NULL
colnames(GeneConversion) <- c("ENSG_ID", "Gene_ID")


degs_list <- lapply(degs_list, function(x){
	left_join(x, GeneConversion, by = c("gene" = "ENSG_ID"))
})


TopMarkers_list <- lapply(names(degs_list), function(x){
	degs_list[[x]] %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
})
names(TopMarkers_list) <- names(degs_list)



print("Making Plots")


DefaultAssay(seurat_integrated) <- "RNA"

pTopDEgenes <- list()
for (res in names(TopMarkers_list)){
	pTopDEgenes[[res]] <- list()
	for (clust in unique(degs_list[[res]]$cluster)){
		print(clust)
			for (gene in TopMarkers_list[[res]]$gene[which(TopMarkers_list[[res]]$cluster == clust)]){
			# print(TopMarkers$gene[which(TopMarkers$cluster == clust)])
				pTopDEgenes[[res]][[paste0("Cluster",clust,"_",unique(TopMarkers_list[[res]][which(TopMarkers_list[[res]]$gene == gene), "Gene_ID"]))]] <- FeaturePlot(seurat_integrated, features = gene) + labs(title = unique(TopMarkers_list[[res]][which(TopMarkers_list[[res]]$gene == gene), "Gene_ID"]))
			}
	}
}


pTopDEgenesNebulosa <- list()
for (res in names(TopMarkers_list)){
	pTopDEgenesNebulosa[[res]] <- list()
	for (clust in unique(degs_list[[res]]$cluster)){
		print(clust)
			for (gene in TopMarkers_list[[res]]$gene[which(TopMarkers_list[[res]]$cluster == clust)]){
			# print(TopMarkers$gene[which(TopMarkers$cluster == clust)])
				pTopDEgenesNebulosa[[res]][[paste0("Cluster",clust,"_",unique(TopMarkers_list[[res]][which(TopMarkers_list[[res]]$gene == gene), "Gene_ID"]))]] <- plot_density(seurat_integrated, gene) + labs(title = unique(TopMarkers_list[[res]][which(TopMarkers_list[[res]]$gene == gene), "Gene_ID"]))
			}
	}
}



print("Saving Plots")
for (res in names(pTopDEgenes)){
	dir.create(paste0(outdir_DEGs, res,"/"))
	lapply(names(pTopDEgenes[[res]]), FUN = function(x){
		ggsave(pTopDEgenes[[res]][[x]], file = paste0(outdir_DEGs,res,"/TopDEgenes",x,".png"))
	})
}



print("Saving Plots")
for (res in names(pTopDEgenesNebulosa)){
	dir.create(paste0(outdir_DEGs, res,"/"))
	lapply(names(pTopDEgenesNebulosa[[res]]), FUN = function(x){
		ggsave(pTopDEgenesNebulosa[[res]][[x]], file = paste0(outdir_DEGs,res,"/Nebulosa_TopDEgenes",x,".png"))
	})
}



##### Generate figures for DE genes at each resolution #####
outdir_DEGs <- paste0(outdir,"/DEGs/")
dir.create(outdir_DEGs)

resolutions <- c(0.28)


resolutions_assignments_list <- lapply(resolutions, function(res){
	read_delim(paste0(outdir_multi_cluster, "/Identities_Resolution_",res,".txt"), delim = "\t")
})

ResolutionsDF <- do.call(cbind, resolutions_assignments_list)

if(!file.exists(paste0(outdir_DEGs,"DEG_dif_res_0.28.rds"))){
	degs_list <- lapply(resolutions, function(x){
		Idents(seurat_integrated) <- ResolutionsDF[,paste0("Resolution_", x)]
		FindAllMarkers(seurat_integrated, assay = "RNA")
	})
	names(degs_list) <- paste0("Resolution_", resolutions)

	saveRDS(degs_list, paste0(outdir_DEGs,"DEG_dif_res_0.28.rds"))
} else {
	degs_list <- readRDS(paste0(outdir_DEGs,"DEG_dif_res_0.28.rds"))
}

degs_list_res <- lapply(names(degs_list), function(x){
	colnames(degs_list[[x]]) <- c(paste0(colnames(degs_list[[x]])[1:(ncol(degs_list[[x]])-1)],"_",x),colnames(degs_list[[x]])[ncol(degs_list[[x]])])
	return(degs_list[[x]])
})
names(degs_list_res) <- paste0("Resolution_", resolutions)


##### Add gene IDs for easy identification downstream #####
GeneConversion1 <- read_delim("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/DRENEA_1/outs/filtered_feature_bc_matrix/features.tsv.gz", col_names = F, delim = "\t")
GeneConversion2 <- read_delim("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/Village_B_1_week/outs/filtered_feature_bc_matrix/features.tsv.gz", col_names = F, delim = "\t")

GeneConversion <- unique(rbind(GeneConversion1, GeneConversion2))
GeneConversion <- GeneConversion[!duplicated(GeneConversion$X1),]
GeneConversion <- left_join(data.frame(rownames(seurat_integrated[["RNA"]])), GeneConversion, by= c("rownames.seurat_integrated...RNA...." = "X1"))
GeneConversion$X3 <- NULL
colnames(GeneConversion) <- c("ENSG_ID", "Gene_ID")


degs_list <- lapply(degs_list, function(x){
	left_join(x, GeneConversion, by = c("gene" = "ENSG_ID"))
})

TopMarkers_list <- lapply(names(degs_list), function(x){
	degs_list[[x]] %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
})
names(TopMarkers_list) <- names(degs_list)



print("Making Plots")


DefaultAssay(seurat_integrated) <- "RNA"

pTopDEgenes <- list()
for (res in names(TopMarkers_list)){
	pTopDEgenes[[res]] <- list()
	for (clust in unique(degs_list[[res]]$cluster)){
		print(clust)
			for (gene in TopMarkers_list[[res]]$gene[which(TopMarkers_list[[res]]$cluster == clust)]){
			# print(TopMarkers$gene[which(TopMarkers$cluster == clust)])
				pTopDEgenes[[res]][[paste0("Cluster",clust,"_",unique(TopMarkers_list[[res]][which(TopMarkers_list[[res]]$gene == gene), "Gene_ID"]))]] <- FeaturePlot(seurat_integrated, features = gene) + labs(title = unique(TopMarkers_list[[res]][which(TopMarkers_list[[res]]$gene == gene), "Gene_ID"]))
			}
	}
}


pTopDEgenesNebulosa <- list()
for (res in names(TopMarkers_list)){
	pTopDEgenesNebulosa[[res]] <- list()
	for (clust in unique(degs_list[[res]]$cluster)){
		print(clust)
			for (gene in TopMarkers_list[[res]]$gene[which(TopMarkers_list[[res]]$cluster == clust)]){
			# print(TopMarkers$gene[which(TopMarkers$cluster == clust)])
				pTopDEgenesNebulosa[[res]][[paste0("Cluster",clust,"_",unique(TopMarkers_list[[res]][which(TopMarkers_list[[res]]$gene == gene), "Gene_ID"]))]] <- plot_density(seurat_integrated, gene) + labs(title = unique(TopMarkers_list[[res]][which(TopMarkers_list[[res]]$gene == gene), "Gene_ID"]))
			}
	}
}



print("Saving Plots")
for (res in names(pTopDEgenes)){
	dir.create(paste0(outdir_DEGs, res,"/"))
	lapply(names(pTopDEgenes[[res]]), FUN = function(x){
		ggsave(pTopDEgenes[[res]][[x]], file = paste0(outdir_DEGs,res,"/TopDEgenes",x,".png"))
	})
}



print("Saving Plots")
for (res in names(pTopDEgenesNebulosa)){
	dir.create(paste0(outdir_DEGs, res,"/"))
	lapply(names(pTopDEgenesNebulosa[[res]]), FUN = function(x){
		ggsave(pTopDEgenesNebulosa[[res]][[x]], file = paste0(outdir_DEGs,res,"/Nebulosa_TopDEgenes",x,".png"))
	})
}








pluri_genes <- read_delim(paste0(dir,"data/pluripotency_genes.tsv"), delim = "\t", col_names = "Gene_ID")
meso_genes <- read_delim(paste0(dir,"data/mesoderm_markers.tsv"), delim = "\t", col_names = "Gene_ID")
endo_genes <- read_delim(paste0(dir,"data/endoderm_markers.tsv"), delim = "\t", col_names = "Gene_ID")
ecto_genes <- read_delim(paste0(dir,"data/ectoderm_markers.tsv"), delim = "\t", col_names = "Gene_ID")


##### Save Pluri Genes #####
pPluriGenes <- list()
for (gene in pluri_genes$Gene_ID){
	pPluriGenes[[gene]] <- FeaturePlot(seurat_integrated, GeneConversion[which(GeneConversion$Gene_ID == gene), "ENSG_ID"]) + labs(title = gene)
}

pPluriGenesNebulosa <- list()
for (gene in pluri_genes$Gene_ID){
	print(gene)
	pPluriGenesNebulosa[[gene]] <- plot_density(seurat_integrated, GeneConversion[which(GeneConversion$Gene_ID == gene), "ENSG_ID"]) + labs(title = gene)
}

dir.create(paste0(outdir_DEGs, "pluripotency_genes/"))

print("Saving Plots")
lapply(names(pPluriGenes), FUN = function(x){
	ggsave(pPluriGenes[[x]], file = paste0(outdir_DEGs, "pluripotency_genes/",x,".png"))
})

lapply(names(pPluriGenesNebulosa), FUN = function(x){
	ggsave(pPluriGenesNebulosa[[x]], file = paste0(outdir_DEGs, "pluripotency_genes/Nebulosa_",x,".png"))
})


##### Save meso Genes #####
pMesoGenes <- list()
for (gene in meso_genes$Gene_ID){
	pMesoGenes[[gene]] <- FeaturePlot(seurat_integrated, GeneConversion[which(GeneConversion$Gene_ID == gene), "ENSG_ID"]) + labs(title = gene)
}

pMesoGenesNebulosa <- list()
for (gene in meso_genes$Gene_ID){
	print(gene)
	pMesoGenesNebulosa[[gene]] <- plot_density(seurat_integrated, GeneConversion[which(GeneConversion$Gene_ID == gene), "ENSG_ID"]) + labs(title = gene)
}

dir.create(paste0(outdir_DEGs, "mesoderm_genes/"))

print("Saving Plots")
lapply(names(pMesoGenes), FUN = function(x){
	ggsave(pMesoGenes[[x]], file = paste0(outdir_DEGs, "mesoderm_genes/",x,".png"))
})

lapply(names(pMesoGenesNebulosa), FUN = function(x){
	ggsave(pMesoGenesNebulosa[[x]], file = paste0(outdir_DEGs, "mesoderm_genes/Nebulosa_",x,".png"))
})


##### Save Endoderm Genes #####
pEndoGenes <- list()
for (gene in endo_genes$Gene_ID){
	pEndoGenes[[gene]] <- FeaturePlot(seurat_integrated, GeneConversion[which(GeneConversion$Gene_ID == gene), "ENSG_ID"]) + labs(title = gene)
}

pEndoGenesNebulosa <- list()
for (gene in endo_genes$Gene_ID){
	print(gene)
	pEndoGenesNebulosa[[gene]] <- plot_density(seurat_integrated, GeneConversion[which(GeneConversion$Gene_ID == gene), "ENSG_ID"]) + labs(title = gene)
}

dir.create(paste0(outdir_DEGs, "endoderm_genes/"))

print("Saving Plots")
lapply(names(pEndoGenes), FUN = function(x){
	ggsave(pEndoGenes[[x]], file = paste0(outdir_DEGs, "endoderm_genes/",x,".png"))
})

lapply(names(pEndoGenesNebulosa), FUN = function(x){
	ggsave(pEndoGenesNebulosa[[x]], file = paste0(outdir_DEGs, "endoderm_genes/Nebulosa_",x,".png"))
})



##### Save Ectoderm Genes #####
pEctoGenes <- list()
for (gene in ecto_genes$Gene_ID){
	pEctoGenes[[gene]] <- FeaturePlot(seurat_integrated, GeneConversion[which(GeneConversion$Gene_ID == gene), "ENSG_ID"]) + labs(title = gene)
}

pEctoGenesNebulosa <- list()
for (gene in ecto_genes$Gene_ID){
	print(gene)
	pEctoGenesNebulosa[[gene]] <- plot_density(seurat_integrated, GeneConversion[which(GeneConversion$Gene_ID == gene), "ENSG_ID"]) + labs(title = gene)
}

dir.create(paste0(outdir_DEGs, "ectoderm_genes/"))

print("Saving Plots")
lapply(names(pEctoGenes), FUN = function(x){
	ggsave(pEctoGenes[[x]], file = paste0(outdir_DEGs, "ectoderm_genes/",x,".png"))
})

lapply(names(pEctoGenesNebulosa), FUN = function(x){
	ggsave(pEctoGenesNebulosa[[x]], file = paste0(outdir_DEGs, "ectoderm_genes/Nebulosa_",x,".png"))
})



##### Delete unused resolutions #####
for (res in c(seq(0,0.27,0.01), seq(0.29,0.5,0.01))){
	print(res)
	seurat_integrated@meta.data[,paste0("integrated_snn_res.",res)] <- NULL
}

seurat_integrated@meta.data[,"Resolution_0.08"] <- NULL
seurat_integrated@meta.data[,"Clusters_Res_0.13"] <- NULL
seurat_integrated@meta.data[,"old.ident"] <- NULL


Idents(seurat_integrated) <- "integrated_snn_res.0.28"


seurat_integrated <- SCTransform(seurat_integrated, verbose = TRUE, vars.to.regress = c("scores.G1", "scores.S", "scores.G2M", "percent.mt"), return.only.var.genes = FALSE)
saveRDS(seurat_integrated, paste0(outdir, "seurat_integrated_all_times_clustered.rds"))


#### Regress out rb % since no reason to think that growing together would result in increased rb % - likely a technical effect?
seurat_integrated <- readRDS(paste0(outdir, "seurat_integrated_all_times_clustered.rds"))
seurat_rb <- SCTransform(seurat_integrated, verbose = TRUE, vars.to.regress = c("scores.G1", "scores.S", "scores.G2M", "percent.mt", "percent.rb"), return.only.var.genes = FALSE)
saveRDS(seurat_rb, paste0(outdir, "seurat_integrated_all_times_clustered_rb_regressed.rds"))


seurat_rb <- readRDS(paste0(outdir, "seurat_integrated_all_times_clustered_rb_regressed.rds"))



