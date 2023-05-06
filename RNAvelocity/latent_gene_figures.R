library(tidyverse)
library(Seurat)
library(data.table)
library(Nebulosa)
library(viridis)
library(gridExtra)
library(ComplexHeatmap)
library('circlize')


##### Setting up Directories
dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
datadir <- paste0(dir,"output/All_data_integrated_remove_bad/")
outdir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/scVelo/velocyto/scvelo_combat_corrected_all2/"

dir.create(outdir, recursive = TRUE)


line_colors = c(FSA0006 = "#F79E29", MBE1006 = "#9B2C99", TOB0421 = "#35369C")
village_colors <- c("Uni-Culture" = "#613246", "Village" = "#A286AA")


##### Read in data #####
seurat <- readRDS(paste0(dir,"output/Distribution_tests/seurat_integrated_all_times_clustered_1pct_expressing.rds"))
meta = fread(paste0(outdir,"metadata.csv"), sep = ",")
gene_meta <- fread(paste0(outdir, "var_meta.csv"), sep = ",")

grep("C9", gene_meta[velocity_genes == "TRUE"]$Gene)

velo_meta <- meta[,c("n_unspliced_counts", "latent_time")]
rownames(velo_meta) <- meta$V1


##### Add metadata time to seurat object
seurat <- AddMetaData(seurat, velo_meta)
DefaultAssay(seurat) <- "SCT"


##### Make plots of velocity genes
#### Conversion dataframe

##### Add gene IDs for easy identification downstream #####
GeneConversion1 <- read_delim("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/DRENEA_1/outs/filtered_feature_bc_matrix/features.tsv.gz", col_names = F, delim = "\t")
GeneConversion2 <- read_delim("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/Village_B_1_week/outs/filtered_feature_bc_matrix/features.tsv.gz", col_names = F, delim = "\t")

GeneConversion <- unique(rbind(GeneConversion1, GeneConversion2))
GeneConversion <- GeneConversion[!duplicated(GeneConversion$X1),]
GeneConversion$X3 <- NULL
colnames(GeneConversion) <- c("ENSG_ID", "Gene_ID")


#### Plots
umap_out_dir <- paste0(outdir, "velocity_gene_umaps/")
dir.create(umap_out_dir, recursive = TRUE)
genes <- gene_meta[velocity_genes == "TRUE"]$Gene

for (gene in genes){
	if (gene %in% GeneConversion$Gene_ID){
		ensg <- GeneConversion$ENSG_ID[which(GeneConversion$Gene_ID == gene)]
		plot <- FeaturePlot(seurat, feature = ensg) + labs(title = gene)
		ggsave(plot, filename = paste0(umap_out_dir,gene,"_test_umap.png"))
	}
}


plot <- FeaturePlot(seurat, feature = "ENSG00000204711")
ggsave(plot, filename = paste0(outdir,"test_umap.png"))

plot2 <- plot_density(seurat, feature = "ENSG00000204711")
ggsave(plot2, filename = paste0(outdir,"test_umap_nebulosa.png"))



##### Make latent time and expression plot #####
### Remove cells without latent time
seurat_noNA <- subset(seurat, subset = latent_time >= 0)

## Make and save table for publication

heatmap_markers <- fread("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/latent_heatmap_markers.tsv", sep = "\t")

heatmap_genes <- t(data.frame(seurat_noNA[["SCT"]]@data[GeneConversion$ENSG_ID[which(GeneConversion$Gene_ID %in% heatmap_markers$Gene)],]))

fig_dt <- data.table(Barcode = colnames(seurat_noNA), 
						seurat_noNA@reductions$umap@cell.embeddings, 
						latent_time = seurat_noNA$latent_time,
						Time = seurat_noNA$Time,
						Cell_Line = seurat_noNA$Final_Assignment,
						Location = seurat_noNA$Site_rep,
						Replicate = seurat_noNA$Site_rep,
						Cryopreserved = seurat_noNA$Pool,
						ENSG00000106153 = seurat_noNA[["SCT"]]@scale.data["ENSG00000106153",])

fig_dt <- cbind(fig_dt, heatmap_genes)

fig_dt$Location <- gsub("Brisbane\\d", "Site 1", fig_dt$Location) %>%
						gsub("Melbourne\\d", "Site 2", .) %>%
						gsub("Sydney\\d", "Site 3", .) 

fig_dt$Replicate <- gsub("Brisbane", "", fig_dt$Replicate) %>%
						gsub("Melbourne", "", .) %>%
						gsub("Sydney", "", .)

fig_dt$Time <- gsub(" Day 4", "", fig_dt$Time) %>%
						gsub("Thawed Village Day 0", "Uni-culture", .) %>%
						gsub("Thawed Village Day 7", "Village", .) %>%
						gsub("Baseline", "Uni-culture", .)

fig_dt$Cryopreserved <- gsub("DRENEA_\\d", "Fresh", fig_dt$Cryopreserved) %>%
						gsub("Village_A_Baseline", "Cryopreserved", .) %>%
						gsub("Village_B_1_week", "Cryopreserved", .)


fwrite(fig_dt, paste0(outdir, "rna_velocity_figures.tsv"), sep = "\t")


pLatent <- FeaturePlot(seurat_noNA, feature = "latent_time", pt.size = 0.1, order = TRUE) +
			theme(text = element_text(size=20)) +
			scale_color_viridis(name = "Pseudo-\ntime")

legend_latent <- cowplot::get_legend(pLatent)


pLatent <- FeaturePlot(seurat_noNA, feature = "latent_time", pt.size = 0.1) + 
			scale_color_viridis(alpha = 0.25) + 
			ggtitle("RNA Velocity Pseudotime") + 
			theme(text = element_text(size=20)) +
			xlab("UMAP 1") +
			ylab("UMAP 2") +
			theme(legend.position = "none")

platent_combined  <- grid.arrange(pLatent, legend_latent, ncol = 2, widths = 2:0.5)
ggsave(platent_combined, filename = paste0(outdir, "umap_latent.png"), height = 5, width = 6.5)
ggsave(platent_combined, filename = paste0(outdir, "umap_latent.pdf"), height = 5, width = 6.5)


pPOU5F1 <- FeaturePlot(seurat_noNA, feature = GeneConversion$ENSG_ID[which(GeneConversion$Gene_ID == "POU5F1")], pt.size = 0.1) + 
			theme(text = element_text(size=20)) +
			scale_color_viridis(option = "inferno", name = "Expression")

legend_POU5F1 <- cowplot::get_legend(pPOU5F1)

pPOU5F1 <- FeaturePlot(seurat_noNA, feature = GeneConversion$ENSG_ID[which(GeneConversion$Gene_ID == "POU5F1")], pt.size = 0.1) +
			scale_color_viridis(option = "inferno", alpha = 0.25) + 
			ggtitle("POU5F1 Expression") + 
			theme(text = element_text(size=20)) +
			xlab("UMAP 1") +
			ylab("UMAP 2") +
			theme(legend.position = "none")

pPOU5F1_combined  <- grid.arrange(pPOU5F1, legend_POU5F1, ncol = 2, widths = 2:0.5)
ggsave(pPOU5F1_combined, filename = paste0(outdir, "umap_POU5F1.png"), height = 5, width = 6.5, res = 600)
ggsave(pPOU5F1_combined, filename = paste0(outdir, "umap_POU5F1.pdf"), height = 5, width = 6.5)


##### For some of the latent genes #####
genes <- c("IDO1", "IRX2","LIX1", "PTN")

for (gene in genes){
	print(gene)
	pUMAP <- FeaturePlot(seurat_noNA, feature = GeneConversion$ENSG_ID[which(GeneConversion$Gene_ID == gene)], pt.size = 0.1) + 
				theme(text = element_text(size=20)) +
				scale_color_viridis(option = "inferno", name = "Expression")

	legend_UMAP <- cowplot::get_legend(pUMAP)

	pUMAP <- FeaturePlot(seurat_noNA, feature = GeneConversion$ENSG_ID[which(GeneConversion$Gene_ID == gene)], pt.size = 0.1) +
				scale_color_viridis(option = "inferno", alpha = 0.25) + 
				ggtitle(paste0(gene," Expression") )+ 
				theme(text = element_text(size=20)) +
				xlab("UMAP 1") +
				ylab("UMAP 2") +
				theme(legend.position = "none")

	pCombined  <- grid.arrange(pUMAP, legend_UMAP, ncol = 2, widths = 2:0.5)
	ggsave(pCombined, filename = paste0(outdir, gene,"_umap.png"), height = 5, width = 6.5, dpi = 600)
	ggsave(pCombined, filename = paste0(outdir, gene,"_umap.pdf"), height = 5, width = 6.5)

}



##### Make figure with quintile on UMAP #####
seurat_noNA@meta.data$Quintile  <- with(seurat_noNA@meta.data, factor(
                            		findInterval(latent_time, c(-Inf,
                               		quantile(latent_time, probs=c(0.2, 0.4, 0.6, 0.8)), Inf)), 
                            		labels=c("Q1","Q2","Q3","Q4","Q5")))


pQuintile <- DimPlot(seurat_noNA, group.by = "Quintile", pt.size = 0.1, order = TRUE) +
			theme(text = element_text(size=20)) +
			scale_color_viridis(name = "Quintile", discrete = TRUE)

legend_quintile <- cowplot::get_legend(pQuintile)


pQuintile <- DimPlot(seurat_noNA, group.by = "Quintile", pt.size = 0.1) + 
			scale_color_viridis(alpha = 0.25, discrete = TRUE) + 
			ggtitle("Quintiles") + 
			theme(text = element_text(size=20)) +
			xlab("UMAP 1") +
			ylab("UMAP 2") +
			theme(legend.position = "none")

pQuintile_combined  <- grid.arrange(pQuintile, legend_quintile, ncol = 2, widths = 2:0.5)
ggsave(pQuintile_combined, filename = paste0(outdir, "umap_quintile.png"), height = 5, width = 6.5)
ggsave(pQuintile_combined, filename = paste0(outdir, "umap_quintile.pdf"), height = 5, width = 6.5)




seurat_noNA@meta.data$Location <- ifelse(grepl("Thawed", seurat_noNA@meta.data$Location), "Sydney-Cryopreserved", gsub("_.+", "", seurat_noNA@meta.data$Location))
seurat_noNA@meta.data$Location <- gsub("Brisbane", "Site 1", seurat_noNA@meta.data$Location) %>% gsub("Melbourne", "Site 2", .) %>% gsub("Sydney", "Site 3", .)
seurat_noNA@meta.data$Village <- gsub("Baseline", "Uni-Culture", seurat_noNA@meta.data$Time) %>% gsub("Village Day 4", "Village", .) %>% gsub("Thawed Village Day 0", "Uni-Culture", .) %>% gsub("Thawed Village Day 7", "Village", .)


seurat_noNA@meta.data$Location_Individual <- paste0(seurat_noNA@meta.data$Location, "-", seurat_noNA@meta.data$Final_Assignment)
seurat_noNA@meta.data$Location <- gsub("Site 3-Cryopreserved", "Site 3\nCryopreserved", seurat_noNA@meta.data$Location)



pLocation_Time_latent_bar <- ggplot(seurat_noNA@meta.data, aes(latent_time, fill = Village)) +
	geom_histogram(alpha = 0.75, position="identity") +
	theme_classic() +
	facet_grid(Location ~ Final_Assignment,scales = "free_y") +
	scale_fill_manual(values = village_colors) +
	scale_x_continuous(breaks=c(0,0.5,1))

ggsave(pLocation_Time_latent_bar, filename = paste0(outdir,"faceted_histogram_latent_bar.png"), width = 4.5, height = 5)
ggsave(pLocation_Time_latent_bar, filename = paste0(outdir,"faceted_histogram_latent_bar.pdf"), width = 4.5, height = 5)




# ##### Make heatmap plot ##### Did in python because better plotting
# heatmap_markers <- fread("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/latent_heatmap_markers.tsv", sep = "\t")

# latent_marker_heatmaps <- DoHeatmap(seurat, features = GeneConversion[c(GeneConversion$Gene_ID %in% heatmap_markers$Gene),]$ENSG_ID, group.bar = FALSE) + scale_fill_viridis()
# # latent_marker_heatmaps <- DoHeatmap(seurat, features = "ENSG00000204711")
# ggsave(latent_marker_heatmaps, filename = paste0(outdir, "latent_marker_heatmap.pdf"), height = 5, width = 15)

# latent_marker_heatmaps$data$Cell <- factor(latent_marker_heatmaps$data$Cell, levels = colnames(seurat)[order(seurat@meta.data$latent_time)])
# ggsave(latent_marker_heatmaps, filename = paste0(outdir, "latent_marker_heatmap_update.pdf"), height = 5, width = 15)

# ### Set up matrix ###
# # seurat_sub <- subset(seurat, cell = colnames(seurat)[order(seurat@meta.data$latent_time)])


# seurat@meta.data$latent_time <- factor(seurat@meta.data$latent_time, levels = colnames(seurat)[order(seurat@meta.data$latent_time))

# mat <- seurat[["SCT"]]@data[GeneConversion[c(GeneConversion$Gene_ID %in% heatmap_markers$Gene),]$ENSG_ID,]

# mat_order <- matrix(mat[,order(seurat@meta.data$latent_time, decreasing = TRUE)])


# ### Set up color ###
# col_fun = colorRamp2(c(-2, 0, 2), c("#440154", "#21908C", "#FDE725"))


# ht = Heatmap(mat_order, cluster_columns = FALSE, col = col_fun)

# pdf(paste0(outdir, "latent_markers_heatmap.pdf"))
# draw(ht)
# dev.off()