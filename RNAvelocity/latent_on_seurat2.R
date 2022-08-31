library(data.table)
library(Seurat)
library(ggplot2)
library(tidyverse)
# library(Nebulosa)
# library(schex)
# library(R.utils)
# library(Hmisc)


##### Setting up Directories
dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
datadir <- paste0(dir,"output/All_data_integrated_remove_bad/")
outdir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/scVelo/velocyto/scvelo_umap2/"

dir.create(outdir, recursive = TRUE)


line_colors = c(FSA0006 = "#F79E29", MBE1006 = "#9B2C99", TOB0421 = "#35369C")
village_colors <- c("Uni-Culture" = "#613246", "Village" = "#A286AA")


##### Read in data #####
seurat = readRDS(paste0(datadir,"seurat_integrated_all_times_clustered.rds"))
gene_df = fread("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/scVelo/velocyto/scvelo_umap/velocity_genes.csv")



latent_files <- list.files(outdir, pattern = ".csv")

meta_list <- lapply(latent_files, function(x){
    fread(paste0(outdir,x), sep = ",")
})


meta_df <- do.call(rbind, meta_list)


meta_latent <- data.table(meta_df[,latent_time])
colnames(meta_latent) <- "latent_time"
rownames(meta_latent) <- meta_df[,V1]

seurat <- AddMetaData(seurat, meta_latent)


### Make a reference from ENSG IDs to Gene IDs ###
GeneConversion1 <- fread("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/DRENEA_1/outs/filtered_feature_bc_matrix/features.tsv.gz", header = F, sep = "\t")
GeneConversion2 <- fread("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/Village_A_Baseline/Village_A_Baseline/outs/filtered_feature_bc_matrix/features.tsv.gz", header = F, sep = "\t")
GeneConversion <- unique(rbind(GeneConversion1, GeneConversion2))
GeneConversion$V3 <- NULL
colnames(GeneConversion) <- c("ENSG_ID", "Gene_ID")

gene_df <- GeneConversion[gene_df, on = c("Gene_ID" = "Gene")]


##### Try making new UMAP with genes used to determine latency #####
DefaultAssay(seurat) <- "integrated"
seurat_new <- ScaleData(seurat, verbose = FALSE, features = GeneConversion$ENSG_ID)
seurat_new <- RunPCA(seurat_new, npcs = 30, verbose = FALSE)
seurat_new <- RunUMAP(seurat_new, reduction = "pca", dims = 1:30)
seurat_new <- FindNeighbors(seurat_new, reduction = "pca", dims = 1:30)
seurat_new <- FindClusters(seurat_new, resolution = 0.5)


pLatent_new <- FeaturePlot(seurat_new, feature = "latent_time") +
            scale_colour_gradient2(low = "#F5AB29", mid = "#8A2F63", high = "#52203D", midpoint = 0.5)
ggsave(pLatent_new, filename = paste0(outdir, "new_latent_umap.png"))

### Try with schex
seurat_new_noNA <- subset(seurat_new, subset = latent_time >= 0)

seurat_new_noNA <- make_hexbin(seurat_new_noNA, nbins = 50, 
    dimension_reduction = "UMAP")

pLatent_umap_hex_new  <- plot_hexbin_meta(seurat_new_noNA, col="latent_time", action="median")
ggsave(pLatent_umap_hex, filename = paste0(outdir, "latent_umap_schex.png"))



latent_hist <- ggplot(seurat_new_noNA@meta.data, aes(latent_time)) +
				geom_density() +
				theme_classic()
ggsave(latent_hist, filename = paste0(outdir, "latent_histogram.png"))



pLatent <- FeaturePlot(seurat, feature = "latent_time") +
            scale_colour_gradient2(low = "#F5AB29", mid = "#8A2F63", high = "#52203D", midpoint = 0.5)
ggsave(pLatent, filename = paste0(outdir, "latent_umap.png"))


seurat_noNA <- subset(seurat, subset = latent_time >= 0)


pLatent_umap <- plot_density(seurat_noNA, "latent_time")
ggsave(pLatent_umap, filename = paste0(outdir, "latent_umap_nebulosa.png"))

### Try with schex
seurat_noNA <- make_hexbin(seurat_noNA, nbins = 50, 
    dimension_reduction = "UMAP")

pLatent_umap_hex  <- plot_hexbin_meta(seurat_noNA, col="latent_time", action="median")
ggsave(pLatent_umap_hex, filename = paste0(outdir, "latent_umap_schex.png"))



df <- data.table("latent_time" = seurat_noNA@meta.data$latent_time, 
                    "POU5F1_integrated" = seurat_noNA[["integrated"]]["ENSG00000204531",][1,],
                    "POU5F1_SCT" = seurat_noNA[["SCT"]]@scale.data["ENSG00000204531",],
                    "POU5F1_SCT_data" = seurat_noNA[["SCT"]]@data["ENSG00000204531",],
                    "POU5F1_RNA" = seurat_noNA[["RNA"]]@data["ENSG00000204531",])

pPseudotime <- ggplot(df, aes(latent_time, POU5F1_integrated)) +
    geom_point(alpha = 0.25) +
    theme_classic()

ggsave(pPseudotime, filename = paste0(outdir,"pseudotime_POU5F1_integrated.png"))


pPseudotime_SCT <- ggplot(df, aes(latent_time, POU5F1_SCT)) +
    geom_point(alpha = 0.25) +
    theme_classic()
    
ggsave(pPseudotime_SCT, filename = paste0(outdir,"pseudotime_POU5F1_SCT.png"))

pPseudotime_SCT_data <- ggplot(df, aes(latent_time, POU5F1_SCT_data)) +
    geom_point(alpha = 0.25) +
    theme_classic()
    
ggsave(pPseudotime_SCT_data, filename = paste0(outdir,"pseudotime_POU5F1_SCT_data.png"))


pPseudotime_RNA <- ggplot(df, aes(latent_time, POU5F1_RNA)) +
    geom_point(alpha = 0.25) +
    theme_classic()
    
ggsave(pPseudotime_RNA, filename = paste0(outdir,"pseudotime_POU5F1_RNA.png"))



#### Do the same for the top variable features ####
outdir_figs <- paste0(outdir,"latent_vs_variable_expression/")
dir.create(outdir_figs)



for (gene in gene_df[51:100,]$Gene){
	if (gene %in% GeneConversion$Gene_ID){
		ensg <- GeneConversion[Gene_ID == gene,]$ENSG_ID
		if (gene %in% rownames(seurat_noNA[["integrated"]])){
			df <- data.table("latent_time" = seurat_noNA@meta.data$latent_time, 
								"integrated" = seurat_noNA[["integrated"]][ensg,][1,],
								"SCT" = seurat_noNA[["SCT"]]@scale.data[ensg,],
								"SCT_data" = seurat_noNA[["SCT"]]@data[ensg,],
								"RNA" = seurat_noNA[["RNA"]]@data[ensg,],
								"Individual" = seurat_noNA@meta.data$Final_Assignment)

		pPseudotime <- ggplot(df, aes(latent_time, integrated, color = Individual)) +
			geom_point(alpha = 0.25) +
			theme_classic() +
			# facet_wrap(vars(Individual), ncol = 1) +
			geom_smooth(aes(group = 1)) +
			scale_color_manual(values = line_colors)

		ggsave(pPseudotime, filename = paste0(outdir_figs,gene,"_pseudotime_integrated.png"), width = 7, height = 6)

		} else {
			df <- data.table("latent_time" = seurat_noNA@meta.data$latent_time, 
								"SCT" = seurat_noNA[["SCT"]]@scale.data[ensg,],
								"SCT_data" = seurat_noNA[["SCT"]]@data[ensg,],
								"RNA" = seurat_noNA[["RNA"]]@data[ensg,],
								"Individual" = seurat_noNA@meta.data$Final_Assignment)
		}


		pPseudotime_SCT <- ggplot(df, aes(latent_time, SCT, color = Individual)) +
			geom_point(alpha = 0.25) +
			theme_classic() +
			# facet_wrap(vars(Individual), ncol = 1) +
			geom_smooth(aes(group = 1)) +
			scale_color_manual(values = line_colors)
			
		ggsave(pPseudotime_SCT, filename = paste0(outdir_figs,gene,"pseudotime_SCT.png"), width = 7, height = 12)

		pPseudotime_SCT_data <- ggplot(df, aes(latent_time, SCT_data, color = Individual)) +
			geom_point(alpha = 0.25) +
			theme_classic() +
			# facet_wrap(vars(Individual), ncol = 1) +
			geom_smooth(aes(group = 1)) +
			scale_color_manual(values = line_colors)
			
		ggsave(pPseudotime_SCT_data, filename = paste0(outdir_figs,gene,"_pseudotime_SCT_data.png"), width = 7, height = 12)


		pPseudotime_RNA <- ggplot(df, aes(latent_time, RNA, color = Individual)) +
			geom_point(alpha = 0.25) +
			theme_classic() +
			# facet_wrap(vars(Individual), ncol = 1) +
			geom_smooth(aes(group = 1)) +
			scale_color_manual(values = line_colors)
			
		ggsave(pPseudotime_RNA, filename = paste0(outdir_figs,gene,"_pseudotime_RNA.png"), width = 7, height = 12)
	}
}













### Plot PC1 vs PC2 colored by latent time
pPCA <- DimPlot(seurat_noNA, reduction = "pca", group.by = "latent_time")


df <- data.table(Embeddings(seurat_noNA)[,c("PC_1", "PC_2")], seurat_noNA@meta.data$latent_time)


pPCA <- ggplot(df, aes(PC_1, PC_2, color = V2)) +
            geom_point() +
            theme_classic()

ggsave(pPCA, filename = paste0(outdir_figs,"PCA_latent.png"), width = 7, height = 12)





pHisto <- ggplot(seurat_noNA@meta.data, aes(latent_time, color = Final_Assignment)) +
            geom_histogram() +
            theme_classic() +
            scale_color_manual(values = line_colors)




NANOG_POU5F1 <- data.table("NANOG_SCT" = seurat[["SCT"]]@scale.data["ENSG00000111704",],
                            "POU5F1_SCT" = seurat[["SCT"]]@scale.data["ENSG00000204531",],
                            "latent" = seurat@meta.data$latent_time)

pNANOG_POU5F1 <- ggplot(NANOG_POU5F1, aes(NANOG_SCT, POU5F1_SCT, color = latent)) +
                        geom_point(alpha = 0.5) +
                        theme_classic() +
                        scale_colour_gradient2(low = "#F5AB29", mid = "#8A2F63", high = "#52203D", midpoint = 0.5)
ggsave(pNANOG_POU5F1, filename = paste0(outdir_figs, "NANOG_POU5F1_scatter_latent.png"))



##### Test different stem cell markers on the latent time #####
primed <- c("ENSG00000109956",
"ENSG00000272398",
"ENSG00000157404",
"ENSG00000198053",
"ENSG00000154096")

naive <- c("ENSG00000128274",
"ENSG00000173762",
"ENSG00000079385",
"ENSG00000143226",
"ENSG00000134352",
"ENSG00000005893",
"ENSG00000073849")

shared <- c("ENSG00000177697",
"ENSG00000117335",
"ENSG00000076706",
"ENSG00000162493")


##### Add gene IDs for easy identification downstream #####
GeneConversion1 <- fread("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/DRENEA_1/outs/filtered_feature_bc_matrix/features.tsv.gz", header = F, sep = "\t")
GeneConversion2 <- fread("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/Village_B_1_week/outs/filtered_feature_bc_matrix/features.tsv.gz", header = F, sep = "\t")

GeneConversion <- unique(rbind(GeneConversion1, GeneConversion2))
GeneConversion <- GeneConversion[!duplicated(GeneConversion$V1),]
GeneConversion$V3 <- NULL
colnames(GeneConversion) <- c("ENSG_ID", "Gene_ID")


outdir_pluri <- paste0(outdir,"pluri_genes/")
dir.create(outdir_pluri)

for (gene in c(primed,naive,shared)){
    gene_id <- GeneConversion[ENSG_ID == gene, Gene_ID]
    df <- data.table("latent_time" = seurat@meta.data$latent_time, 
                        "SCT" = seurat[["SCT"]]@scale.data[gene,],
                        "RNA" = seurat[["RNA"]]@data[gene,],
                        "Individual" = seurat@meta.data$Final_Assignment)


    pPseudotime_SCT <- ggplot(df, aes(latent_time, SCT, color = Individual)) +
        geom_point(alpha = 0.25) +
        theme_classic() +
        # facet_wrap(vars(Individual), ncol = 1) +
        geom_smooth(aes(group = 1)) +
        scale_color_manual(values = line_colors)
        
    ggsave(pPseudotime_SCT, filename = paste0(outdir_pluri,gene_id,"pseudotime_SCT.png"), width = 7, height = 12)


    pPseudotime_RNA <- ggplot(df, aes(latent_time, RNA, color = Individual)) +
        geom_point(alpha = 0.25) +
        theme_classic() +
        # facet_wrap(vars(Individual), ncol = 1) +
        geom_smooth(aes(group = 1)) +
        scale_color_manual(values = line_colors)
    ggsave(pPseudotime_RNA, filename = paste0(outdir_pluri,gene_id,"_pseudotime_RNA.png"), width = 7, height = 12)
        

}


mean(seurat@meta.data$latent_time[which(seurat[["RNA"]]@data["ENSG00000204531",] == 0 & seurat[["RNA"]]@data["ENSG00000154096",] == 0)], na.rm = TRUE)
mean(seurat@meta.data$latent_time[which(seurat[["RNA"]]@data["ENSG00000204531",] > 0 & seurat[["RNA"]]@data["ENSG00000154096",] > 0)], na.rm = TRUE)
mean(seurat@meta.data$latent_time[which(seurat[["RNA"]]@data["ENSG00000204531",] > 0 & seurat[["RNA"]]@data["ENSG00000154096",] == 0)], na.rm = TRUE)
mean(seurat@meta.data$latent_time[which(seurat[["RNA"]]@data["ENSG00000204531",] == 0 & seurat[["RNA"]]@data["ENSG00000154096",] > 0)], na.rm = TRUE)



table(seurat[["RNA"]]@data["ENSG00000204531",] == 0, seurat[["RNA"]]@data["ENSG00000154096",] == 0)



THY1_POU5F1 <- data.table("THY1_SCT" = seurat[["RNA"]]@data["ENSG00000154096",],
                            "POU5F1_SCT" = seurat[["RNA"]]@data["ENSG00000204531",],
                            "latent" = seurat@meta.data$latent_time)

pTHY1_POU5F1 <- ggplot(THY1_POU5F1, aes(THY1_SCT, POU5F1_SCT, color = latent)) +
                        geom_point(alpha = 0.5) +
                        theme_classic() +
                        scale_colour_gradient2(low = "#F5AB29", mid = "#8A2F63", high = "#52203D", midpoint = 0.5)
ggsave(pTHY1_POU5F1, filename = paste0(outdir_figs, "THY1_POU5F1_scatter_latent.png"))






CCNB1_POU5F1 <- data.table("CCNB1_SCT" = seurat[["SCT"]]@scale.data["ENSG00000134057",],
                            "POU5F1_SCT" = seurat[["SCT"]]@scale.data["ENSG00000204531",],
                            "latent" = seurat@meta.data$latent_time)

pCCNB1_POU5F1 <- ggplot(CCNB1_POU5F1, aes(CCNB1_SCT, POU5F1_SCT, color = latent)) +
                        geom_point(alpha = 0.5) +
                        theme_classic() +
                        scale_colour_gradient2(low = "#F5AB29", mid = "#8A2F63", high = "#52203D", midpoint = 0.5)
ggsave(pCCNB1_POU5F1, filename = paste0(outdir_figs, "CCNB1_POU5F1_scatter_latent.png"))



C9orf135_POU5F1 <- data.table("C9orf135_SCT" = seurat[["SCT"]]@scale.data["ENSG00000204711",],
                            "POU5F1_SCT" = seurat[["SCT"]]@scale.data["ENSG00000204531",],
                            "latent" = seurat@meta.data$latent_time)

pC9orf135_POU5F1 <- ggplot(C9orf135_POU5F1, aes(C9orf135_SCT, POU5F1_SCT, color = latent)) +
                        geom_point(alpha = 0.5) +
                        theme_classic() +
                        scale_colour_gradient2(low = "#F5AB29", mid = "#8A2F63", high = "#52203D", midpoint = 0.5)
ggsave(pC9orf135_POU5F1, filename = paste0(outdir_figs, "C9orf135_POU5F1_scatter_latent.png"))




##### Try splitting latent time into 10 bins #####
latent_split <- data.table(latent_time = seurat@meta.data$latent_time, 
							equal_groups = as.numeric(cut2(seurat@meta.data$latent_time, g=10)), 
							POU5F1 = seurat[["RNA"]]@data["ENSG00000204531",],
							equal_bin_width = ifelse(seurat@meta.data$latent_time < 0.1, 1,
												ifelse(seurat@meta.data$latent_time < 0.2, 2,
													ifelse(seurat@meta.data$latent_time < 0.3, 3,
														ifelse(seurat@meta.data$latent_time < 0.4, 4,
															ifelse(seurat@meta.data$latent_time < 0.5, 5,
																ifelse(seurat@meta.data$latent_time < 0.6, 6,
																	ifelse(seurat@meta.data$latent_time < 0.7, 7,
																		ifelse(seurat@meta.data$latent_time < 0.8, 8,
																			ifelse(seurat@meta.data$latent_time < 0.9, 9,
																				ifelse(seurat@meta.data$latent_time <= 1, 10, NA)))))))))))

rownames(latent_split) <- colnames(seurat)

latent_split <- latent_split[!is.na(latent_time)]


bin_box_pou5f1 <- ggplot(latent_split, aes(factor(equal_groups), POU5F1)) +
	geom_boxplot(outlier.size = 0.5) +
	theme_classic()
ggsave(bin_box_pou5f1, filename = paste0(outdir_figs, "bin_box_pou5f1.png"))


summary(aov(POU5F1 ~ equal_groups, data = latent_split))
summary(lm(POU5F1 ~ equal_groups, data = latent_split))



##### Make histograms of each location, colored by time, faceted by location #####
seurat_noNA@meta.data$Location <- ifelse(grepl("Thawed", seurat_noNA@meta.data$Location), "Sydney-Cryopreserved", gsub("_.+", "", seurat_noNA@meta.data$Location))
seurat_noNA@meta.data$Location <- gsub("Brisbane", "Site 1", seurat_noNA@meta.data$Location) %>% gsub("Melbourne", "Site 2", .) %>% gsub("Sydney", "Site 3", .)
seurat_noNA@meta.data$Village <- gsub("Baseline", "Uni-Culture", seurat_noNA@meta.data$Time) %>% gsub("Village Day 4", "Village", .) %>% gsub("Thawed Village Day 0", "Uni-Culture", .) %>% gsub("Thawed Village Day 7", "Village", .)


seurat_noNA@meta.data$Location_Individual <- paste0(seurat_noNA@meta.data$Location, "-", seurat_noNA@meta.data$Final_Assignment)


pLocation_Time_latent <- ggplot(seurat_noNA@meta.data, aes(latent_time, fill = Village)) +
	geom_density(alpha = 0.75) +
	theme_classic() +
	facet_grid(Location ~ Final_Assignment) +
	scale_fill_manual(values = village_colors)
ggsave(pLocation_Time_latent, filename = paste0(outdir,"faceted_histogram_latent.png"))
ggsave(pLocation_Time_latent, filename = paste0(outdir,"faceted_histogram_latent.pdf"))



pLocation_Time_latent_line <- ggplot(seurat_noNA@meta.data, aes(latent_time, fill = Village)) +
	geom_density(alpha = 0.75) +
	theme_classic() +
	facet_grid(vars(Final_Assignment)) +
	scale_fill_manual(values = village_colors)
ggsave(pLocation_Time_latent_line, filename = paste0(outdir,"faceted_histogram_latent_line.png"), width = 4, height = 3)
ggsave(pLocation_Time_latent_line, filename = paste0(outdir,"faceted_histogram_latent_line.pdf"), width = 4, height = 3)

