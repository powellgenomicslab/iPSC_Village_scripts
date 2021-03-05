##### Read in Arguments #####
print("Reading and assigning input arguments")

args <- commandArgs(TRUE)
arguments <- read.table(args, header = F)
stage <- arguments[1,]
outdir <- arguments[2,]
dir <- arguments[3,]
sge <- arguments[4,]
print(sge)
sge <- as.numeric(as.character(sge))
print(sge)
print(outdir)



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



print(sessionInfo())

### Set colors
cell_line_colors <- c(ppalette[c(2,4,5)], c(gpalette[c(2,3,4)]))
names(cell_line_colors) <- c("FSA0006","MBE1006","TOB0421","doublet","unassigned","combination_fail_remove")

PoolColors <- brewer.pal(n = 8, name = "Dark2")

### Set names ###
pools <- dir("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/")
old_name <- c("22_FSA", "29_MBE", "36_TOB00421_i_E8")
new_name <- c("FSA0006","MBE1006","TOB0421")
cell_key <- data.frame(old_name, new_name)


# dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
# outbase <- paste0(dir,"output/Brisbane_Clustering_Testing/")
# dir.create(outbase)


# stage = "SCTnormalization_no_covatiates"
# outdir <- paste0(outbase,stage,"/")
# dir.create(outdir)


##### Read in the data that has been filtered (only filtered on mt %)
seurat <- readRDS(paste0(dir,"output/Seurat_w_Hash_freeze_thaw/CellCycle/Brisbane_Baseline_seurat.rds"))


if (stage == "SCTnormalization_no_covatiates"){
    seurat <- SCTransform(seurat, verbose = TRUE, return.only.var.genes = FALSE)
    seurat <- RunPCA(seurat, npcs = 100)
    seurat <- RunUMAP(seurat, dims = 1:100, verbose = TRUE)
    seurat <- FindNeighbors(seurat, dims = 1:100, verbose = TRUE)


    saveRDS(seurat, paste0(outdir, "Brisbane_Baseline_seurat_SCT_noCov.rds"))
    # seurat <- readRDS(paste0(outdir, "seurat_SCT_noCov.rds"))




    plots <- DimPlot(seurat, group.by = c("Time", "Final_Assignment"), combine = FALSE)
    plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 4, 
        byrow = TRUE, override.aes = list(size = 2.5))))
    ggsave(CombinePlots(plots, ncol = 3), filename = paste0(outdir,"umap.png"),width = 20, height = 9)


    umap <- DimPlot(seurat, reduction = "umap", group.by = "Pool", pt.size = .01, split.by = 'Pool')
    ggsave(umap, filename = paste0(outdir, "umap_pool.png"), height = 4, width = 10)

    umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .01, split.by = 'Time')
    ggsave(umap, filename = paste0(outdir, "umap_Time.png"), height = 4, width = 10)

    umap <- DimPlot(seurat, reduction = "umap", group.by = "MULTI_ID", pt.size = .01, split.by = 'MULTI_ID')
    ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate.png"), height = 4, width = 10)

    umap <- DimPlot(seurat, reduction = "umap", group.by = "Final_Assignment", pt.size = .01, split.by = 'Final_Assignment')
    ggsave(umap, filename = paste0(outdir, "umap_CellLine.png"), height = 4, width = 10)

    umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'MULTI_ID')
    ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate_Time.png"), height = 4, width = 10)

    umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'Final_Assignment')
    ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time.png"), height = 4, width = 10)

    seurat$Time_Final_Assignment <- paste0(seurat$Time, seurat$Final_Assignment)
    umap <- DimPlot(seurat, reduction = "umap", group.by = "Time_Final_Assignment", pt.size = .01, split.by = 'MULTI_ID')
    ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time_Site.png"), height = 4, width = 10)

    umap <- DimPlot(seurat,  pt.size = .01, group.by = c("phases"))  
    ggsave(umap, filename = paste0(outdir, "umap_CellCycle.png"))



} else if (stage == "SCTnormalization_individual"){
    seurat <- SCTransform(seurat, verbose = TRUE, return.only.var.genes = FALSE, vars.to.regress = "Final_Assignment")
    seurat <- RunPCA(seurat, npcs = 30)
    seurat <- RunUMAP(seurat, dims = 1:30, verbose = TRUE)
    seurat <- FindNeighbors(seurat, dims = 1:30, verbose = TRUE)


    saveRDS(seurat, paste0(outdir, "Brisbane_Baseline_seurat_SCT_Individual_cov.rds"))
    # seurat <- readRDS(paste0(outdir, "Brisbane_Baseline_seurat_SCT_Individual_cov.rds"))




    plots <- DimPlot(seurat, group.by = c("Final_Assignment"))
    ggsave(plots, filename = paste0(outdir,"umap.png"),width = 20, height = 9)

    umap <- DimPlot(seurat, reduction = "umap", group.by = "Pool", pt.size = .01, split.by = 'Pool')
    ggsave(umap, filename = paste0(outdir, "umap_pool.png"), height = 4, width = 10)

    umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .01, split.by = 'Time')
    ggsave(umap, filename = paste0(outdir, "umap_Time.png"), height = 4, width = 10)

    umap <- DimPlot(seurat, reduction = "umap", group.by = "MULTI_ID", pt.size = .01, split.by = 'MULTI_ID')
    ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate.png"), height = 4, width = 10)

    umap <- DimPlot(seurat, reduction = "umap", group.by = "Final_Assignment", pt.size = .01, split.by = 'Final_Assignment')
    ggsave(umap, filename = paste0(outdir, "umap_CellLine.png"), height = 4, width = 10)

    umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'MULTI_ID')
    ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate_Time.png"), height = 4, width = 10)

    umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'Final_Assignment')
    ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time.png"), height = 4, width = 10)

    seurat$Time_Final_Assignment <- paste0(seurat$Time, seurat$Final_Assignment)
    umap <- DimPlot(seurat, reduction = "umap", group.by = "Time_Final_Assignment", pt.size = .01, split.by = 'MULTI_ID')
    ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time_Site.png"), height = 4, width = 10)

    umap <- DimPlot(seurat,  pt.size = .01, group.by = c("phases"))  
    ggsave(umap, filename = paste0(outdir, "umap_CellCycle.png"))

    mt_umap <- FeaturePlot(seurat, features = "percent.mt")
    ggsave(mt_umap, filename = paste0(outdir,"mt_percent_umap_seurat.png"))

    rb_umap <- FeaturePlot(seurat, features = "percent.rb")
    ggsave(rb_umap, filename = paste0(outdir,"rb_percent_umap_seurat.png"))

    libsize_umap <- FeaturePlot(seurat, features = "nCount_RNA")
    ggsave(libsize_umap, filename = paste0(outdir,"libsize_umap_seurat.png"))

    Ngene_umap <- FeaturePlot(seurat, features = "nFeature_RNA")
    ggsave(Ngene_umap, filename = paste0(outdir,"Ngene_umap_seurat.png"))





} else if (stage == "SCTnormalization_CellCycle"){
    seurat <- SCTransform(seurat, verbose = TRUE, return.only.var.genes = FALSE, vars.to.regress = c("scores.G1", "scores.S", "scores.G2M"))
    seurat <- RunPCA(seurat, npcs = 30)
    seurat <- RunUMAP(seurat, dims = 1:30, verbose = TRUE)
    seurat <- FindNeighbors(seurat, dims = 1:30, verbose = TRUE)


    saveRDS(seurat, paste0(outdir, "Brisbane_seurat_SCT_CellCycle_cov.rds"))
    # seurat <- readRDS(paste0(outdir, "Brisbane_seurat_SCT_CellCycle_cov.rds"))



    plots <- DimPlot(seurat, group.by = c("Final_Assignment"))
    ggsave(plots, filename = paste0(outdir,"umap.png"),width = 20, height = 9)

    umap <- DimPlot(seurat, reduction = "umap", group.by = "Pool", pt.size = .01, split.by = 'Pool')
    ggsave(umap, filename = paste0(outdir, "umap_pool.png"), height = 4, width = 10)

    umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .01, split.by = 'Time')
    ggsave(umap, filename = paste0(outdir, "umap_Time.png"), height = 4, width = 10)

    umap <- DimPlot(seurat, reduction = "umap", group.by = "MULTI_ID", pt.size = .01, split.by = 'MULTI_ID')
    ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate.png"), height = 4, width = 10)

    umap <- DimPlot(seurat, reduction = "umap", group.by = "Final_Assignment", pt.size = .01, split.by = 'Final_Assignment')
    ggsave(umap, filename = paste0(outdir, "umap_CellLine.png"), height = 4, width = 10)

    umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'MULTI_ID')
    ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate_Time.png"), height = 4, width = 10)

    umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'Final_Assignment')
    ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time.png"), height = 4, width = 10)

    seurat$Time_Final_Assignment <- paste0(seurat$Time, seurat$Final_Assignment)
    umap <- DimPlot(seurat, reduction = "umap", group.by = "Time_Final_Assignment", pt.size = .01, split.by = 'MULTI_ID')
    ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time_Site.png"), height = 4, width = 10)

    umap <- DimPlot(seurat,  pt.size = .01, group.by = c("Phase"))  
    ggsave(umap, filename = paste0(outdir, "umap_CellCycle.png"))

    mt_umap <- FeaturePlot(seurat, features = "percent.mt")
    ggsave(mt_umap, filename = paste0(outdir,"mt_percent_umap_seurat.png"))

    rb_umap <- FeaturePlot(seurat, features = "percent.rb")
    ggsave(rb_umap, filename = paste0(outdir,"rb_percent_umap_seurat.png"))

    libsize_umap <- FeaturePlot(seurat, features = "nCount_RNA")
    ggsave(libsize_umap, filename = paste0(outdir,"libsize_umap_seurat.png"))

    Ngene_umap <- FeaturePlot(seurat, features = "nFeature_RNA")
    ggsave(Ngene_umap, filename = paste0(outdir,"Ngene_umap_seurat.png"))




} else if (stage == "SCTnormalization_CellCycle_AllTime"){

    seurat_4d <- readRDS(paste0(dir,"output/Seurat_w_Hash_freeze_thaw/CellCycle/Brisbane_Village_Day_4_seurat.rds"))
    seurat <- merge( seurat, seurat_4d, project = "Brisbane_clustering")

    seurat <- SCTransform(seurat, verbose = TRUE, return.only.var.genes = FALSE, vars.to.regress = c("scores.G1", "scores.S", "scores.G2M"))
    seurat <- RunPCA(seurat, npcs = 30)
    seurat <- RunUMAP(seurat, dims = 1:30, verbose = TRUE)
    seurat <- FindNeighbors(seurat, dims = 1:30, verbose = TRUE)


    saveRDS(seurat, paste0(outdir, "Brisbane_seurat_SCT_CellCycle_cov_allTime.rds"))
    # seurat <- readRDS(paste0(outdir, "Brisbane_seurat_SCT_CellCycle_cov_allTime.rds"))



    plots <- DimPlot(seurat, group.by = c("Final_Assignment"))
    ggsave(plots, filename = paste0(outdir,"umap.png"),width = 20, height = 9)

    umap <- DimPlot(seurat, reduction = "umap", group.by = "Pool", pt.size = .01, split.by = 'Pool')
    ggsave(umap, filename = paste0(outdir, "umap_pool.png"), height = 4, width = 10)

    umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .01, split.by = 'Time')
    ggsave(umap, filename = paste0(outdir, "umap_Time.png"), height = 4, width = 10)

    umap <- DimPlot(seurat, reduction = "umap", group.by = "MULTI_ID", pt.size = .01, split.by = 'MULTI_ID')
    ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate.png"), height = 4, width = 10)

    umap <- DimPlot(seurat, reduction = "umap", group.by = "Final_Assignment", pt.size = .01, split.by = 'Final_Assignment')
    ggsave(umap, filename = paste0(outdir, "umap_CellLine.png"), height = 4, width = 10)

    umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'MULTI_ID')
    ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate_Time.png"), height = 4, width = 10)

    umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'Final_Assignment')
    ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time.png"), height = 4, width = 10)

    seurat$Time_Final_Assignment <- paste0(seurat$Time, seurat$Final_Assignment)
    umap <- DimPlot(seurat, reduction = "umap", group.by = "Time_Final_Assignment", pt.size = .01, split.by = 'MULTI_ID')
    ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time_Site.png"), height = 4, width = 10)

    umap <- DimPlot(seurat,  pt.size = .01, group.by = c("phases"))  
    ggsave(umap, filename = paste0(outdir, "umap_CellCycle.png"))

    mt_umap <- FeaturePlot(seurat, features = "percent.mt")
    ggsave(mt_umap, filename = paste0(outdir,"mt_percent_umap_seurat.png"))

    rb_umap <- FeaturePlot(seurat, features = "percent.rb")
    ggsave(rb_umap, filename = paste0(outdir,"rb_percent_umap_seurat.png"))

    libsize_umap <- FeaturePlot(seurat, features = "nCount_RNA")
    ggsave(libsize_umap, filename = paste0(outdir,"libsize_umap_seurat.png"))

    Ngene_umap <- FeaturePlot(seurat, features = "nFeature_RNA")
    ggsave(Ngene_umap, filename = paste0(outdir,"Ngene_umap_seurat.png"))





} else if (stage == "SCTnormalization_CellCycleScores_AllTime"){

    seurat_4d <- readRDS(paste0(dir,"output/Seurat_w_Hash_freeze_thaw/CellCycle/Brisbane_Village_Day_4_seurat.rds"))
    seurat <- merge( seurat, seurat_4d, project = "Brisbane_clustering")

    seurat <- SCTransform(seurat, verbose = TRUE, return.only.var.genes = FALSE, vars.to.regress = c("scores.G1", "scores.S", "scores.G2M"))
    seurat <- RunPCA(seurat, npcs = 30)
    seurat <- RunUMAP(seurat, dims = 1:30, verbose = TRUE)
    seurat <- FindNeighbors(seurat, dims = 1:30, verbose = TRUE)


    saveRDS(seurat, paste0(outdir, "Brisbane_seurat_SCT_CellCycle_cov_allTime.rds"))
    # seurat <- readRDS(paste0(outdir, "Brisbane_seurat_SCT_CellCycle_cov_allTime.rds"))



    plots <- DimPlot(seurat, group.by = c("Final_Assignment"))
    ggsave(plots, filename = paste0(outdir,"umap.png"),width = 20, height = 9)

    umap <- DimPlot(seurat, reduction = "umap", group.by = "Pool", pt.size = .01, split.by = 'Pool')
    ggsave(umap, filename = paste0(outdir, "umap_pool.png"), height = 4, width = 10)

    umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .01, split.by = 'Time')
    ggsave(umap, filename = paste0(outdir, "umap_Time.png"), height = 4, width = 10)

    umap <- DimPlot(seurat, reduction = "umap", group.by = "MULTI_ID", pt.size = .01, split.by = 'MULTI_ID')
    ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate.png"), height = 4, width = 10)

    umap <- DimPlot(seurat, reduction = "umap", group.by = "Final_Assignment", pt.size = .01, split.by = 'Final_Assignment')
    ggsave(umap, filename = paste0(outdir, "umap_CellLine.png"), height = 4, width = 10)

    umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'MULTI_ID')
    ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate_Time.png"), height = 4, width = 10)

    umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'Final_Assignment')
    ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time.png"), height = 4, width = 10)

    seurat$Time_Final_Assignment <- paste0(seurat$Time, seurat$Final_Assignment)
    umap <- DimPlot(seurat, reduction = "umap", group.by = "Time_Final_Assignment", pt.size = .01, split.by = 'MULTI_ID')
    ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time_Site.png"), height = 4, width = 10)

    umap <- DimPlot(seurat,  pt.size = .01, group.by = c("phases"))  
    ggsave(umap, filename = paste0(outdir, "umap_CellCycle.png"))

    mt_umap <- FeaturePlot(seurat, features = "percent.mt")
    ggsave(mt_umap, filename = paste0(outdir,"mt_percent_umap_seurat.png"))

    rb_umap <- FeaturePlot(seurat, features = "percent.rb")
    ggsave(rb_umap, filename = paste0(outdir,"rb_percent_umap_seurat.png"))

    libsize_umap <- FeaturePlot(seurat, features = "nCount_RNA")
    ggsave(libsize_umap, filename = paste0(outdir,"libsize_umap_seurat.png"))

    Ngene_umap <- FeaturePlot(seurat, features = "nFeature_RNA")
    ggsave(Ngene_umap, filename = paste0(outdir,"Ngene_umap_seurat.png"))


} else if (stage == "SCTnormalization_CellCycle_integrated_Time"){
    seurat_4d <- readRDS(paste0(dir,"output/Seurat_w_Hash_freeze_thaw/CellCycle/Brisbane_Village_Day_4_seurat.rds"))

    seurat_list <- list()
        seurat_list[[paste0("Baseline")]] <- seurat
        seurat_list[[paste0("Day4")]] <- seurat_4d
    
    message("The seurat objects for integration:")
    print(seurat_list)

    seurat_list <- lapply(seurat_list, function(x) {
        SCTransform(x, verbose = TRUE, vars.to.regress = c("scores.G1", "scores.S", "scores.G2M"))
    })


    seurat_features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)
    seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = seurat_features)
    seurat_anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT", 
        anchor.features = seurat_features)
    seurat_integrated <- IntegrateData(anchorset = seurat_anchors, normalization.method = "SCT")

    seurat_integrated <- RunPCA(object = seurat_integrated, verbose = FALSE)
    seurat_integrated <- RunUMAP(object = seurat_integrated, dims = 1:30)

    seurat_integrated <- FindNeighbors(seurat_integrated, reduction = "pca", assay = "integrated", dims = 1:20)

    saveRDS(seurat_integrated, paste0(outdir, "Brisbane_seurat_SCT_CellCycle_cov_time_integration.rds"))
    # seurat <- readRDS(paste0(outdir, "Brisbane_seurat_SCT_CellCycle_cov_Indiv_integration.rds"))


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
    ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time_Site.png"), height = 4, width = 10)

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



} else if (stage == "SCTnormalization_CellCycle_integrated_Time_DE"){
    seurat <- readRDS(paste0(dir, "output/Brisbane_Clustering_Testing/SCTnormalization_CellCycle_integrated_Time/Brisbane_seurat_SCT_CellCycle_cov_time_integration.rds"))

    Idents(seurat) <- "Final_Assignment"

    degs <- FindAllMarkers(seurat, assay = "RNA")

    ##### Add gene IDs for easy identification downstream #####
    GeneConversion1 <- read_delim("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/DRENEA_1/outs/filtered_feature_bc_matrix/features.tsv.gz", col_names = F, delim = "\t")
    GeneConversion2 <- read_delim("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/Village_B_1_week/outs/filtered_feature_bc_matrix/features.tsv.gz", col_names = F, delim = "\t")

    GeneConversion <- unique(rbind(GeneConversion1, GeneConversion2))
    GeneConversion <- GeneConversion[!duplicated(GeneConversion$X1),]
    GeneConversion <- left_join(data.frame(rownames(seurat)), GeneConversion, by= c("rownames.seurat." = "X1"))
    GeneConversion$X3 <- NULL
    colnames(GeneConversion) <- c("ENSG_ID", "Gene_ID")


    ##### Add Gene IDs #####
    degs <- left_join(degs, GeneConversion, by = c("gene" = "ENSG_ID"))


    TopMarkers <- degs %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
    print(dim(TopMarkers))
    print(head(TopMarkers))

    print("Making Plots")
    pTopDEgenes <- list()
    for (clust in unique(degs$cluster)){
        print(clust)
        for (gene in TopMarkers$gene[which(TopMarkers$cluster == clust)]){
        # print(TopMarkers$gene[which(TopMarkers$cluster == clust)])
            pTopDEgenes[[paste0("Cluster",clust,"_",unique(TopMarkers[which(TopMarkers$gene == gene), "Gene_ID"]))]] <- FeaturePlot(seurat, features = gene) + labs(title = unique(TopMarkers[which(TopMarkers$gene == gene), "Gene_ID"]))
        }
    }

    print("Saving Plots")
    lapply(names(pTopDEgenes), FUN = function(x){
        ggsave(pTopDEgenes[[x]], file = paste0(outdir,"TopDEgenesCluster",x,".png"))
    })




} else if (stage=="SCTnormalization_CellCycle_integrated_Time_Multi_Resolution"){
    seurat <- readRDS(paste0(dir, "output/Brisbane_Clustering_Testing/SCTnormalization_CellCycle_integrated_Time/Brisbane_seurat_SCT_CellCycle_cov_time_integration.rds"))

    ### Find the clusters for this resolution ###
    res <- sge*0.01-0.01
    print(paste0("Finding Clusters for resolution=",res))
    seurat <- FindClusters(seurat, resolution = res)

    ### Plot the UMAP for this resolution ###
    UMAP <- DimPlot(seurat, reduction = "umap")
    ggsave(file = paste0(outdir,"UMAPresolution",res,".png"), plot = UMAP)

    ### Save the identities for this resolution ###
    identities <- as.data.frame(Idents(seurat))
    colnames(identities) <- c(paste0("Resolution_",res))

    write_delim(identities, paste0(outdir,"Identities_Resolution_",res,".txt"), delim = "\t")



} else if (stage=="SCTnormalization_CellCycle_integrated_Time_Clustree"){
    ResolutionsFileList <- list.files(paste0(dir,"output/Brisbane_Clustering_Testing/SCTnormalization_CellCycle_integrated_Time_Multi_Resolution/"), pattern = "Identities_Resolution_")

    ResolutionsList <- lapply(ResolutionsFileList, FUN = function(x){
        read_delim(paste0(dir,"output/Brisbane_Clustering_Testing/SCTnormalization_CellCycle_integrated_Time_Multi_Resolution/",x), delim = "\t")
    })
    
    ResolutionsDF <- do.call(cbind, ResolutionsList)

    ##### CLUSTREE #####
    print("Using Clustree to diagram changes with increasing resolutions")
    pClustree <- clustree(ResolutionsDF, prefix = "Resolution_")
    ggsave(file = paste0(outdir, "ClusteringClustree.png"), plot = pClustree)




} else if (stage == "SCTnormalization_CellCycle_integrated_Individual"){
    seurat_4d <- readRDS(paste0(dir,"output/Seurat_w_Hash_freeze_thaw/CellCycle/Brisbane_Village_Day_4_seurat.rds"))

    seurat_list <- list()
    for (line in unique(seurat@meta.data$Final_Assignment)){
        seurat_list[[paste0("Baseline_",line)]] <- subset(seurat, subset = Final_Assignment == line)
        seurat_list[[paste0("Day4_",line)]] <- subset(seurat_4d, subset = Final_Assignment == line)
    }
    message("The seurat objects for integration:")
    print(seurat_list)

    seurat_list <- lapply(seurat_list, function(x) {
        SCTransform(x, verbose = TRUE, vars.to.regress = c("scores.G1", "scores.S", "scores.G2M"))
    })


    seurat_features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)
    seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = seurat_features)
    seurat_anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT", 
        anchor.features = seurat_features)
    seurat_integrated <- IntegrateData(anchorset = seurat_anchors, normalization.method = "SCT")

    seurat_integrated <- RunPCA(object = seurat_integrated, verbose = FALSE)
    seurat_integrated <- RunUMAP(object = seurat_integrated, dims = 1:30)

    seurat_integrated <- FindNeighbors(seurat_integrated, reduction = "pca", assay = "integrated", dims = 1:20)

    saveRDS(seurat_integrated, paste0(outdir, "Brisbane_seurat_SCT_CellCycle_cov_Indiv_integration.rds"))
    # seurat <- readRDS(paste0(outdir, "Brisbane_seurat_SCT_CellCycle_cov_Indiv_integration.rds"))


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
    ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time_Site.png"), height = 4, width = 10)

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



} else if (stage=="SCTnormalization_CellCycle_integrated_Individual_Multi_Resolution"){
    seurat <- readRDS(paste0(dir, "output/Brisbane_Clustering_Testing/SCTnormalization_CellCycle_integrated_Individual/Brisbane_seurat_SCT_CellCycle_cov_Indiv_integration.rds"))

    ### Find the clusters for this resolution ###
    res <- sge*0.01-0.01
    print(paste0("Finding Clusters for resolution=",res))
    seurat <- FindClusters(seurat, resolution = res)

    ### Plot the UMAP for this resolution ###
    UMAP <- DimPlot(seurat, reduction = "umap")
    ggsave(file = paste0(outdir,"UMAPresolution",res,".png"), plot = UMAP)

    ### Save the identities for this resolution ###
    identities <- as.data.frame(Idents(seurat))
    colnames(identities) <- c(paste0("Resolution_",res))

    write_delim(identities, paste0(outdir,"Identities_Resolution_",res,".txt"), delim = "\t")



} else if (stage=="SCTnormalization_CellCycle_integrated_Individual_Clustree"){
    ResolutionsFileList <- list.files(paste0(dir,"output/Brisbane_Clustering_Testing/SCTnormalization_CellCycle_integrated_Individual_Multi_Resolution/"), pattern = "Identities_Resolution_")

    ResolutionsList <- lapply(ResolutionsFileList, FUN = function(x){
        read_delim(paste0(dir,"output/Brisbane_Clustering_Testing/SCTnormalization_CellCycle_integrated_Individual_Multi_Resolution/",x), delim = "\t")
    })
    
    ResolutionsDF <- do.call(cbind, ResolutionsList)

    ##### CLUSTREE #####
    print("Using Clustree to diagram changes with increasing resolutions")
    pClustree <- clustree(ResolutionsDF, prefix = "Resolution_")
    ggsave(file = paste0(outdir, "ClusteringClustree.png"), plot = pClustree)






} else if (stage == "SCTnormalization_CellCycle_integrated_Individual_DE"){
    seurat <- readRDS(paste0(dir, "output/Brisbane_Clustering_Testing/SCTnormalization_CellCycle_integrated_Time/Brisbane_seurat_SCT_CellCycle_cov_Indiv_integration.rds"))
    resolutions <- c(0.02,0.04,0.08,0.09,0.12,0.14,0.2)
    
    resolutions_assignments_list <- lapply(resolutions, function(res){
        read_delim(paste0(dir,"output/Brisbane_Clustering_Testing/SCTnormalization_CellCycle_integrated_Individual_Multi_Resolution/Identities_Resolution_",res,".txt"), delim = "\t")
    })

    ResolutionsDF <- do.call(cbind, resolutions_assignments_list)

    degs_list <- lapply(resolutions, function(x){
        Idents(seurat) <- ResolutionsDF[,paste0("Resolution_", x)]
        FindAllMarkers(seurat, assay = "RNA")
    })
    names(degs_list) <- paste0("Resolution_", resolutions)

    degs_list <- lapply(names(degs_list), function(x){
        colnames(degs_list[[x]]) <- c(paste0(colnames(degs_list[[x]])[,1:(ncol(degs_list[[x]])-1)],"_",x),colnames(degs_list[[x]][,ncol(degs_list[[x]])))
        return(degs_list[[x]])
    })

    degs <- join_all(degs_list, by = "gene", type = "full")

    ##### Add gene IDs for easy identification downstream #####
    GeneConversion1 <- read_delim("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/DRENEA_1/outs/filtered_feature_bc_matrix/features.tsv.gz", col_names = F, delim = "\t")
    GeneConversion2 <- read_delim("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/Village_B_1_week/outs/filtered_feature_bc_matrix/features.tsv.gz", col_names = F, delim = "\t")

    GeneConversion <- unique(rbind(GeneConversion1, GeneConversion2))
    GeneConversion <- GeneConversion[!duplicated(GeneConversion$X1),]
    GeneConversion <- left_join(data.frame(rownames(seurat)), GeneConversion, by= c("rownames.seurat." = "X1"))
    GeneConversion$X3 <- NULL
    colnames(GeneConversion) <- c("ENSG_ID", "Gene_ID")


    ##### Add Gene IDs #####
    degs <- left_join(degs, GeneConversion, by = c("gene" = "ENSG_ID"))

    degs$count <- rowSums(!is.na(degs[,grep("", colnames(degs))]))


    TopMarkers <- degs %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
    print(dim(TopMarkers))
    print(head(TopMarkers))

    print("Making Plots")
    pTopDEgenes <- list()
    for (clust in unique(degs$cluster)){
        print(clust)
            for (gene in TopMarkers$gene[which(TopMarkers$cluster == clust)]){
            # print(TopMarkers$gene[which(TopMarkers$cluster == clust)])
                pTopDEgenes[[paste0("Cluster",clust,"_",unique(TopMarkers[which(TopMarkers$gene == gene), "Gene_ID"]))]] <- FeaturePlot(seurat, features = gene) + labs(title = unique(TopMarkers[which(TopMarkers$gene == gene), "Gene_ID"]))
            }
        # + labs(title = TopMarkers[which(TopMarkers$gene == )]
    }

    print("Saving Plots")
    lapply(names(pTopDEgenes), FUN = function(x){
        ggsave(pTopDEgenes[[x]], file = paste0(outdir,"TopDEgenesCluster",x,".png"))
    })





}


