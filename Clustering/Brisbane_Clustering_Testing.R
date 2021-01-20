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


dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
outbase <- paste0(dir,"output/Brisbane_Clustering_Testing/")
dir.create(outbase)


stage = "SCTnormalization_no_covatiates"
outdir <- paste0(outbase,stage,"/")
dir.create(outdir)


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

    umap <- DimPlot(seurat,  pt.size = .01, group.by = c("Phase"))  
    ggsave(umap, filename = paste0(outdir, "umap_CellCycle.png"))



stage = "SCTnormalization_individual"
outdir <- paste0(outbase,stage,"/")
dir.create(outdir)


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
}


