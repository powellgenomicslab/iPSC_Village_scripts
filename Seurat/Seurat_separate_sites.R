##### Read in libraries #####   
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(Seurat)
library(jcolors)
library(cowplot)
library(RColorBrewer)


##### Set up directories #####
dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
outdir <- paste0(dir,"output/Seurat_w_Hash/Seurat_separate_sites/")
datadir <- paste0(dir,"/data/PBMC_scRNA/")


##### Create 
if (!file.exists(paste0(dir,"output/Seurat_w_Hash/Seurat_separate_sites/seurat_norm_NoOutliers_by_site.RDS"))){
    seurat <- readRDS(paste0(dir,"output/Seurat_w_Hash/QC/seurat_norm_NoOutliers.rds"))
    seurat@meta.data$Site <- gsub("\\d","", seurat@meta.data$MULTI_ID)
    assigned <- read.table(paste0(dir,"output/Seurat_w_Hash/CellCycle/CellCycleProportions.txt"), sep = "\t")
    seurat <- AddMetaData(seurat, assigned)
    seurat$Time_FinalAssignment <- paste0(seurat$Time, seurat$FinalAssignment)

    seurat_sites <- list()
    seurat_sites[["Brisbane"]] <- subset(seurat, subset = Site == "Brisbane")
    seurat_sites[["Melbourne"]] <- subset(seurat, subset = Site == "Melbourne")
    seurat_sites[["Sydney"]] <- subset(seurat, subset = Site == "Sydney")

    saveRDS(seurat_sites, paste0(dir,"output/Seurat_w_Hash/Seurat_separate_sites/seurat_norm_NoOutliers_by_site.RDS"))

} else {
    readRDS(paste0(dir,"output/Seurat_w_Hash/Seurat_separate_sites/seurat_norm_NoOutliers_by_site.RDS"))
}


if ( stage == "SCTnormalization_MtRbRegressed"){

    if (!file.exists(paste0(outdir, "seurat_SCT_MtRbRegressed_separate_sites.rds")))
        seurat_sites <- lapply(seurat_sites, function(x){
            SCTransform(x, vars.to.regress = c("percent.mt","percent.rb"), verbose = TRUE, return.only.var.genes = FALSE)
        })

        seurat_sites <- lapply(seurat_sites, function(x){
            RunPCA(x, npcs = 100)
        })

        seurat_sites <- lapply(seurat_sites, function(x){
            RunUMAP(x, dims = 1:100, verbose = TRUE)
        })

        seurat_sites <- lapply(seurat_sites, function(x){
            FindNeighbors(x, dims = 1:100, verbose = TRUE)
        })

        saveRDS(seurat_sites, paste0(outdir, "seurat_SCT_MtRbRegressed_separate_sites.rds"))
    else {
        seurat_sites <- readRDS(paste0(outdir, "seurat_SCT_MtRbRegressed_separate_sites.rds"))
    }
    

    ##### Make Plots #####
    plots <- lapply(seurat_sites, function(x){
        DimPlot(x, group.by = c("Time", "FinalAssignment"), combine = FALSE)
    }

    plots <- lapply(plots, function(x){
        lapply(x, FUN = function(y) {
            y + theme(legend.position = "top") + guides(color = guide_legend(nrow = 4, 
                byrow = TRUE, override.aes = list(size = 2.5)))
        })
    })

    for (site in names(plots)){
        ggsave(CombinePlots(plots[[site]], ncol = 3), filename = paste0(outdir,site,"_umap_TimeIndividual.png"),width = 20, height = 9)
    }


    umaps_split <- list()

    umaps <- lapply(names(seurat_sites), function(x){
        umaps[["pool"]] <- DimPlot(seurat_sites[[x]], reduction = "umap", group.by = "Pool", pt.size = .01, split.by = 'Pool')
        umaps[["Time"]] <- DimPlot(seurat_sites[[x]], reduction = "umap", group.by = "Time", pt.size = .01, split.by = 'Time')
        umaps[["SampleReplicate"]] <- DimPlot(seurat_sites[[x]], reduction = "umap", group.by = "MULTI_ID", pt.size = .01, split.by = 'MULTI_ID')
        umaps[["CellLine"]] <- DimPlot(seurat_sites[[x]], reduction = "umap", group.by = "FinalAssignment", pt.size = .01, split.by = 'FinalAssignment')
        umaps[["SampleReplicate_Time"]] <- DimPlot(seurat_sites[[x]], reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'MULTI_ID')
        umaps[["CellLine_Time"]] <- DimPlot(seurat_sites[[x]], reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'FinalAssignment')
        umaps[["CellLine_Time_Replicate"]] <- DimPlot(seurat_sites[[x]], reduction = "umap", group.by = "Time_FinalAssignment", pt.size = .01, split.by = 'MULTI_ID')
        umaps[["mt_percent"]] <- FeaturePlot(seurat_sites[[x]], reduction = "umap", feature = "percent.mt")
        umaps[["rb_percent"]] <- FeaturePlot(seurat_sites[[x]], reduction = "umap", feature = "percent.rb")
        umaps[["nUMI"]] <- FeaturePlot(seurat_sites[[x]], reduction = "umap", feature = "nCount_RNA")
        umaps[["nFeatures"]] <- FeaturePlot(seurat_sites[[x]], reduction = "umap", feature = "nFeature_RNA")
        umaps[["CellCycle"]] <- DimPlot(seurat, group.by = c("phases"))  

        return(umaps)
    })

    for (name in names(umaps)){
        for (site in names(umaps[[name]])){
            ggsave(umaps[[name]][[site]], filename = paste0(outdir, "umap_",site,"_",name,".png"), height = 4, width = 10)
        }
    }
    ggsave(umap, filename = paste0(outdir, "umap_pool.png"), height = 4, width = 10)
    ggsave(umap, filename = paste0(outdir, "umap_Time.png"), height = 4, width = 10)
    ggsave(umap, filename = paste0(outdir, "umap_SampleReplicate.png"), height = 4, width = 10)
    ggsave(umap, filename = paste0(outdir, "umap_CellLine.png"), height = 4, width = 10)
    ggsave(umap, filename = paste0(outdir, "umap_SampleReplicate_Time.png"), height = 4, width = 10)
    ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time.png"), height = 4, width = 10)
    ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time_Site.png"), height = 4, width = 10)
    ggsave(UMAPmt, filename = paste0(outdir,"UMAPmt.png"))
    ggsave(UMAPrb, filename = paste0(outdir,"UMAPrb.png"))
    ggsave(umap, filename = paste0(outdir, "umap_CellCycle.png"))



