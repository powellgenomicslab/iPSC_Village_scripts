library(dplyr, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/harmony/lib/R/library")
library(tidyr, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/harmony/lib/R/library")
library(tidyverse, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/harmony/lib/R/library")
library(ggplot2, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/harmony/lib/R/library")
library(clustree, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/harmony/lib/R/library")
library(reticulate, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/harmony/lib/R/library")
library(Seurat, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/harmony/lib/R/library")
library(schex, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/harmony/lib/R/library")
library(BiocParallel, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/harmony/lib/R/library")
library(coop, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/harmony/lib/R/library")
library(harmony, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/harmony/lib/R/library")
library(cowplot, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/harmony/lib/R/library")
library(RColorBrewer, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/harmony/lib/R/library")
library(readr, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/harmony/lib/R/library")
library(purrr, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/harmony/lib/R/library")

print(sessionInfo())


mad_function <- function(seurat, column, number_mad){
    mad <- mad(seurat@meta.data[,column])
    low <- median(seurat@meta.data[,column]) - number_mad*mad
    high <- median(seurat@meta.data[,column]) + number_mad*mad
    print("The mad is:")
    print(mad)
    print("The lower bound is:")
    print(low)
    print("The upper bound is:")
    print(high)
    seurat@meta.data[,paste0(column,"_mad")] <- ifelse((seurat@meta.data[,column] > low & seurat@meta.data[,column] < high),"NotOutlier", "Outlier")
    return(seurat)
}

PoolColors <- brewer.pal(n = 6, name = "Dark2")

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


if ( stage == "CheckHashtags"){
    pools <- dir("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/", pattern = "DRENEA")
    dirs_hash <- paste0("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Barcodes_200213_A00152_0206_AHK3HYDRXX/Hashing/", pools, "/umi_count/")
    dirs10x <- paste0("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/", pools, "/outs/filtered_feature_bc_matrix/")
    seurat <- readRDS("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/QC/AllCellsSeurat.rds")

    ### Read in the hashing data ###
    hash_list <- lapply(dirs_hash, function(x){
        Read10X(x, gene.column =1)
    })
    names(hash_list) <- pools

    hash_list <- lapply(hash_list, function(x){
        x <- as.data.frame(t(as.data.frame(x)))
        x$Barcode <- rownames(x)
        return(x)
    })
    names(hash_list) <- pools

    hash_list <- lapply(names(hash_list), function(x){
        left_join(barcodes[[x]], hash_list[[x]], by = c("Barcode"))
    })
    names(hash_list) <- pools

    hash_list <- lapply(names(hash_list), function(x){
        tmp <- as.data.frame(t(hash_list[[x]]))
        colnames(tmp) <- barcodes[[x]]$Barcode
        tmp <- tmp[2:5,]
    })
    names(hash_list) <- pools

    ### Add pool names to barcodes so all are unique ###
    hash_list <- lapply(names(hash_list), function(x){
        colnames(hash_list[[x]]) <- paste0(x, "_", colnames(hash_list[[x]]))
        return(hash_list[[x]])
    })
    names(hash_list) <- pools


    hash_list <- lapply(hash_list, function(x){
        x$Hashtag <- rownames(x)
        return(x)
    })
    names(hash_list) <- pools

    ### Combine hashing tables and correct the direction
    hash <- hash_list %>% reduce(full_join, by = "Hashtag")
    hash_barcode <- colnames(hash)
    hash_t_colnames <- hash$Hashtag
    hash_t <- as.data.frame(t(hash))
    colnames(hash_t) <- hash_t_colnames
    hash_t <- hash_t[-c(which(rownames(hash_t) == "Hashtag")),]
    hash_t <- as.data.frame(data.matrix(hash_t))
    hash_t[is.na(hash_t)] <- 0
    hash_t$Barcode <- colnames(hash)[-c(which(colnames(hash) == "Hashtag"))]

    write.table(hash_t, paste0(outdir,"hashtag_dataframe.txt"), sep = "\t")
    hash_t <- read.table(paste0(outdir,"hashtag_dataframe.txt"), sep = "\t")

    barcodes <- as.data.frame(colnames(seurat))
    colnames(barcodes) <- "Barcode"
    hash_t <- left_join(barcodes, hash_t, by = "Barcode")
    rownames(hash_t) <- hash_t$Barcode
    hash_t$Barcode <- NULL

    hash_t_t <- as.data.frame(t(hash_t))

    seurat[["HTO"]] <- CreateAssayObject(counts = hash_t_t)
    seurat <- NormalizeData(seurat, assay = "HTO", normalization.method = "CLR")
    table(seurat$HTO_classification.global)




} else if (stage[1]=="QC"){
    pools <- dir("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/", pattern = "DRENEA")

    ## Read in data
    counts_list <- lapply(dirs10x, function(x){
        Read10X(x)
    })
    names(counts_list) <- pools

    ## Add poolnames to cell names so can easily match back if there are any cells with same barcoe
    counts_list <- lapply(names(counts_list), function(x){
        colnames(counts_list[[x]]) <- paste0(x, "_", colnames(counts_list[[x]]))
        return(counts_list[[x]])
    })
    names(counts_list) <- pools



    ## Combine the matrices together
    counts <- do.call(cbind, counts_list)


    ## Make seurat object
    seurat <- CreateSeuratObject(counts = counts, meta.data = meta)

    ## Add QC metrics
    RbGeneList <- read.delim(file = "/directflow/SCCGGroupShare/projects/DrewNeavin/References/RibosomalGeneList.txt",header = F)
    print("Calculating Mt %")
    seurat <- PercentageFeatureSet(seurat, pattern = "^MT-", col.name = "percent.mt")
    print("Calculating Rb %")
    seurat <- PercentageFeatureSet(seurat, features = RbGeneList$V1, col.name = "percent.rb")

    write_delim(seurat@meta.data, paste0(outdir,"AllCells_metadata.txt"), delim = "\t")
    saveRDS(seurat, paste0(outdir,"AllCellsSeurat.rds"))
    seurat <- readRDS(paste0(outdir,"AllCellsSeurat.rds"))


    Barcodes <- as.data.frame(colnames(seurat))
    colnames(Barcodes) <- "Barcode"

    ## Get the metadata
    demultiplex_meta <- read_delim(paste0(dir,"output/CompareDemultiplexing/doublet_metadata.txt"), delim = "\t")

    ## Join with the barcode list ##
    meta <- left_join(Barcodes, demultiplex_meta, by = c("Barcode"))
    rownames(meta) <- meta$Barcode
    print(dim(meta))

    seurat@meta.data <- meta
    saveRDS(seurat, paste0(outdir,"AllCellsSeurat_meta.rds"))
    seurat <- readRDS(paste0(outdir,"AllCellsSeurat_meta.rds"))

    ## Add the hashing data ##
    seurat[["HTO"]] <- CreateAssayObject(counts = hash)
    seurat <- NormalizeData(seurat, assay = "HTO", normalization.method = "CLR")

    ## Scrublet vs demultiplexing tables ##
    print(prop.table(table(seurat@meta.data$FinalAssignment, gsub(TRUE,"scrublet_doublet",seurat@meta.data$ScrubletDoublet) %>% gsub(FALSE,"scrublet_singlet",.)), margin = 1))
    print(prop.table(table(seurat@meta.data$FinalAssignment, gsub(TRUE,"scrublet_doublet",seurat@meta.data$ScrubletDoublet) %>% gsub(FALSE,"scrublet_singlet",.)), margin = 2))

    print("The seurat object before removing doublets:")
    print(seurat)
    Idents(seurat) <- "ScrubletDoublet"
    seurat <- subset(seurat, idents = FALSE)
    Idents(seurat) <- "FinalAssignment"
    seurat <- subset(seurat,  idents = c("FSA006", "TOB421","MBE1006"))
    print("The seurat object after removing doublets")
    print(seurat)

    ## Calculate MADs
    seurat <- mad_function(seurat, "percent.mt",3)
    seurat <- mad_function(seurat, "percent.rb",3)
    seurat <- mad_function(seurat, "nCount_RNA", 3)
    seurat <- mad_function(seurat, "nFeature_RNA", 3)


    plot_mt_pct <- VlnPlot(seurat, features = c( "percent.mt"), group.by = "Pool", pt.size = 0, cols = PoolColors) +
                    geom_hline(yintercept=0, linetype="dashed", color = "black") +
                    geom_hline(yintercept=13.732771, linetype="dashed", color = "black")
    ggsave(plot_mt_pct, filename = paste0(outdir,"Mt_pct_vln.png"))

    plot_rb_pct <- VlnPlot(seurat, features = c( "percent.rb"), group.by = "Pool", pt.size = 0, cols = PoolColors) +
                    geom_hline(yintercept=15.42803, linetype="dashed", color = "black") +
                    geom_hline(yintercept=37.07026, linetype="dashed", color = "black")
    ggsave(plot_rb_pct, filename = paste0(outdir,"Rb_pct_vln.png"))

    plot_n_count <- VlnPlot(seurat, features = c( "nCount_RNA"), group.by = "Pool", pt.size = 0, cols = PoolColors) +
                    geom_hline(yintercept=0, linetype="dashed", color = "black") +
                    geom_hline(yintercept=44548.25, linetype="dashed", color = "black")
    ggsave(plot_n_count, filename = paste0(outdir,"N_count_vln.png"))

    plot_nFeature_RNA <- VlnPlot(seurat, features = c( "nFeature_RNA"), group.by = "Pool", pt.size = 0, cols = PoolColors) +
                    geom_hline(yintercept=1192.719, linetype="dashed", color = "black") +
                    geom_hline(yintercept=7535.281, linetype="dashed", color = "black")
    ggsave(plot_nFeature_RNA, filename = paste0(outdir,"nFeature_RNA_vln.png"))

    lib_mt <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Pool", cols = alpha(brewer.pal(n = 6, name = "Dark2"), 0.25)) +
                    geom_vline(xintercept=0, linetype="dashed", color = "black") +
                    geom_vline(xintercept=44548.25, linetype="dashed", color = "black") +
                    geom_hline(yintercept=0, linetype="dashed", color = "black") +
                    geom_hline(yintercept=13.732771, linetype="dashed", color = "black") 
    ggsave(lib_mt, filename = paste0(outdir,"NoDoublets_lib_mt.png"))
    lib_genes <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Pool", cols = alpha(brewer.pal(n = 6, name = "Dark2"), 0.25)) +
                    geom_vline(xintercept=0, linetype="dashed", color = "black") +
                    geom_vline(xintercept=44548.25, linetype="dashed", color = "black") +
                    geom_hline(yintercept=1192.719, linetype="dashed", color = "black") +
                    geom_hline(yintercept=7535.281, linetype="dashed", color = "black")
    ggsave(lib_genes, filename = paste0(outdir,"NoDoublets_lib_genes.png"))

    ##### Remove the outliers #####
    seurat <- subset(seurat, subset = percent.mt_mad == "NotOutlier") 
    seurat <- subset(seurat, subset = percent.rb_mad == "NotOutlier")
    seurat <- subset(seurat, subset = nCount_RNA_mad == "NotOutlier")
    seurat <- subset(seurat, subset = nFeature_RNA_mad == "NotOutlier")
    saveRDS(seurat, paste0(outdir,"NoOutliers.rds"))


    plot_mt_pct <- VlnPlot(seurat, features = c( "percent.mt"), group.by = "Pool", pt.size = 0, cols = PoolColors)
    ggsave(plot_mt_pct, filename = paste0(outdir,"Mt_pct_vln_NoOutliers.png"))

    plot_rb_pct <- VlnPlot(seurat, features = c( "percent.rb"), group.by = "Pool", pt.size = 0, cols = PoolColors)
    ggsave(plot_rb_pct, filename = paste0(outdir,"Rb_pct_vln_NoOutliers.png"))

    plot_n_count <- VlnPlot(seurat, features = c( "nCount_RNA"), group.by = "Pool", pt.size = 0, cols = PoolColors)
    ggsave(plot_n_count, filename = paste0(outdir,"N_count_vln_NoOutliers.png"))

    plot_nFeature_RNA <- VlnPlot(seurat, features = c( "nFeature_RNA"), group.by = "Pool", pt.size = 0, cols = PoolColors)
    ggsave(plot_nFeature_RNA, filename = paste0(outdir,"nFeature_RNA_vln_NoOutliers.png"))

    lib_mt <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Pool", cols = alpha(brewer.pal(n = 6, name = "Dark2"), 0.25))
    ggsave(lib_mt, filename = paste0(outdir,"NoDoublets_lib_mt_NoOutliers.png"))
    lib_genes <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Pool", cols = alpha(brewer.pal(n = 6, name = "Dark2"), 0.25))
    ggsave(lib_genes, filename = paste0(outdir,"NoDoublets_lib_genes_NoOutliers.png"))




} else if ( stage == "SCTnormalization_no_covatiates"){
    seurat <- readRDS(paste0(dir,"output/Seurat/QC/NoOutliers.rds"))
    seurat <- SCTransform(seurat, verbose = TRUE)
    seurat <- RunPCA(seurat, npcs = 100)
    seurat <- RunUMAP(seurat, dims = 1:100, verbose = TRUE)
    seurat <- FindNeighbors(seurat, dims = 1:100, verbose = TRUE)
    saveRDS(seurat, paste0(outdir, "seurat_SCT_noCov.rds"))

} else if ( stage == "SCTnormalization_PoolRegressed" ){
    seurat <- readRDS(paste0(dir,"output/Seurat/QC/NoOutliers.rds"))
    seurat <- SCTransform(seurat, vars.to.regress = "Pool", verbose = TRUE)
    seurat <- RunPCA(seurat, npcs = 100)
    seurat <- RunUMAP(seurat, dims = 1:100, verbose = TRUE)
    seurat <- FindNeighbors(seurat, dims = 1:100, verbose = TRUE)
    saveRDS(seurat, paste0(outdir, "seurat_SCT_PoolRegressed.rds"))

} else if ( stage == "SCTnormalization_PoolDayRegressed"){
    seurat <- readRDS(paste0(dir,"output/Seurat/QC/NoOutliers.rds"))
    seurat@meta.data$Day <- NA
    seurat@meta.data$Day <- ifelse((seurat@meta.data$Pool == "DRENEA_1" | seurat@meta.data$Pool == "DRENEA_2" | seurat@meta.data$Pool == "DRENEA_3"), "Day0", "Day4")
    print(head(seurat@meta.data$Day))
    print(tail(seurat@meta.data$Day))
    seurat <- SCTransform(seurat, vars.to.regress = c("Pool", "Day"), verbose = TRUE)
    seurat <- RunPCA(seurat, npcs = 100)
    seurat <- RunUMAP(seurat, dims = 1:100, verbose = TRUE)
    seurat <- FindNeighbors(seurat, dims = 1:100, verbose = TRUE)
    saveRDS(seurat, paste0(outdir, "seurat_SCT_PoolDayRegressed.rds"))

} else if (stage[1]=="Harmony_Pool"){
    seurat <- readRDS(paste0(dir,"output/Seurat/QC/NoOutliers.rds"))

    seurat <- NormalizeData(seurat, verbose = FALSE) %>%
        FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
        ScaleData(verbose = FALSE) %>% 
        RunPCA(pc.genes = seurat@var.genes, npcs = 20, verbose = FALSE)
    seurat <- RunHarmony(seurat,"Pool")

    ### Plot the metrics for the results ###
    p1 <- DimPlot(object = seurat, reduction = "harmony", pt.size = .1, group.by = "Pool", do.return = TRUE)
    ggsave(p1, filename = paste0(outdir,"harmony_PCA.png"))

    p2 <- VlnPlot(object = seurat, features = "harmony_1", group.by = "Pool", pt.size = .1)
    ggsave(p2, filename = paste0(outdir,"harmony_violin.png"))

    seurat <- RunUMAP(seurat, reduction = "harmony", dims = 1:20)
    seurat <- FindNeighbors(seurat, reduction = "harmony", dims = 1:20)

    umap <- DimPlot(seurat, reduction = "umap", group.by = "Pool", pt.size = .01, split.by = 'Pool')
    ggsave(umap, filename = paste0(outdir, "umap.png"), height = 4, width = 10)

    saveRDS(seurat, paste0(outdir,"NoOutliers_Pool_harmony.rds"))

} else if (stage[1]=="Harmony_Day"){
    seurat <- readRDS(paste0(dir,"output/Seurat/QC/NoOutliers.rds"))
    seurat@meta.data$Day <- NA
    seurat@meta.data$Day <- ifelse((seurat@meta.data$Pool == "DRENEA_1" | seurat@meta.data$Pool == "DRENEA_2" | seurat@meta.data$Pool == "DRENEA_3"), "Day0", "Day4")
    print(head(seurat@meta.data$Day))
    print(tail(seurat@meta.data$Day))

    seurat <- NormalizeData(seurat, verbose = FALSE) %>%
        FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
        ScaleData(verbose = FALSE) %>% 
        RunPCA(pc.genes = seurat@var.genes, npcs = 20, verbose = FALSE)
    seurat <- RunHarmony(seurat,"Day", plot_convergence = TRUE)

    ### Plot the metrics for the results ###
    p1 <- DimPlot(object = seurat, reduction = "harmony", pt.size = .1, group.by = "Day")
    ggsave(p1, filename = paste0(outdir,"harmony_PCA.png"))

    p2 <- VlnPlot(object = seurat, features = "harmony_1", group.by = "Day", pt.size = .1)
    ggsave(p2, filename = paste0(outdir,"harmony_violin.png"))


    seurat <- RunUMAP(seurat, reduction = "harmony", dims = 1:20)
    seurat <- FindNeighbors(seurat, reduction = "harmony", dims = 1:20)

    umap <- DimPlot(seurat, reduction = "umap", group.by = "Pool", pt.size = .01, split.by = 'Day')
    ggsave(umap, filename = paste0(outdir, "umap.png"), height = 4, width = 10)

    saveRDS(seurat, paste0(outdir,"NoOutliers_Day_harmony.rds"))

} else if (stage[1]=="Harmony_DayPatient"){
    seurat <- readRDS(paste0(dir,"output/Seurat/QC/NoOutliers.rds"))
    seurat@meta.data$Day <- NA
    seurat@meta.data$Day <- ifelse((seurat@meta.data$Pool == "DRENEA_1" | seurat@meta.data$Pool == "DRENEA_2" | seurat@meta.data$Pool == "DRENEA_3"), "Day0", "Day4")
    print(head(seurat@meta.data$Day))
    print(tail(seurat@meta.data$Day))

    seurat <- SCTransform(seurat, vars.to.regress = "Pool", verbose = TRUE)
    seurat <- RunPCA(seurat, npcs = 100)
    seurat <- RunHarmony(seurat,c("Day", "FinalAssignment"), plot_convergence = TRUE, assay.use = "SCT")

    ### Plot the metrics for the results ###
    p1 <- DimPlot(object = seurat, reduction = "harmony", pt.size = .1, group.by = "Day")
    ggsave(p1, filename = paste0(outdir,"harmony_PCA.png"))

    p2 <- VlnPlot(object = seurat, features = "harmony_1", group.by = "Day", pt.size = .1)
    ggsave(p2, filename = paste0(outdir,"harmony_violin.png"))

    seurat <- RunUMAP(seurat, reduction = "harmony", dims = 1:20)
    seurat <- FindNeighbors(seurat, reduction = "harmony", dims = 1:20)

    umap <- DimPlot(seurat, reduction = "umap", group.by = "Pool", pt.size = .01, split.by = 'Day')
    ggsave(umap, filename = paste0(outdir, "umap.png"), height = 4, width = 10)

    saveRDS(seurat, paste0(outdir,"NoOutliers_DayPatient_harmony.rds"))


} else if (stage[1]=="Integration_DayPatient"){
    # seurat <- readRDS(paste0(dir,"output/Seurat/QC/NoOutliers.rds"))
    # seurat@meta.data$Day <- NA
    # seurat@meta.data$Day <- ifelse((seurat@meta.data$Pool == "DRENEA_1" | seurat@meta.data$Pool == "DRENEA_2" | seurat@meta.data$Pool == "DRENEA_3"), "Day0", "Day4")
    # seurat@meta.data$Day_Patient <- paste0(seurat@meta.data$Day, "_", seurat@meta.data$FinalAssignment)
    # print(unique( seurat@meta.data$Day_Patient))

    # seurat_list <- SplitObject(seurat, split.by = "Day_Patient")

    # seurat_list <- lapply(seurat_list, function(x) {
    #     SCTransform(x, verbose = TRUE)
    # })

    # seurat_features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)
    # seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = seurat_features)
    # seurat_anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT", 
    #     anchor.features = seurat_features)
    # seurat_integrated <- IntegrateData(anchorset = seurat_anchors, normalization.method = "SCT")

    # seurat_integrated <- RunPCA(object = seurat_integrated, verbose = FALSE)
    # seurat_integrated <- RunUMAP(object = seurat_integrated, dims = 1:30)

    seurat_integrated <- readRDS(paste0(outdir,"NoOutliers_DayPatient_SeuratIntegratee.rds"))

    seurat_integrated <- FindNeighbors(seurat_integrated, reduction = "pca", assay = "integrated", dims = 1:20)


    # plots <- DimPlot(seurat_integrated, group.by = c("Day", "FinalAssignment"), combine = FALSE)
    # plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 4, 
    #     byrow = TRUE, override.aes = list(size = 2.5))))
    # ggsave(CombinePlots(plots), filename = paste0(outdir,"umap.png"),width = 12, height = 9)

    saveRDS(seurat_integrated, paste0(outdir,"NoOutliers_DayPatient_SeuratIntegratee.rds"))


} else if (stage[1]=="Integration_Day"){
    seurat <- readRDS(paste0(dir,"output/Seurat/QC/NoOutliers.rds"))
    seurat@meta.data$Day <- NA
    seurat@meta.data$Day <- ifelse((seurat@meta.data$Pool == "DRENEA_1" | seurat@meta.data$Pool == "DRENEA_2" | seurat@meta.data$Pool == "DRENEA_3"), "Day0", "Day4")

    seurat_list <- SplitObject(seurat, split.by = "Day")

    seurat_list <- lapply(seurat_list, function(x) {
        SCTransform(x, verbose = TRUE)
    })

    seurat_features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)
    seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = seurat_features)
    seurat_anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT", 
        anchor.features = seurat_features)
    seurat_integrated <- IntegrateData(anchorset = seurat_anchors, normalization.method = "SCT")

    seurat_integrated <- RunPCA(object = seurat_integrated, verbose = FALSE)
    seurat_integrated <- RunUMAP(object = seurat_integrated, dims = 1:30)

    seurat_integrated <- FindNeighbors(seurat_integrated, reduction = "pca", assay = "integrated", dims = 1:20)


    plots <- DimPlot(seurat_integrated, group.by = c("Day", "FinalAssignment"), combine = FALSE)
    plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 4, 
        byrow = TRUE, override.aes = list(size = 2.5))))
    ggsave(CombinePlots(plots), filename = paste0(outdir,"umap.png"),width = 12, height = 9)

    saveRDS(seurat_integrated, paste0(outdir,"NoOutliers_Day_SeuratIntegratee.rds"))


} else if (stage[1]=="Integration_Day_then_Individual"){
    seurat <- readRDS(paste0(dir,"output/Seurat/QC/NoOutliers.rds"))
    seurat@meta.data$Day <- NA
    seurat@meta.data$Day <- ifelse((seurat@meta.data$Pool == "DRENEA_1" | seurat@meta.data$Pool == "DRENEA_2" | seurat@meta.data$Pool == "DRENEA_3"), "Day0", "Day4")

    seurat_list <- SplitObject(seurat, split.by = "Day")

    seurat_list <- lapply(seurat_list, function(x) {
        SCTransform(x, verbose = TRUE)
    })

    seurat_features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)
    seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = seurat_features)
    seurat_anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT", 
        anchor.features = seurat_features)
    seurat_integrated <- IntegrateData(anchorset = seurat_anchors, normalization.method = "SCT")

    ### Second round of integration ###
    seurat_integrated_list <- SplitObject(seurat_integrated, split.by = "FinalAssignment")
    seurat_integrated_list <- lapply(seurat_integrated_list, function(x) {
        SCTransform(x, verbose = TRUE)
    })

    seurat_features <- SelectIntegrationFeatures(object.list = seurat_integrated_list, nfeatures = 3000, assay = c("integrated","integrated","integrated"))
    seurat_integrated_list <- PrepSCTIntegration(object.list = seurat_integrated_list, anchor.features = seurat_features, assay = "integrated")
    seurat_anchors <- FindIntegrationAnchors(object.list = seurat_integrated_list, normalization.method = "SCT", 
        anchor.features = seurat_features, assay = c("integrated","integrated","integrated"))
    seurat_integrated_2 <- IntegrateData(anchorset = seurat_anchors, normalization.method = "SCT", new.assay.name = "integrated_2")


    seurat_integrated_2 <- RunPCA(object = seurat_integrated_2, verbose = FALSE, assay = "integrated_2")
    seurat_integrated_2 <- RunUMAP(object = seurat_integrated_2, dims = 1:30, assay = "integrated_2")


    seurat_integrated_2 <- FindNeighbors(seurat_integrated_2, reduction = "pca", assay = "integrated_2", dims = 1:20)


    plots <- DimPlot(seurat_integrated_2, group.by = c("Day", "FinalAssignment"), combine = FALSE)
    plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 4, 
        byrow = TRUE, override.aes = list(size = 2.5))))
    ggsave(CombinePlots(plots), filename = paste0(outdir,"umap.png"),width = 12, height = 9)

    saveRDS(seurat_integrated_2, paste0(outdir,"NoOutliers_Day_SeuratIntegratee.rds"))



} else if (stage[1]=="Multi_Resolution_no_covatiates"){
    seurat <- readRDS(paste0(dir,"output/Seurat/SCTnormalization_no_covatiates/seurat_SCT_noCov.rds"))

    ### Find the clusters for this resolution ###
    res <- sge*0.1-0.1
    print(paste0("Finding Clusters for resolution=",res))
    seurat <- FindClusters(seurat, resolution = res)

    ### Plot the UMAP for this resolution ###
    UMAP <- DimPlot(seurat, reduction = "umap")
    ggsave(file = paste0(outdir,"UMAPresolution",res,".png"), plot = UMAP)

    ### Save the identities for this resolution ###
    identities <- as.data.frame(Idents(seurat))
    colnames(identities) <- c(paste0("Resolution_",res))

    write_delim(identities, paste0(outdir,"Identities_Resolution_",res,".txt"), delim = "\t")

} else if (stage[1]=="Multi_Resolution_PoolRegressed"){
    seurat <- readRDS(paste0(dir,"output/Seurat/SCTnormalization_PoolRegressed/seurat_SCT_PoolRegressed.rds"))

    ### Find the clusters for this resolution ###
    res <- sge*0.1-0.1
    print(paste0("Finding Clusters for resolution=",res))
    seurat <- FindClusters(seurat, resolution = res)

    ### Plot the UMAP for this resolution ###
    UMAP <- DimPlot(seurat, reduction = "umap")
    ggsave(file = paste0(outdir,"UMAPresolution",res,".png"), plot = UMAP)

    ### Save the identities for this resolution ###
    identities <- as.data.frame(Idents(seurat))
    colnames(identities) <- c(paste0("Resolution_",res))

    write_delim(identities, paste0(outdir,"Identities_Resolution_",res,".txt"), delim = "\t")


} else if (stage[1]=="Multi_Resolution_PoolDayRegressed"){
    seurat <- readRDS(paste0(dir,"output/Seurat/SCTnormalization_PoolDayRegressed/seurat_SCT_PoolDayRegressed.rds"))

    ### Find the clusters for this resolution ###
    res <- sge*0.1-0.1
    print(paste0("Finding Clusters for resolution=",res))
    seurat <- FindClusters(seurat, resolution = res)

    ### Plot the UMAP for this resolution ###
    UMAP <- DimPlot(seurat, reduction = "umap")
    ggsave(file = paste0(outdir,"UMAPresolution",res,".png"), plot = UMAP)

    ### Save the identities for this resolution ###
    identities <- as.data.frame(Idents(seurat))
    colnames(identities) <- c(paste0("Resolution_",res))

    write_delim(identities, paste0(outdir,"Identities_Resolution_",res,".txt"), delim = "\t")


} else if (stage[1]=="Multi_Resolution_Harmony_Pool"){
    seurat <- readRDS(paste0(dir,"output/Seurat/Harmony_Pool/NoOutliers_Pool_harmony.rds"))

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


} else if (stage[1]=="Multi_Resolution_Harmony_Day"){
    seurat <- readRDS(paste0(dir,"output/Seurat/Harmony_Day/NoOutliers_Day_harmony.rds"))

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



} else if (stage[1]=="Multi_Resolution_Harmony_DayPatient"){
    seurat <- readRDS(paste0(dir,"output/Seurat/Harmony_DayPatient/NoOutliers_DayPatient_harmony.rds"))

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


} else if (stage[1]=="Multi_Resolution_Integration_Day_then_Individual"){
    seurat <- readRDS(paste0(dir,"output/Seurat/Integration_Day_then_Individual/NoOutliers_Day_SeuratIntegratee.rds"))

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



} else if (stage[1]=="Clustree_no_covatiates"){
    ResolutionsFileList <- list.files(paste0(dir,"output/Seurat/Multi_Resolution_no_covatiates/"), pattern = "Identities_Resolution_")

    ResolutionsList <- lapply(ResolutionsFileList, FUN = function(x){
        read_delim(paste0(dir,"output/Seurat/Multi_Resolution_no_covatiates/",x), delim = "\t")
    })
    
    ResolutionsDF <- do.call(cbind, ResolutionsList)

    ##### CLUSTREE #####
    print("Using Clustree to diagram changes with increasing resolutions")
    pClustree <- clustree(ResolutionsDF, prefix = "Resolution_")
    ggsave(file = paste0(outdir, "ClusteringClustree.png"), plot = pClustree)


} else if (stage[1]=="Clustree_PoolRegressed"){
    ResolutionsFileList <- list.files(paste0(dir,"output/Seurat/Multi_Resolution_PoolRegressed/"), pattern = "Identities_Resolution_")

    ResolutionsList <- lapply(ResolutionsFileList, FUN = function(x){
        read_delim(paste0(dir,"output/Seurat/Multi_Resolution_PoolRegressed/",x), delim = "\t")
    })
    
    ResolutionsDF <- do.call(cbind, ResolutionsList)

    ##### CLUSTREE #####
    print("Using Clustree to diagram changes with increasing resolutions")
    pClustree <- clustree(ResolutionsDF, prefix = "Resolution_")
    ggsave(file = paste0(outdir, "ClusteringClustree.png"), plot = pClustree)


} else if (stage[1]=="Clustree_PoolDayRegressed"){
    ResolutionsFileList <- list.files(paste0(dir,"output/Seurat/Multi_Resolution_PoolDayRegressed/"), pattern = "Identities_Resolution_")

    ResolutionsList <- lapply(ResolutionsFileList, FUN = function(x){
        read_delim(paste0(dir,"output/Seurat/Multi_Resolution_PoolDayRegressed/",x), delim = "\t")
    })
    
    ResolutionsDF <- do.call(cbind, ResolutionsList)

    ##### CLUSTREE #####
    print("Using Clustree to diagram changes with increasing resolutions")
    pClustree <- clustree(ResolutionsDF, prefix = "Resolution_")
    ggsave(file = paste0(outdir, "ClusteringClustree.png"), plot = pClustree)


} else if (stage[1]=="Clustree_Harmony_Pool"){
    ResolutionsFileList <- list.files(paste0(dir,"output/Seurat/Multi_Resolution_Harmony_Pool/"), pattern = "Identities_Resolution_")

    ResolutionsList <- lapply(ResolutionsFileList, FUN = function(x){
        read_delim(paste0(dir,"output/Seurat/Multi_Resolution_Harmony_Pool/",x), delim = "\t")
    })
    
    ResolutionsDF <- do.call(cbind, ResolutionsList)

    ##### CLUSTREE #####
    print("Using Clustree to diagram changes with increasing resolutions")
    pClustree <- clustree(ResolutionsDF, prefix = "Resolution_")
    ggsave(file = paste0(outdir, "ClusteringClustree.png"), plot = pClustree)


} else if (stage[1]=="Clustree_Harmony_Day"){
    ResolutionsFileList <- list.files(paste0(dir,"output/Seurat/Multi_Resolution_Harmony_Day/"), pattern = "Identities_Resolution_")

    ResolutionsList <- lapply(ResolutionsFileList, FUN = function(x){
        read_delim(paste0(dir,"output/Seurat/Multi_Resolution_Harmony_Day/",x), delim = "\t")
    })
    
    ResolutionsDF <- do.call(cbind, ResolutionsList)

    ##### CLUSTREE #####
    print("Using Clustree to diagram changes with increasing resolutions")
    pClustree <- clustree(ResolutionsDF, prefix = "Resolution_")
    ggsave(file = paste0(outdir, "ClusteringClustree.png"), plot = pClustree)


} else if (stage[1]=="Clustree_Harmony_DayPatient"){
    ResolutionsFileList <- list.files(paste0(dir,"output/Seurat/Multi_Resolution_Harmony_DayPatient/"), pattern = "Identities_Resolution_")

    ResolutionsList <- lapply(ResolutionsFileList, FUN = function(x){
        read_delim(paste0(dir,"output/Seurat/Multi_Resolution_Harmony_DayPatient/",x), delim = "\t")
    })
    
    ResolutionsDF <- do.call(cbind, ResolutionsList)

    ##### CLUSTREE #####
    print("Using Clustree to diagram changes with increasing resolutions")
    pClustree <- clustree(ResolutionsDF, prefix = "Resolution_")
    ggsave(file = paste0(outdir, "ClusteringClustree.png"), plot = pClustree)



} else if (stage[1]=="Clustree_Integration_DayPatient"){
    ResolutionsFileList <- list.files(paste0(dir,"output/Seurat/Multi_Resolution_Integration_DayPatient/"), pattern = "Identities_Resolution_")

    ResolutionsList <- lapply(ResolutionsFileList, FUN = function(x){
        read_delim(paste0(dir,"output/Seurat/Multi_Resolution_Integration_DayPatient/",x), delim = "\t")
    })
    
    ResolutionsDF <- do.call(cbind, ResolutionsList)

    ##### CLUSTREE #####
    print("Using Clustree to diagram changes with increasing resolutions")
    pClustree <- clustree(ResolutionsDF, prefix = "Resolution_")
    ggsave(file = paste0(outdir, "ClusteringClustree.png"), plot = pClustree)


} else if (stage[1]=="Clustree_Integration_Day_then_Individual"){
    ResolutionsFileList <- list.files(paste0(dir,"output/Seurat/Multi_Resolution_Integration_Day_then_Individual/"), pattern = "Identities_Resolution_")

    ResolutionsList <- lapply(ResolutionsFileList, FUN = function(x){
        read_delim(paste0(dir,"output/Seurat/Multi_Resolution_Integration_Day_then_Individual/",x), delim = "\t")
    })
    
    ResolutionsDF <- do.call(cbind, ResolutionsList)

    ##### CLUSTREE #####
    print("Using Clustree to diagram changes with increasing resolutions")
    pClustree <- clustree(ResolutionsDF, prefix = "Resolution_")
    ggsave(file = paste0(outdir, "ClusteringClustree.png"), plot = pClustree)




} else if (stage[1]=="QCfigurew_no_covatiates"){
    seurat <- readRDS(paste0(dir,"output/Seurat/SCTnormalization_no_covatiates/seurat_SCT_noCov.rds"))
    UMAPpool <- DimPlot(seurat, reduction = "umap", group.by = "Pool")
    ggsave(UMAPpool, filename = paste0(outdir,"UMAPpool.png"))
    UMAPmt <- FeaturePlot(seurat, reduction = "umap", feature = "percent.mt")
    ggsave(UMAPmt, filename = paste0(outdir,"UMAPmt.png"))
    UMAPrb <- FeaturePlot(seurat, reduction = "umap", feature = "percent.rb")
    ggsave(UMAPrb, filename = paste0(outdir,"UMAPrb.png"))
    UMAPlib <- FeaturePlot(seurat, reduction = "umap", feature = "nCount_RNA")
    ggsave(UMAPlib, filename = paste0(outdir,"UMAPlib.png"))
    UMAPgenes <- FeaturePlot(seurat, reduction = "umap", feature = "nFeature_RNA")
    ggsave(UMAPgenes, filename = paste0(outdir,"UMAPgenes.png"))
    UMAPlibSCT <- FeaturePlot(seurat, reduction = "umap", feature = "nCount_SCT")
    ggsave(UMAPlibSCT, filename = paste0(outdir,"UMAPlibSCT.png"))
    UMAPgenesSCT <- FeaturePlot(seurat, reduction = "umap", feature = "nFeature_SCT")
    ggsave(UMAPgenesSCT, filename = paste0(outdir,"UMAPgenesSCT.png"))
    UMAPdonor_souporcell <- DimPlot(seurat, reduction = "umap", group.by = "common_assignment_souporcell")
    ggsave(UMAPdonor_souporcell, filename = paste0(outdir,"UMAPdonor_souporcell.png"))
    UMAPdonor_scSplit <- DimPlot(seurat, reduction = "umap", group.by = "common_assignment_scSplit")
    ggsave(UMAPdonor_scSplit, filename = paste0(outdir,"UMAPdonor_scSplit.png"))
    UMAPdonor_freemuxlet <- DimPlot(seurat, reduction = "umap", group.by = "common_assignment_freemuxlet")
    ggsave(UMAPdonor_freemuxlet, filename = paste0(outdir,"UMAPdonor_freemuxlet.png"))

} else if (stage[1]=="QCfigurew_PoolRegressed"){
    seurat <- readRDS(paste0(dir,"output/Seurat/SCTnormalization_PoolRegressed/seurat_SCT_PoolRegressed.rds"))
    UMAPpool <- DimPlot(seurat, reduction = "umap", group.by = "Pool")
    ggsave(UMAPpool, filename = paste0(outdir,"UMAPpool.png"))
    UMAPmt <- FeaturePlot(seurat, reduction = "umap", feature = "percent.mt")
    ggsave(UMAPmt, filename = paste0(outdir,"UMAPmt.png"))
    UMAPrb <- FeaturePlot(seurat, reduction = "umap", feature = "percent.rb")
    ggsave(UMAPrb, filename = paste0(outdir,"UMAPrb.png"))
    UMAPlib <- FeaturePlot(seurat, reduction = "umap", feature = "nCount_RNA")
    ggsave(UMAPlib, filename = paste0(outdir,"UMAPlib.png"))
    UMAPgenes <- FeaturePlot(seurat, reduction = "umap", feature = "nFeature_RNA")
    ggsave(UMAPgenes, filename = paste0(outdir,"UMAPgenes.png"))
    UMAPlibSCT <- FeaturePlot(seurat, reduction = "umap", feature = "nCount_SCT")
    ggsave(UMAPlibSCT, filename = paste0(outdir,"UMAPlibSCT.png"))
    UMAPgenesSCT <- FeaturePlot(seurat, reduction = "umap", feature = "nFeature_SCT")
    ggsave(UMAPgenesSCT, filename = paste0(outdir,"UMAPgenesSCT.png"))
    UMAPdonor_souporcell <- DimPlot(seurat, reduction = "umap", group.by = "common_assignment_souporcell")
    ggsave(UMAPdonor_souporcell, filename = paste0(outdir,"UMAPdonor_souporcell.png"))
    UMAPdonor_scSplit <- DimPlot(seurat, reduction = "umap", group.by = "common_assignment_scSplit")
    ggsave(UMAPdonor_scSplit, filename = paste0(outdir,"UMAPdonor_scSplit.png"))
    UMAPdonor_freemuxlet <- DimPlot(seurat, reduction = "umap", group.by = "common_assignment_freemuxlet")
    ggsave(UMAPdonor_freemuxlet, filename = paste0(outdir,"UMAPdonor_freemuxlet.png"))


} else if (stage[1]=="QCfigurew_PoolDayRegressed"){
    seurat <- readRDS(paste0(dir,"output/Seurat/SCTnormalization_PoolDayRegressed/seurat_SCT_PoolDayRegressed.rds"))
    UMAPpool <- DimPlot(seurat, reduction = "umap", group.by = "Pool")
    ggsave(UMAPpool, filename = paste0(outdir,"UMAPpool.png"))
    UMAPmt <- FeaturePlot(seurat, reduction = "umap", feature = "percent.mt")
    ggsave(UMAPmt, filename = paste0(outdir,"UMAPmt.png"))
    UMAPrb <- FeaturePlot(seurat, reduction = "umap", feature = "percent.rb")
    ggsave(UMAPrb, filename = paste0(outdir,"UMAPrb.png"))
    UMAPlib <- FeaturePlot(seurat, reduction = "umap", feature = "nCount_RNA")
    ggsave(UMAPlib, filename = paste0(outdir,"UMAPlib.png"))
    UMAPgenes <- FeaturePlot(seurat, reduction = "umap", feature = "nFeature_RNA")
    ggsave(UMAPgenes, filename = paste0(outdir,"UMAPgenes.png"))
    UMAPlibSCT <- FeaturePlot(seurat, reduction = "umap", feature = "nCount_SCT")
    ggsave(UMAPlibSCT, filename = paste0(outdir,"UMAPlibSCT.png"))
    UMAPgenesSCT <- FeaturePlot(seurat, reduction = "umap", feature = "nFeature_SCT")
    ggsave(UMAPgenesSCT, filename = paste0(outdir,"UMAPgenesSCT.png"))
    UMAPdonor_souporcell <- DimPlot(seurat, reduction = "umap", group.by = "common_assignment_souporcell")
    ggsave(UMAPdonor_souporcell, filename = paste0(outdir,"UMAPdonor_souporcell.png"))
    UMAPdonor_scSplit <- DimPlot(seurat, reduction = "umap", group.by = "common_assignment_scSplit")
    ggsave(UMAPdonor_scSplit, filename = paste0(outdir,"UMAPdonor_scSplit.png"))
    UMAPdonor_freemuxlet <- DimPlot(seurat, reduction = "umap", group.by = "common_assignment_freemuxlet")
    ggsave(UMAPdonor_freemuxlet, filename = paste0(outdir,"UMAPdonor_freemuxlet.png"))


} else if (stage[1]=="QCfigures_Harmony_Pool"){
    seurat <- readRDS(paste0(dir,"output/Seurat/Harmony_Pool/NoOutliers_Pool_harmony.rds"))
    UMAPpool <- DimPlot(seurat, reduction = "umap", group.by = "Pool")
    ggsave(UMAPpool, filename = paste0(outdir,"UMAPpool.png"))
    UMAPmt <- FeaturePlot(seurat, reduction = "umap", feature = "percent.mt")
    ggsave(UMAPmt, filename = paste0(outdir,"UMAPmt.png"))
    UMAPrb <- FeaturePlot(seurat, reduction = "umap", feature = "percent.rb")
    ggsave(UMAPrb, filename = paste0(outdir,"UMAPrb.png"))
    UMAPlib <- FeaturePlot(seurat, reduction = "umap", feature = "nCount_RNA")
    ggsave(UMAPlib, filename = paste0(outdir,"UMAPlib.png"))
    UMAPgenes <- FeaturePlot(seurat, reduction = "umap", feature = "nFeature_RNA")
    ggsave(UMAPgenes, filename = paste0(outdir,"UMAPgenes.png"))
    UMAPdonor_souporcell <- DimPlot(seurat, reduction = "umap", group.by = "common_assignment_souporcell")
    ggsave(UMAPdonor_souporcell, filename = paste0(outdir,"UMAPdonor_souporcell.png"))
    UMAPdonor_scSplit <- DimPlot(seurat, reduction = "umap", group.by = "common_assignment_scSplit")
    ggsave(UMAPdonor_scSplit, filename = paste0(outdir,"UMAPdonor_scSplit.png"))
    UMAPdonor_freemuxlet <- DimPlot(seurat, reduction = "umap", group.by = "common_assignment_freemuxlet")
    ggsave(UMAPdonor_freemuxlet, filename = paste0(outdir,"UMAPdonor_freemuxlet.png"))


} else if (stage[1]=="QCfigures_Harmony_Day"){
    seurat <- readRDS(paste0(dir,"output/Seurat/Harmony_Day/NoOutliers_Day_harmony.rds"))
    UMAPpool <- DimPlot(seurat, reduction = "umap", group.by = "Pool")
    ggsave(UMAPpool, filename = paste0(outdir,"UMAPpool.png"))
    UMAPmt <- FeaturePlot(seurat, reduction = "umap", feature = "percent.mt")
    ggsave(UMAPmt, filename = paste0(outdir,"UMAPmt.png"))
    UMAPrb <- FeaturePlot(seurat, reduction = "umap", feature = "percent.rb")
    ggsave(UMAPrb, filename = paste0(outdir,"UMAPrb.png"))
    UMAPlib <- FeaturePlot(seurat, reduction = "umap", feature = "nCount_RNA")
    ggsave(UMAPlib, filename = paste0(outdir,"UMAPlib.png"))
    UMAPgenes <- FeaturePlot(seurat, reduction = "umap", feature = "nFeature_RNA")
    ggsave(UMAPgenes, filename = paste0(outdir,"UMAPgenes.png"))
    UMAPdonor_souporcell <- DimPlot(seurat, reduction = "umap", group.by = "common_assignment_souporcell")
    ggsave(UMAPdonor_souporcell, filename = paste0(outdir,"UMAPdonor_souporcell.png"))
    UMAPdonor_scSplit <- DimPlot(seurat, reduction = "umap", group.by = "common_assignment_scSplit")
    ggsave(UMAPdonor_scSplit, filename = paste0(outdir,"UMAPdonor_scSplit.png"))
    UMAPdonor_freemuxlet <- DimPlot(seurat, reduction = "umap", group.by = "common_assignment_freemuxlet")
    ggsave(UMAPdonor_freemuxlet, filename = paste0(outdir,"UMAPdonor_freemuxlet.png"))


} else if (stage[1]=="QCfigures_Harmony_DayPatient"){
    seurat <- readRDS(paste0(dir,"output/Seurat/Harmony_DayPatient/NoOutliers_DayPatient_harmony.rds"))
    UMAPpool <- DimPlot(seurat, reduction = "umap", group.by = "Pool")
    ggsave(UMAPpool, filename = paste0(outdir,"UMAPpool.png"))
    UMAPmt <- FeaturePlot(seurat, reduction = "umap", feature = "percent.mt")
    ggsave(UMAPmt, filename = paste0(outdir,"UMAPmt.png"))
    UMAPrb <- FeaturePlot(seurat, reduction = "umap", feature = "percent.rb")
    ggsave(UMAPrb, filename = paste0(outdir,"UMAPrb.png"))
    UMAPlib <- FeaturePlot(seurat, reduction = "umap", feature = "nCount_RNA")
    ggsave(UMAPlib, filename = paste0(outdir,"UMAPlib.png"))
    UMAPgenes <- FeaturePlot(seurat, reduction = "umap", feature = "nFeature_RNA")
    ggsave(UMAPgenes, filename = paste0(outdir,"UMAPgenes.png"))
    UMAPdonor_souporcell <- DimPlot(seurat, reduction = "umap", group.by = "FinalAssignment")
    ggsave(UMAPdonor_souporcell, filename = paste0(outdir,"UMAPdonor_souporcell.png"))



} else if (stage[1]=="QCfigures_Integration_DayPatient"){
    seurat <- readRDS(paste0(dir,"output/Seurat/Integration_DayPatient/NoOutliers_DayPatient_SeuratIntegratee.rds"))
    UMAPpool <- DimPlot(seurat, reduction = "umap", group.by = "Pool")
    ggsave(UMAPpool, filename = paste0(outdir,"UMAPpool.png"))
    UMAPmt <- FeaturePlot(seurat, reduction = "umap", feature = "percent.mt")
    ggsave(UMAPmt, filename = paste0(outdir,"UMAPmt.png"))
    UMAPrb <- FeaturePlot(seurat, reduction = "umap", feature = "percent.rb")
    ggsave(UMAPrb, filename = paste0(outdir,"UMAPrb.png"))
    UMAPlib <- FeaturePlot(seurat, reduction = "umap", feature = "nCount_RNA")
    ggsave(UMAPlib, filename = paste0(outdir,"UMAPlib.png"))
    UMAPgenes <- FeaturePlot(seurat, reduction = "umap", feature = "nFeature_RNA")
    ggsave(UMAPgenes, filename = paste0(outdir,"UMAPgenes.png"))
    UMAPdonor_souporcell <- DimPlot(seurat, reduction = "umap", group.by = "FinalAssignment")
    ggsave(UMAPdonor_souporcell, filename = paste0(outdir,"UMAPdonor_souporcell.png"))


} else if (stage[1]=="QCfigures_Integration_Day_then_Individual"){
    seurat <- readRDS(paste0(dir,"output/Seurat/Integration_Day_then_Individual/NoOutliers_Day_SeuratIntegratee.rds"))
    UMAPpool <- DimPlot(seurat, reduction = "umap", group.by = "Pool")
    ggsave(UMAPpool, filename = paste0(outdir,"UMAPpool.png"))
    UMAPmt <- FeaturePlot(seurat, reduction = "umap", feature = "percent.mt")
    ggsave(UMAPmt, filename = paste0(outdir,"UMAPmt.png"))
    UMAPrb <- FeaturePlot(seurat, reduction = "umap", feature = "percent.rb")
    ggsave(UMAPrb, filename = paste0(outdir,"UMAPrb.png"))
    UMAPlib <- FeaturePlot(seurat, reduction = "umap", feature = "nCount_RNA")
    ggsave(UMAPlib, filename = paste0(outdir,"UMAPlib.png"))
    UMAPgenes <- FeaturePlot(seurat, reduction = "umap", feature = "nFeature_RNA")
    ggsave(UMAPgenes, filename = paste0(outdir,"UMAPgenes.png"))
    UMAPdonor_souporcell <- DimPlot(seurat, reduction = "umap", group.by = "FinalAssignment")
    ggsave(UMAPdonor_souporcell, filename = paste0(outdir,"UMAPdonor_souporcell.png"))





} else if (stage[1]=="Compare_Harmony_Seurat_Clusters"){
    dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Seurat/"
    outdir <- paste0(dir,"Harmony_vs_Integration/")
    dir.create(outdir)
    seurat_harmony <- readRDS("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Seurat/Harmony_DayPatient/NoOutliers_DayPatient_harmony.rds")
    seurat_integration <- readRDS("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Seurat/Integration_Day_then_Individual/NoOutliers_Day_SeuratIntegratee.rds")
    

    harmony_resolution_file_list <- list.files(paste0(dir,"/Multi_Resolution_Harmony_DayPatient/"), pattern = "Identities_Resolution_")
    seurat_harmony_res_list <- lapply(harmony_resolution_file_list, FUN = function(x){
        read_delim(paste0(dir,"Multi_Resolution_Harmony_DayPatient/",x), delim = "\t")
    })
    
    harmony_ResolutionsDF <- do.call(cbind, seurat_harmony_res_list)
    colnames(harmony_ResolutionsDF) <- paste0("Harmony_", colnames(harmony_ResolutionsDF))
    rownames(harmony_ResolutionsDF) <- colnames(seurat_harmony)
    seurat_integration <- AddMetaData(seurat_integration, harmony_ResolutionsDF)

    for (res in colnames(seurat_integration@meta.data)[grep("Resolution", colnames(seurat_integration@meta.data))]){
        UMAP <- DimPlot(seurat_integration, reduction = "umap", group.by = res)
        ggsave(UMAP, filename = paste0(outdir,"seuratUMAP_harmony",res,".png"))
    }


    integration_resolution_file_list <- list.files(paste0(dir,"Multi_Resolution_Integration_Day_then_Individual/"), pattern = "Identities_Resolution_")
    seurat_integration_res_list <- lapply(integration_resolution_file_list, FUN = function(x){
        read_delim(paste0(dir,"Multi_Resolution_Integration_Day_then_Individual/",x), delim = "\t")
    })
    
    integration_ResolutionsDF <- do.call(cbind, seurat_integration_res_list)
    colnames(integration_ResolutionsDF) <- paste0("Seurat_",colnames(integration_ResolutionsDF))
    rownames(integration_ResolutionsDF) <- colnames(seurat_integration)
    seurat_harmony <- AddMetaData(seurat_harmony, integration_ResolutionsDF)

    for (res in colnames(seurat_harmony@meta.data)[grep("Resolution", colnames(seurat_harmony@meta.data))]){
        UMAP <- DimPlot(seurat_harmony, reduction = "umap", group.by = res)
        ggsave(UMAP, filename = paste0(outdir,"harmonyUMAP_seurat",res,".png"))
    }

    integration_ResolutionsDF$Barcode <- rownames(integration_ResolutionsDF)
    harmony_ResolutionsDF$Barcode <- rownames(harmony_ResolutionsDF)
    ResolutionsDF <- left_join(integration_ResolutionsDF,harmony_ResolutionsDF)

    plotlist <- list()

for (col in 1:ncol(ResolutionsDF)){
    ResolutionsDF[,col] <- as.factor(ResolutionsDF[,col])
}

for (res in seq(0,0.24,0.01)){
    seurat <- paste0("Seurat_Resolution_",res)
    harmony <- paste0("Harmony_Resolution_",res)
            plotlist[[paste0(seurat,harmony)]] <- ggplot(ResolutionsDF, aes_string(seurat)) +
                        geom_bar(position = "stack", aes_string(fill = harmony)) +
                        theme_classic()

}
    lapply(names(plotlist), function(x){ggsave(plotlist[[x]], filename = paste0(outdir,x,".png"))})
    
    lapply(plotlist1, function(x){ggsave(x, filename = paste0(outdir,"seuratClusters_fillHarmony_",res1,"_",res2,".png"))})
    lapply(plotlist2, function(x){ggsave(x, filename = paste0(outdir,"HarmonyClusters_fillSeurat_",res1,"_",res2,".png"))})




}


