library(Seurat)
library(data.table)
library(dplyr)
library(tidyverse)
library(colorspace)
library(Nebulosa)


##### Set up directories #####
dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
datadir <- paste0(dir, "data/multi-passage/GIMR_GWCCG_211213_JOSPOW_10x_3p_results/GIMR_GWCCG_211213_JOSPOW_10x_3p/analyses/")
outdir <- paste0(dir, "output/multi-passage/preQC/")



##### Set up variables #####
pools <- dir(datadir)


## Read in expression data
counts_list <- lapply(pools, function(x){
    print(x)
    Read10X(paste0(datadir,x, "/GE/",x,"/outs/per_sample_outs/",x,"/count/sample_feature_bc_matrix"), gene.column = 1)
})
names(counts_list) <- pools


## Add poolnames to cell names so can easily match back if there are any cells with same barcoe
counts_list <- lapply(names(counts_list), function(x){
    colnames(counts_list[[x]]) <- gsub("-1", "", paste0(x, "_", colnames(counts_list[[x]])))
    return(counts_list[[x]])
})
names(counts_list) <- pools


## Make seurat object ##
seurat_list <- lapply(counts_list, function(x){
    CreateSeuratObject(counts = x)
})
names(seurat_list) <- pools


seurat <- merge(seurat_list[[1]], y = seurat_list[2:length(seurat_list)],project = "Village_Phase1_multi-passage")



### Add gene ID data to gene level metadata ###
features <- fread(paste0(datadir,pools[1], "/GE/",pools[1],"/outs/per_sample_outs/",pools[1],"/count/sample_feature_bc_matrix/features.tsv.gz"), col.names = c("ENSG", "Gene_ID", "Assay"))
features$Assay <- NULL

features_df <- data.frame(features)
rownames(features_df) <- features$ENSG
features_df$ENSG <- NULL

seurat[["RNA"]] <- AddMetaData(seurat[["RNA"]], features_df)


## Add QC metrics
RbGeneList <- read.delim(file = "/directflow/SCCGGroupShare/projects/DrewNeavin/References/RibosomalGeneList_GeneID_ENSG.txt")
MtGeneList <- read.delim(file = "/directflow/SCCGGroupShare/projects/DrewNeavin/References/MtGeneList_GeneID_ENSG.txt")
print("Calculating Mt %")
seurat <- PercentageFeatureSet(seurat, features = MtGeneList$ENSG, col.name = "percent.mt")
print("Calculating Rb %")
RbGeneList <- RbGeneList[which(RbGeneList$ENSG %in% rownames(seurat)),]
seurat <- PercentageFeatureSet(seurat, features = RbGeneList$ENSG, col.name = "percent.rb")


### Make pre-QC figures ###
seurat <- NormalizeData(seurat, verbose = TRUE)
seurat <- FindVariableFeatures(seurat, selection.method = "mean.var.plot")
seurat <- ScaleData(seurat, features = VariableFeatures(seurat))
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))
seurat <- FindNeighbors(seurat, dims = 1:10)
seurat <- FindClusters(seurat, resolution = 0.5)
seurat <- RunUMAP(seurat, dims = 1:10)


seurat$Pool <- gsub("_[ATCG]+", "", colnames(seurat))


### QC Figures pre filtering ###
plot_mt_pct <- VlnPlot(seurat, features = c( "percent.mt"), group.by = "Pool", pt.size = 0) +
                    scale_fill_discrete_sequential(palette = "SunsetDark")
ggsave(plot_mt_pct, filename = paste0(outdir,"Mt_pct_vln.png"))

plot_rb_pct <- VlnPlot(seurat, features = c( "percent.rb"), group.by = "Pool", pt.size = 0) +
                    scale_fill_discrete_sequential(palette = "SunsetDark")
ggsave(plot_rb_pct, filename = paste0(outdir,"Rb_pct_vln.png"))

plot_n_count <- VlnPlot(seurat, features = c( "nCount_RNA"), group.by = "Pool", pt.size = 0) +
                    scale_fill_discrete_sequential(palette = "SunsetDark")
ggsave(plot_n_count, filename = paste0(outdir,"N_count_vln.png"))

plot_nFeature_RNA <- VlnPlot(seurat, features = c( "nFeature_RNA"), group.by = "Pool", pt.size = 0) +
                    scale_fill_discrete_sequential(palette = "SunsetDark")
ggsave(plot_nFeature_RNA, filename = paste0(outdir,"nFeature_RNA_vln.png"))

lib_mt <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Pool") +
                    scale_color_discrete_sequential(palette = "SunsetDark", alpha = 0.25)
ggsave(lib_mt, filename = paste0(outdir,"lib_mt.png"))

lib_genes <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Pool") +
                    scale_color_discrete_sequential(palette = "SunsetDark", alpha = 0.25)
ggsave(lib_genes, filename = paste0(outdir,"lib_genes.png"))

nebulosa_mt_umap <- plot_density(seurat, "percent.mt", pal = "plasma")
ggsave(nebulosa_mt_umap, filename = paste0(outdir,"mt_percent_umap.png"))

nebulosa_rb_umap <- plot_density(seurat, "percent.rb", pal = "plasma")
ggsave(nebulosa_rb_umap, filename = paste0(outdir,"rb_percent_umap.png"))

umap_Pool <- DimPlot(seurat, group.by = "Pool") +
                    scale_color_discrete_sequential(palette = "SunsetDark")
ggsave(umap_Pool, filename = paste0(outdir,"pool_umap.png"))

mt_umap <- FeaturePlot(seurat, features = "percent.mt") +
                    scale_color_continuous_sequential(palette = "RedPurple")
ggsave(mt_umap, filename = paste0(outdir,"mt_percent_umap_seurat.png"))

