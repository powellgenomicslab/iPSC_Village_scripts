library(Seurat)
library(data.table)
library(dplyr)
library(tidyverse)
library(colorspace)
library(Nebulosa)
library(RColorBrewer)
library(scran)
library(ggbreak) 


##### Set up directories #####
dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
datadir <- paste0(dir, "data/multi-passage/resequencing/220405/")
outdir <- paste0(dir, "output/multi-passage/preQC/")



##### Set up variables #####
pools <- dir(datadir)


## Read in expression data
counts_list <- lapply(pools, function(x){
    print(x)
    Read10X(paste0(datadir,x, "/GE/summary/sample_feature_bc_matrix/"), gene.column = 1)
})
names(counts_list) <- pools


## Add poolnames to cell names so can easily match back if there are any cells with same barcoe
counts_list <- lapply(names(counts_list), function(x){
    colnames(counts_list[[x]]) <- gsub("-1", "", paste0(x, "_", colnames(counts_list[[x]])))
    return(counts_list[[x]])
})
names(counts_list) <- pools


## Make seurat object ##
seurat_list <- lapply(names(counts_list), function(x){
    temp <- CreateSeuratObject(counts = counts_list[[x]])
    temp@meta.data$Pool <- x
    return(temp)
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


## Add Demultiplexing Results ##
demultiplexing_list <- lapply(pools, function(pool){
    temp <- data.frame(fread(paste0(dir, "output/multi-passage/demultiplexed/with_demuxlet/updated_2022_06_26/atleasthalf_singlet/", pool, "/atleasthalf_singlet_w_combined_assignments.tsv")))
    rownames(temp) <- gsub("-1", "", paste0(pool, "_", temp$Barcode))
    return(temp)
})

demultiplexing_dt <- do.call(rbind,demultiplexing_list)


seurat <- AddMetaData(seurat, demultiplexing_dt)


print("Starting cell cycle step")
hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))

countsENSG <- seurat[["RNA"]]@counts # Create a counts matrix that has the ENSG gene clasifiers so can determine cell cycle phase
print("The size of the counts object in genes x cells is:")
print(dim(countsENSG))

### Run the cell cycle identification ###
print("Starting cell cycle determination")
assigned <- cyclone(countsENSG, pairs=hs.pairs)  #Note, this takes hours
table(assigned$phases)
write.table(assigned, file = paste0(outdir,"CellCycleProportions.txt"), quote = F, sep = "\t") #Save so that can read in and don't have to wait to recompute again


assigned <- read.table(paste0(outdir,"CellCycleProportions.txt"), sep = "\t")
assigned <- as.data.frame(assigned)
rownames(assigned) <- colnames(seurat)
write.table(assigned, file = paste0(outdir,"CellCycleProportions.txt"), quote = F, sep = "\t") #Save so that can read in and don't have to wait to recompute again
assigned <- read.table(paste0(outdir,"CellCycleProportions.txt"), sep = "\t")

AddMetaData(seurat, assigned)
saveRDS(seurat, paste0(outdir,"seurat_all_cells_cell_cycle.rds"))
seurat <- readRDS(paste0(outdir,"seurat_all_cells_cell_cycle.rds"))



### Make pre-QC figures ###
seurat <- NormalizeData(seurat, verbose = TRUE)
seurat <- FindVariableFeatures(seurat, selection.method = "mean.var.plot")
seurat <- ScaleData(seurat, features = VariableFeatures(seurat))
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))
seurat <- FindNeighbors(seurat, dims = 1:10)
seurat <- FindClusters(seurat, resolution = 0.5)
seurat <- RunUMAP(seurat, dims = 1:10)



### QC Figures pre filtering ###
plot_mt_pct <- VlnPlot(seurat, features = c( "percent.mt"), group.by = "Pool", pt.size = 0) +
                    scale_fill_discrete_sequential(palette = "SunsetDark") +
                    theme(plot.background = element_rect(fill = "white"))
ggsave(plot_mt_pct, filename = paste0(outdir,"Mt_pct_vln.png"))

plot_rb_pct <- VlnPlot(seurat, features = c( "percent.rb"), group.by = "Pool", pt.size = 0) +
                    scale_fill_discrete_sequential(palette = "SunsetDark") +
                    theme(plot.background = element_rect(fill = "white"))
ggsave(plot_rb_pct, filename = paste0(outdir,"Rb_pct_vln.png"))

plot_n_count <- VlnPlot(seurat, features = c( "nCount_RNA"), group.by = "Pool", pt.size = 0) +
                    scale_fill_discrete_sequential(palette = "SunsetDark") +
                    theme(plot.background = element_rect(fill = "white"))
ggsave(plot_n_count, filename = paste0(outdir,"N_count_vln.png"))

plot_nFeature_RNA <- VlnPlot(seurat, features = c( "nFeature_RNA"), group.by = "Pool", pt.size = 0) +
                    scale_fill_discrete_sequential(palette = "SunsetDark") +
                    theme(plot.background = element_rect(fill = "white"))
ggsave(plot_nFeature_RNA, filename = paste0(outdir,"nFeature_RNA_vln.png"))

lib_mt <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Pool") +
                    scale_color_discrete_sequential(palette = "SunsetDark", alpha = 0.25) +
                    theme(plot.background = element_rect(fill = "white"))
ggsave(lib_mt, filename = paste0(outdir,"lib_mt.png"))

lib_genes <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Pool") +
                    scale_color_discrete_sequential(palette = "SunsetDark", alpha = 0.25) +
                    theme(plot.background = element_rect(fill = "white"))
ggsave(lib_genes, filename = paste0(outdir,"lib_genes.png"))

nebulosa_mt_umap <- plot_density(seurat, "percent.mt", pal = "plasma") +
                    theme(plot.background = element_rect(fill = "white"))
ggsave(nebulosa_mt_umap, filename = paste0(outdir,"mt_percent_umap.png"))

nebulosa_rb_umap <- plot_density(seurat, "percent.rb", pal = "plasma") +
                    theme(plot.background = element_rect(fill = "white"))
ggsave(nebulosa_rb_umap, filename = paste0(outdir,"rb_percent_umap.png"))

umap_Pool <- DimPlot(seurat, group.by = "Pool") +
                    scale_color_discrete_sequential(palette = "SunsetDark") +
                    theme(plot.background = element_rect(fill = "white"))
ggsave(umap_Pool, filename = paste0(outdir,"pool_umap.png"))

mt_umap <- FeaturePlot(seurat, features = "percent.mt") +
                    scale_color_continuous_sequential(palette = "RedPurple") +
                    theme(plot.background = element_rect(fill = "white"))
ggsave(mt_umap, filename = paste0(outdir,"mt_percent_umap_seurat.png"))

colourCount = length(unique(seurat@meta.data$AtLeastHalfSinglet_Individual_Assignment))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

umap_individual <- DimPlot(seurat, group.by = "AtLeastHalfSinglet_Individual_Assignment") +
                    scale_color_manual(values = getPalette(colourCount)) +
                    theme(plot.background = element_rect(fill = "white"))
ggsave(umap_individual, filename = paste0(outdir,"pool_individual.png"))

umap_droplet_type <- DimPlot(seurat, group.by = "AtLeastHalfSinglet_DropletType") +
                    scale_color_discrete_sequential(palette = "Red-Blue", rev = FALSE) +
                    theme(plot.background = element_rect(fill = "white"))
ggsave(umap_droplet_type, filename = paste0(outdir,"pool_droplet_type.png"))

umap_phases <- DimPlot(seurat, group.by = "phases") +
                    scale_color_discrete_sequential(palette = "Sunset", rev = FALSE) +
                    theme(plot.background = element_rect(fill = "white"))
ggsave(umap_phases, filename = paste0(outdir,"cell_cycle_type.png"))



##### Subset seurat object by cells that have no info for cell assignments (removed by dropletQC) #####
seurat_sub <- subset(seurat, subset = AtLeastHalfSinglet_Individual_Assignment != "doublet")



### Make pre-QC figures ###
seurat_sub <- NormalizeData(seurat_sub, verbose = TRUE)
seurat_sub <- FindVariableFeatures(seurat_sub, selection.method = "mean.var.plot")
seurat_sub <- ScaleData(seurat_sub, features = VariableFeatures(seurat_sub))
seurat_sub <- RunPCA(seurat_sub, features = VariableFeatures(object = seurat_sub))
seurat_sub <- FindNeighbors(seurat_sub, dims = 1:10)
seurat_sub <- FindClusters(seurat_sub, resolution = 0.5)
seurat_sub <- RunUMAP(seurat_sub, dims = 1:10)



### QC Figures pre filtering ###
plot_mt_pct <- VlnPlot(seurat_sub, features = c( "percent.mt"), group.by = "Pool", pt.size = 0) +
                    scale_fill_discrete_sequential(palette = "SunsetDark") +
                    theme(plot.background = element_rect(fill = "white"))
ggsave(plot_mt_pct, filename = paste0(outdir,"Mt_pct_vln_singlets.png"))

plot_rb_pct <- VlnPlot(seurat_sub, features = c( "percent.rb"), group.by = "Pool", pt.size = 0) +
                    scale_fill_discrete_sequential(palette = "SunsetDark") +
                    theme(plot.background = element_rect(fill = "white"))
ggsave(plot_rb_pct, filename = paste0(outdir,"Rb_pct_vln_singlets.png"))

plot_n_count <- VlnPlot(seurat_sub, features = c( "nCount_RNA"), group.by = "Pool", pt.size = 0) +
                    scale_fill_discrete_sequential(palette = "SunsetDark") +
                    theme(plot.background = element_rect(fill = "white"))
ggsave(plot_n_count, filename = paste0(outdir,"N_count_vln_singlets.png"))

plot_nFeature_RNA <- VlnPlot(seurat_sub, features = c( "nFeature_RNA"), group.by = "Pool", pt.size = 0) +
                    scale_fill_discrete_sequential(palette = "SunsetDark") +
                    theme(plot.background = element_rect(fill = "white"))
ggsave(plot_nFeature_RNA, filename = paste0(outdir,"nFeature_RNA_vln_singlets.png"))

lib_mt <- FeatureScatter(seurat_sub, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Pool") +
                    scale_color_discrete_sequential(palette = "SunsetDark", alpha = 0.25) +
                    theme(plot.background = element_rect(fill = "white"))
ggsave(lib_mt, filename = paste0(outdir,"lib_mt_singlets.png"))

lib_genes <- FeatureScatter(seurat_sub, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Pool") +
                    scale_color_discrete_sequential(palette = "SunsetDark", alpha = 0.25) +
                    theme(plot.background = element_rect(fill = "white"))
ggsave(lib_genes, filename = paste0(outdir,"lib_genes_singlets.png"))

nebulosa_mt_umap <- plot_density(seurat_sub, "percent.mt", pal = "plasma") +
                    theme(plot.background = element_rect(fill = "white"))
ggsave(nebulosa_mt_umap, filename = paste0(outdir,"mt_percent_umap_singlets.png"))

nebulosa_rb_umap <- plot_density(seurat_sub, "percent.rb", pal = "plasma") +
                    theme(plot.background = element_rect(fill = "white"))
ggsave(nebulosa_rb_umap, filename = paste0(outdir,"rb_percent_umap_singlets.png"))

umap_Pool <- DimPlot(seurat_sub, group.by = "Pool") +
                    scale_color_discrete_sequential(palette = "SunsetDark") +
                    theme(plot.background = element_rect(fill = "white"))
ggsave(umap_Pool, filename = paste0(outdir,"pool_umap_singlets.png"))

mt_umap <- FeaturePlot(seurat_sub, features = "percent.mt") +
                    scale_color_continuous_sequential(palette = "RedPurple") +
                    theme(plot.background = element_rect(fill = "white"))
ggsave(mt_umap, filename = paste0(outdir,"mt_percent_umap_seurat_singlets.png"))

colourCount = length(unique(seurat_sub@meta.data$AtLeastHalfSinglet_Individual_Assignment))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

umap_individual <- DimPlot(seurat_sub, group.by = "AtLeastHalfSinglet_Individual_Assignment") +
                    scale_color_manual(values = getPalette(colourCount)) +
                    theme(plot.background = element_rect(fill = "white"))
ggsave(umap_individual, filename = paste0(outdir,"pool_individual_singlets.png"))

umap_droplet_type <- DimPlot(seurat_sub, group.by = "AtLeastHalfSinglet_DropletType") +
                    scale_color_discrete_sequential(palette = "Red-Blue", rev = FALSE) +
                    theme(plot.background = element_rect(fill = "white"))
ggsave(umap_droplet_type, filename = paste0(outdir,"pool_droplet_type_singlets.png"))

umap_phases <- DimPlot(seurat_sub, group.by = "phases") +
                    scale_color_discrete_sequential(palette = "Sunset", rev = FALSE) +
                    theme(plot.background = element_rect(fill = "white"))
ggsave(umap_phases, filename = paste0(outdir,"cell_cycle_type_singlets.png"))





##### Remove high mitochondrial % cells #####
seurat_sub_mt <- subset(seurat_sub, subset = percent.mt < 25)


### Integrate different time-points
seurat_sub_mt_list <- SplitObject(seurat_sub_mt, split.by = "Pool")
seurat_sub_mt_list <- lapply(X = seurat_sub_mt_list, function(x) SCTransform(x, vars.to.regress = c("scores.G1", "scores.S", "scores.G2M", "percent.mt", "percent.rb"), return.only.var.genes = FALSE))
# features <- SelectIntegrationFeatures(object.list = seurat_sub_mt_list, nfeatures = 3000)
features <- Reduce(intersect, lapply(seurat_sub_mt_list, rownames))  
seurat_sub_mt_list_prep <- PrepSCTIntegration(object.list = seurat_sub_mt_list, anchor.features = features)

saveRDS(seurat_sub_mt_list, paste0(outdir,"seurat_sub_mt_list.rds"))
seurat_sub_mt_list <- readRDS(paste0(outdir,"seurat_sub_mt_list.rds"))

saveRDS(seurat_sub_mt_list_prep, paste0(outdir,"seurat_sub_mt_list_prep.rds"))
seurat_sub_mt_list_prep <- readRDS(paste0(outdir,"seurat_sub_mt_list_prep.rds"))


### Check if prepsctintegration impact the sct scale.data, if not, just merge them together

seurat_sub_mt_anchors <- FindIntegrationAnchors(object.list = seurat_sub_mt_list_prep, normalization.method = "SCT",
    anchor.features = features)
combined_sct <- IntegrateData(anchorset = seurat_sub_mt_anchors, normalization.method = "SCT")

combined_sct <- RunPCA(combined_sct, verbose = FALSE)
combined_sct <- RunUMAP(combined_sct, reduction = "pca", dims = 1:30)


### QC Figures  ###
plot_mt_pct <- VlnPlot(combined_sct, features = c( "percent.mt"), group.by = "Pool", pt.size = 0) +
                    scale_fill_discrete_sequential(palette = "SunsetDark") +
                    theme(plot.background = element_rect(fill = "white"))
ggsave(plot_mt_pct, filename = paste0(outdir,"Mt_pct_vln_singlets_filtered.png"))

plot_rb_pct <- VlnPlot(combined_sct, features = c( "percent.rb"), group.by = "Pool", pt.size = 0) +
                    scale_fill_discrete_sequential(palette = "SunsetDark") +
                    theme(plot.background = element_rect(fill = "white"))
ggsave(plot_rb_pct, filename = paste0(outdir,"Rb_pct_vln_singlets_filtered.png"))

plot_n_count <- VlnPlot(combined_sct, features = c( "nCount_RNA"), group.by = "Pool", pt.size = 0) +
                    scale_fill_discrete_sequential(palette = "SunsetDark") +
                    theme(plot.background = element_rect(fill = "white"))
ggsave(plot_n_count, filename = paste0(outdir,"N_count_vln_singlets_filtered.png"))

plot_nFeature_RNA <- VlnPlot(combined_sct, features = c( "nFeature_RNA"), group.by = "Pool", pt.size = 0) +
                    scale_fill_discrete_sequential(palette = "SunsetDark") +
                    theme(plot.background = element_rect(fill = "white"))
ggsave(plot_nFeature_RNA, filename = paste0(outdir,"nFeature_RNA_vln_singlets_filtered.png"))

lib_mt <- FeatureScatter(combined_sct, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Pool") +
                    scale_color_discrete_sequential(palette = "SunsetDark", alpha = 0.25) +
                    theme(plot.background = element_rect(fill = "white"))
ggsave(lib_mt, filename = paste0(outdir,"lib_mt_singlets_filtered.png"))

lib_genes <- FeatureScatter(combined_sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Pool") +
                    scale_color_discrete_sequential(palette = "SunsetDark", alpha = 0.25) +
                    theme(plot.background = element_rect(fill = "white"))
ggsave(lib_genes, filename = paste0(outdir,"lib_genes_singlets_filtered.png"))

nebulosa_mt_umap <- plot_density(combined_sct, "percent.mt", pal = "plasma") +
                    theme(plot.background = element_rect(fill = "white"))
ggsave(nebulosa_mt_umap, filename = paste0(outdir,"mt_percent_umap_singlets_filtered.png"))

nebulosa_rb_umap <- plot_density(combined_sct, "percent.rb", pal = "plasma") +
                    theme(plot.background = element_rect(fill = "white"))
ggsave(nebulosa_rb_umap, filename = paste0(outdir,"rb_percent_umap_singlets_filtered.png"))

umap_Pool <- DimPlot(combined_sct, group.by = "Pool") +
                    scale_color_discrete_sequential(palette = "SunsetDark") +
                    theme(plot.background = element_rect(fill = "white"))
ggsave(umap_Pool, filename = paste0(outdir,"pool_umap_singlets_filtered.png"))

mt_umap <- FeaturePlot(combined_sct, features = "percent.mt") +
                    scale_color_continuous_sequential(palette = "RedPurple") +
                    theme(plot.background = element_rect(fill = "white"))
ggsave(mt_umap, filename = paste0(outdir,"mt_percent_umap_seurat_singlets_filtered.png"))

colourCount = length(unique(combined_sct@meta.data$AtLeastHalfSinglet_Individual_Assignment))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

umap_individual <- DimPlot(combined_sct, group.by = "AtLeastHalfSinglet_Individual_Assignment") +
                    scale_color_manual(values = getPalette(colourCount)) +
                    theme(plot.background = element_rect(fill = "white"))
ggsave(umap_individual, filename = paste0(outdir,"pool_individual_singlets_filtered.png"))

umap_droplet_type <- DimPlot(combined_sct, group.by = "AtLeastHalfSinglet_DropletType") +
                    scale_color_discrete_sequential(palette = "Red-Blue", rev = FALSE) +
                    theme(plot.background = element_rect(fill = "white"))
ggsave(umap_droplet_type, filename = paste0(outdir,"pool_droplet_type_singlets_filtered.png"))

umap_phases <- DimPlot(combined_sct, group.by = "phases") +
                    scale_color_discrete_sequential(palette = "Sunset", rev = FALSE) +
                    theme(plot.background = element_rect(fill = "white"))
ggsave(umap_phases, filename = paste0(outdir,"cell_cycle_type_singlets_filtered.png"))


### Update individual IDs to match what used for Nona's plots ###
combined_sct$Assignment <- gsub("^0_", "", combined_sct$AtLeastHalfSinglet_Individual_Assignment) %>%
                                gsub("D-", "", .) %>%
                                    gsub("\\.\\d\\.", "", .) %>%
                                        gsub("N-", "", .) %>%
                                            gsub("-P36", "", .)  %>%
                                                gsub("-", "", .)

seurat_sub$Assignment <- gsub("^0_", "", seurat_sub$AtLeastHalfSinglet_Individual_Assignment) %>%
                                gsub("D-", "", .) %>%
                                    gsub("\\.\\d\\.", "", .) %>%
                                        gsub("N-", "", .) %>%
                                            gsub("-P36", "", .)  %>%
                                                gsub("-", "", .)


##### Make proportion plots (area plot) #####
village_summary <- data.table(prop.table(table(seurat_sub@meta.data[,c("Assignment", "Pool")]), margin = 2))
village_summary$Assignment <- factor(village_summary$Assignment, levels = rev(village_summary[Pool == "Village_P8"]$Assignment[order(village_summary[Pool == "Village_P8"]$N)]))

colors_original <- readRDS("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Nona_multiome/line_colors")
colors_original <- colors_original[levels(village_summary$Assignment)]


p_stacked_area <- ggplot(village_summary, aes(x = as.numeric(as.character(gsub("Village_P","", Pool))), y = N, fill = factor(Assignment), group = Assignment)) +
    geom_area(alpha=0.6 , size=0.5, colour="black") +
    theme_classic() +
    scale_fill_manual(values = colors_original) +
    xlab("Passage") +
    ylab("Proportion of Cells")
ggsave(p_stacked_area, filename = paste0(outdir,"stacked_area.png"), width = 7, height = 4)
ggsave(p_stacked_area, filename = paste0(outdir,"stacked_area.pdf"), width = 7, height = 4)


village_summary_singlets <- data.table(prop.table(table(combined_sct@meta.data[,c("Assignment", "Pool")]), margin = 2))
village_summary_singlets$Assignment <- factor(village_summary_singlets$Assignment, levels = rev(village_summary_singlets[Pool == "Village_P8"]$Assignment[order(village_summary_singlets[Pool == "Village_P8"]$N)]))
village_summary_singlets$Pool <- as.numeric(gsub("Village_P", "", village_summary_singlets$Pool))
colnames(village_summary_singlets) <- gsub("Pool", "Day", colnames(village_summary_singlets))



colors <- readRDS("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Nona_multiome/line_colors")
colors <- colors[levels(village_summary_singlets$Assignment)]


p_stacked_area_filt <- ggplot(village_summary_singlets, aes(x = Day, y = N, fill = factor(Assignment), group = Assignment)) +
    geom_area(alpha=0.6 , size=0.1, colour="black") +
    theme_classic() +
    scale_fill_manual(values = colors) +
    xlab("Passage") +
    ylab("Proportion of Cells")

ggsave(p_stacked_area_filt, filename = paste0(outdir,"stacked_area_filtered.png"), width = 4, height = 2)
ggsave(p_stacked_area_filt, filename = paste0(outdir,"stacked_area_filtered.pdf"), width = 4, height = 2)






saveRDS(combined_sct, paste0(outdir,"time-integrated_filtered_seurat.rds"))
combined_sct <- readRDS(paste0(outdir,"time-integrated_filtered_seurat.rds"))



##### Read in Nona's original data and add pre-cryopreserved numbers with break #####
non_dir <- "/directflow/SCCGGroupShare/projects/nonfar/analysis/cardiac_multiome_directflow/demux_obj/"


##### Get a list of the village pools #####
villages <- list.files(non_dir, pattern = "Village")
villages <- grep("DemuxALL", villages, value = TRUE)



##### Get the singlets from the file #####
village_id_list <- lapply(villages, function(x){
    print(x)
    # tmp <- fread(paste0(datadir,x,"/CombinedResults/Final_Assignments_demultiplexing_doublets_new_edit.txt"), sep = "\t")
    tmp <- readRDS(paste0(non_dir,x))
    dt <- data.table(tmp@meta.data)
    dt$Pool_ID <- gsub("_DemuxALL.rds", "",x)
    dt$Day <- as.numeric(as.character(gsub("Village_Day", "", dt$Pool_ID)))
    return(dt)
})

village_id <- do.call(rbind, village_id_list)

village_id$Pool_ID_updated <- gsub("Day7$", "Day5b", village_id$Pool_ID) %>%
                                gsub("Day5$", "Day7b", .) %>% 
                                    gsub("Day15$", "Day4b", .) %>%
                                        gsub("Day4$", "Day15b", .) %>%
                                            gsub("b", "", .)

village_id$Day_updated <- as.numeric(as.character(gsub("Village_Day", "", village_id$Pool_ID_updated)))

### Update columnes ###
village_id$Day <- factor(village_id$Day, levels = c(0,1,2,3,4,5,7,15))
village_id$Assignment <- gsub("^0_", "", village_id$Assignment) %>%
                                gsub("D-", "", .) %>%
                                    gsub("\\.\\d\\.", "", .) %>%
                                        gsub("N-", "", .) %>%
                                            gsub("-P36", "", .)  %>%
                                                gsub("-", "", .)

village_id$DropletType <- ifelse(village_id$DropletType != "singlet", "doublet", village_id$DropletType)

data.table(prop.table(table(village_id[,c("DropletType", "Day_updated")]), margin = 2))

village_summary <- data.table(prop.table(table(village_id[,c("Assignment", "Day_updated")]), margin = 2))

village_summary_singlets_non <- data.table(prop.table(table(village_id[Assignment != "unassigned" & Assignment != "doublet",c("Assignment", "Day_updated")]), margin = 2))

village_summary_singlets_non$Assignment <- factor(village_summary_singlets$Assignment, levels = rev(village_summary_singlets[Day_updated == 15]$Assignment[order(village_summary_singlets[Day_updated == 15]$N)]))

colnames(village_summary_singlets_non) <- gsub("Day_updated", "Day", colnames(village_summary_singlets_non))

fwrite(village_summary_singlets_non, paste0(outdir, "cardiac_diff_prop_lines.tsv"), sep = "\t")
fwrite(village_summary_singlets, paste0(outdir, "multi_passage_prop_lines.tsv"), sep = "\t")

village_summary_singlets_non <- fread(paste0(outdir, "cardiac_diff_prop_lines.tsv"), sep = "\t")
village_summary_singlets <- fread( paste0(outdir, "multi_passage_prop_lines.tsv"), sep = "\t")

village_summary_combined <- rbind(village_summary_singlets_non[Day == 0], village_summary_singlets)
village_summary_combined$Day <- gsub(0,-1,village_summary_combined$Day)


p_stacked_area_filt_combined <- ggplot(village_summary_combined, aes(x = as.numeric(Day), y = N, fill = factor(Assignment, levels = names(colors)), group = factor(Assignment, levels = names(colors)))) +
    geom_area(alpha=0.6 , size=0.1, colour="black") +
    theme_classic() +
    scale_fill_manual(values = colors) +
    xlab("Passage") +
    ylab("Proportion of Cells") +
    theme(legend.position = "none")


ggsave(p_stacked_area_filt_combined, filename = paste0(outdir,"stacked_area_filtered_w_cryo.png"), width = 4, height = 2)
ggsave(p_stacked_area_filt_combined, filename = paste0(outdir,"stacked_area_filtered_w_cryo.pdf"), width = 4, height = 2)




##### Make Number at each time plot #####
village_summary_combined_n <- data.table(Day = unique(village_summary_combined$Day), N = as.numeric(NA))

for (day in village_summary_combined_n$Day){
    village_summary_combined_n[Day == day]$N <- nrow(village_summary_combined[Day == day & N > 0])
}


p_N_line <- ggplot(village_summary_combined_n, aes(x = as.numeric(as.character(Day)), y = N)) +
    geom_point(colour="black", size = 0.8) +
    geom_line() +
    # geom_smooth(method = "lm", se = FALSE) +
    theme_classic() +
    xlab("Passage") +
    ylab("# hiPSC\nLines") +
    ylim(0,18)

ggsave(p_N_line, filename = paste0(outdir,"line_N.png"), width = 4.4, height = 1.5)
ggsave(p_N_line, filename = paste0(outdir,"line_N.pdf"), width = 4.4, height = 1.5)






### Add number of cells at each time into metadata and 1/n for models
freq_dt <- as.data.frame(table(combined_sct$Assignment, combined_sct$Pool))
freq_dt$Ncov <- ifelse(freq_dt$Freq == 0, NA, 1/freq_dt$Freq)

meta <- left_join(combined_sct@meta.data, freq_dt, by = c("Assignment" = "Var1", "Pool" = "Var2"))
rownames(meta) <- rownames(combined_sct@meta.data)

combined_sct_meta <- AddMetaData(combined_sct, meta)


### SCT the dataset
feats <- rownames(combined_sct_meta[["SCT"]]@counts)[which((rowSums(combined_sct_meta[["SCT"]]@counts > 0)/ncol(combined_sct_meta[["SCT"]]@counts)) >= 0.01)]
combined_sct_filt <- subset(combined_sct_meta, features = feats)

saveRDS(combined_sct_filt, paste0(outdir, "time-integrated_filtered_seurat_1pct_expressing.rds"))


DefaultAssay(combined_sct_filt) <- "RNA"

combined_sct_filt <- SCTransform(combined_sct_filt, verbose = TRUE, vars.to.regress = c("scores.G1", "scores.S", "scores.G2M", "percent.mt", "percent.rb"), return.only.var.genes = FALSE)

saveRDS(combined_sct_filt, paste0(outdir, "seurat_joint_SCT_1pct_expressing.rds"))

fwrite(data.table(Gene = rownames(combined_sct_filt)), paste0(outdir, "genes_1pct_expressing.tsv"), sep = "\t")
