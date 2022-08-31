##### Read in Arguments #####
print("Reading and assigning input arguments")


dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
outdir <- paste0(dir,"output/All_Integrated_QC/")
datadir <- paste0(dir,"output/All_data_integrated_remove_bad/")

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
library(ggpubr)
library(ggforce)


##### Make function for saving figures in different formats #####
save_figs <- function(plot, basename, width = 17, height = 17, units = "cm"){
    ggsave(plot, filename = paste0(basename,".png"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".pdf"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".eps"), height = height, width = width, units = units)
}


##### Read in Data #####
seurat <- readRDS(paste0(datadir,"seurat_integrated_all_times_clustered.rds"))
seurat@meta.data$Location <- gsub("_.+", "", seurat@meta.data$Location_Time)
seurat@meta.data$phases <- factor(seurat@meta.data$phases, levels = c("G1", "S", "G2M"))
seurat@meta.data$Village <- gsub("Baseline", "Uni-Culture", seurat@meta.data$Time) %>% gsub("Village Day 4", "Village", .) %>% gsub("Thawed Village Day 0", "Uni-Culture", .) %>% gsub("Thawed Village Day 7", "Village", .)

##### Set up Colors #####
cell_line_colors <- c("FSA0006" = "#F79E29", "MBE1006" = "#9B2C99", "TOB0421"= "#35369C")
replicate_colors <- c("1" = "#ACD39E", "2" = "#5AAE61", "3" = "#1B7837")
time_colors <- c("Baseline" = "#b9cee4", "Village Day 4" = "#8a92bb", "Thawed Village Day 0" = "#7d57a0", "Thawed Village Day 7" = "#853786")
village_colors <- c("Uni-Culture" = "#613246", "Village" = "#A286AA")
cycle_colors <- c("G1" = "#4393C3", "S" = "#92C5DE", "G2M" = "#D1E5F0")
site_colors <- c("Brisbane" = "#536CB4", "Sydney" = "#62BD67", "Melbourne" = "#D24F72")
cluster_colors <- c("0" = "#C9D8EA", "1" = "#928CC5", "2" = "#A7C9A9", "3" = "#179085", "4" = "#F79F9F", "5" = "#C35B76", "6" = "#F4C893", "7" = "#F6AA4B")


##### Make UMAPS #####
## With Clusters ##
UMAP_cell_line <- DimPlot(seurat, reduction = "umap", group.by = c("integrated_snn_res.0.28"), cols = cluster_colors) + labs(color="Cluster") + ggtitle(NULL) 
save_figs(UMAP_cell_line, basename = paste0(outdir,"Cluster_umap"),width = 15, height = 15)

## With Village ##
UMAP_village <- DimPlot(seurat, reduction = "umap", group.by = c("Village"), cols = village_colors) + labs(color="Village") + ggtitle(NULL) 
save_figs(UMAP_village, basename = paste0(outdir,"Village_umap"),width = 15, height = 15)

## By Cell Line ##
UMAP_cell_line <- DimPlot(seurat, reduction = "umap", group.by = c("Final_Assignment"), cols = cell_line_colors) + labs(color="Cell Line") + ggtitle(NULL) 
save_figs(UMAP_cell_line, basename = paste0(outdir,"Cell_Line_umap"),width = 15, height = 15)

## By Location ##
UMAP_site <- DimPlot(seurat, reduction = "umap", group.by = c("Location"), cols = alpha(site_colors, 0.7)) + labs(color="Location") + ggtitle(NULL)
UMAP_site[[1]]$layers[[1]]$aes_params$shape = 16
UMAP_site[[1]]$layers[[1]]$aes_params$size = 0.3
save_figs(UMAP_site, basename = paste0(outdir,"Location_umap"),width = 15, height = 15)

## By Cell Cycle Phase ##
Idents(seurat) <- "phases"
UMAP_cell_cycle <- DimPlot(seurat, reduction = "umap", group.by = c("phases"), cols = cycle_colors, order = c("G2M", "S", "G1")) + labs(color="Cell Cycle\nPhase") + ggtitle(NULL)
UMAP_cell_cycle[[1]]$layers[[1]]$aes_params$shape = 16
UMAP_cell_cycle[[1]]$layers[[1]]$aes_params$size = 0.3
save_figs(UMAP_cell_cycle, basename = paste0(outdir,"CellCycle_umap"),width = 15, height = 15)

## With mt % ##
Idents(seurat) <- "phases"
UMAP_mt_pct <- FeaturePlot(seurat, reduction = "umap", features = c("percent.mt")) + labs(color="Mitochondrial\nPercent") + ggtitle(NULL) + scale_colour_gradient(low = "grey97", high = "red4")
save_figs(UMAP_mt_pct, basename = paste0(outdir,"MTpercent_umap"),width = 16, height = 15)



##### Make Bar Charts and Box Plots #####
## Cluster per cell line ##
cluster_per_line <- ggplot(seurat@meta.data, aes(integrated_snn_res.0.28, fill = Final_Assignment)) +
	geom_bar(position = "fill") +
	theme_classic() +
	scale_fill_manual(values = cell_line_colors) +
	labs(fill="Cell Line")+
	xlab("Cluster")
save_figs(cluster_per_line, basename = paste0(outdir,"cluster_per_line"))


## Cell line per cluster ##
line_per_cluster <- ggplot(seurat@meta.data, aes(Final_Assignment, fill = integrated_snn_res.0.28)) +
	geom_bar(position = "fill") +
	theme_classic() +
	scale_fill_manual(values = cluster_colors) +
	labs(fill="Cluster")+
	xlab("Cell Line")
save_figs(line_per_cluster, basename = paste0(outdir,"line_per_cluster"), width = 10)


## Cluster per Cell Cycle ##
cluster_per_cell_cycle <- ggplot(seurat@meta.data, aes(integrated_snn_res.0.28, fill = phases)) +
	geom_bar(position = "fill") +
	theme_classic() +
	scale_fill_manual(values = cycle_colors) +
	labs(fill="Cell Cycle") +
	xlab("Cluster")
save_figs(cluster_per_cell_cycle, basename = paste0(outdir,"cluster_per_cell_cycle"))


## Cell Cycle per cluster ##
cell_cycle_per_cluster <- ggplot(seurat@meta.data, aes(phases, fill = integrated_snn_res.0.28)) +
	geom_bar(position = "fill") +
	theme_classic() +
	scale_fill_manual(values = cluster_colors) +
	labs(fill="Cluster") +
	xlab("Cell Cycle")
save_figs(cell_cycle_per_cluster, basename = paste0(outdir,"cell_cycle_per_cluster"), width = 10)


## Cluster per Location ##
cluster_per_cell_cycle <- ggplot(seurat@meta.data, aes(integrated_snn_res.0.28, fill = Location)) +
	geom_bar(position = "fill") +
	theme_classic() +
	scale_fill_manual(values = site_colors) +
	labs(fill="Cell Cycle") +
	xlab("Location")
save_figs(cluster_per_cell_cycle, basename = paste0(outdir,"cluster_per_location"))


## Location per cluster ##
cell_cycle_per_cluster <- ggplot(seurat@meta.data, aes(Location, fill = integrated_snn_res.0.28)) +
	geom_bar(position = "fill") +
	theme_classic() +
	scale_fill_manual(values = cluster_colors) +
	labs(fill="Cluster") +
	xlab("Location")
save_figs(cell_cycle_per_cluster, basename = paste0(outdir,"location_per_cluster"), width = 10)


## Cluster per Mt % ##
mt_pct_per_cluster <- ggplot(seurat@meta.data, aes(integrated_snn_res.0.28, percent.mt, fill = integrated_snn_res.0.28)) +
	geom_boxplot() +
	theme_classic() +
	scale_fill_manual(values = cluster_colors) +
	labs(fill="Cluster") +
	ylab("Mitochondrial Percent") +
	xlab("Cluster")
save_figs(mt_pct_per_cluster, basename = paste0(outdir,"mt_pct_cluster_per"))





##### Bar Plot with Clusters by location, time, facet line #####
table(seurat@meta.data$Location_Time,seurat@meta.data$integrated_snn_res.0.28)
prop.table(table(seurat@meta.data$Location_Time,seurat@meta.data$integrated_snn_res.0.28), margin = 1)
table(paste0(seurat@meta.data$Location_Time,seurat@meta.data$Site_rep) ,seurat@meta.data$integrated_snn_res.0.28)
prop.table(table(paste0(seurat@meta.data$Location_Time,seurat@meta.data$Site_rep) ,seurat@meta.data$integrated_snn_res.0.28), margin = 1)
table(paste0(seurat@meta.data$Final_Assignment, seurat@meta.data$Location_Time,seurat@meta.data$Site_rep) ,seurat@meta.data$phases)
prop_df <- prop.table(table(paste0(seurat@meta.data$Final_Assignment, ".",seurat@meta.data$Location_Time,".", seurat@meta.data$Site_rep) ,seurat@meta.data$integrated_snn_res.0.28), margin = 1)

prop_df_long <- as.data.frame(prop_df)
colnames(prop_df_long) <- c("Category", "Cluster", "Proportion")
prop_df_long <- separate(prop_df_long, sep = "\\.", col = "Category", into = c("CellLine", "Location_Time", "Site_rep"))
prop_df_long$Location_Time <- gsub("Brisbane_", "Brisbane.", prop_df_long$Location_Time) %>%
                                gsub("Sydney_", "Sydney.", .) %>%
                                gsub("Melbourne_", "Melbourne.", .)
prop_df_long$Location_Time <- factor(prop_df_long$Location_Time, levels = c("Brisbane.Baseline","Brisbane.Village_Day_4", "Melbourne.Baseline", "Melbourne.Village_Day_4", "Sydney.Baseline", "Sydney.Village_Day_4", "Sydney.Thawed_Village_Day_0", "Sydney.Thawed_Village_Day_7"))

prop_df_long <- separate(prop_df_long, col = "Location_Time", sep = "\\.", into = c("Location", "Time"), remove = FALSE)
prop_df_long$Cluster <- factor(prop_df_long$Cluster, levels = c(0,1,2,3,4,5,6,7))


cluster_loc_time_line <- ggbarplot(prop_df_long, "Location_Time", "Proportion", add = c("mean_se"), fill = "Cluster", facet.by = c("CellLine")) +
    rotate_x_text(45) +
    scale_fill_manual(values = cluster_colors)

save_figs(cluster_loc_time_line, paste0(outdir,"cell_cycle_stacked_bar_location_time_line"))
