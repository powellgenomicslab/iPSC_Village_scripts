#### Run with R from `conda activate baseR402` #####

library(tidyverse)
library(ggplot2)
library(Seurat)
library(scran)
library(SingleCellExperiment)
library(ggpubr)
library(rstatix)
library(ggsci)
library(ggforce)


save_figs <- function(plot, basename, width = 17, height = 17, units = "cm"){
    ggsave(plot, filename = paste0(basename,".png"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".pdf"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".eps"), height = height, width = width, units = units)
}

##### Set up colors #####
cell_line_colors <- c("FSA0006" = "#DDAA33", "MBE1006" = "#BB5566", "TOB0421"= "#004488")
replicate_colors <- c("1" = "#ACD39E", "2" = "#5AAE61", "3" = "#1B7837")
time_colors <- c("Baseline" = "#b9cee4", "Village Day 4" = "#8a92bb", "Thawed Village Day 0" = "#7d57a0", "Thawed Village Day 7" = "#853786")


dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
outdir <- paste0(dir,"/output/QC_metric_figs/")
dir.create(outdir)


##### Read in the data that has been filtered (only filtered on mt %)
seurat <- readRDS(paste0(dir,"output/Seurat_w_Hash_freeze_thaw/CellCycle/seurat_filtered_cell_cycle.rds.rds"))





table(seurat@meta.data$Location_Time,seurat@meta.data$phases)
prop.table(table(seurat@meta.data$Location_Time,seurat@meta.data$phases), margin = 1)
table(paste0(seurat@meta.data$Location_Time,seurat@meta.data$Site_rep) ,seurat@meta.data$phases)
prop.table(table(paste0(seurat@meta.data$Location_Time,seurat@meta.data$Site_rep) ,seurat@meta.data$phases), margin = 1)
table(paste0(seurat@meta.data$Final_Assignment, seurat@meta.data$Location_Time,seurat@meta.data$Site_rep) ,seurat@meta.data$phases)
prop_df <- prop.table(table(paste0(seurat@meta.data$Final_Assignment, ".",seurat@meta.data$Location_Time,".", seurat@meta.data$Site_rep) ,seurat@meta.data$phases), margin = 1)

prop_df_long <- as.data.frame(prop_df)
colnames(prop_df_long) <- c("Category", "CellCycle", "Proportion")
prop_df_long <- separate(prop_df_long, sep = "\\.", col = "Category", into = c("CellLine", "Location_Time", "Site_rep"))
prop_df_long$Location_Time <- gsub("Brisbane_", "Brisbane.", prop_df_long$Location_Time) %>%
                                gsub("Sydney_", "Sydney.", .) %>%
                                gsub("Melbourne_", "Melbourne.", .)
prop_df_long$Location_Time <- factor(prop_df_long$Location_Time, levels = c("Brisbane.Baseline","Brisbane.Village_Day_4", "Melbourne.Baseline", "Melbourne.Village_Day_4", "Sydney.Baseline", "Sydney.Village_Day_4", "Sydney.Thawed_Village_Day_0", "Sydney.Thawed_Village_Day_7"))

prop_df_long <- separate(prop_df_long, col = "Location_Time", sep = "\\.", into = c("Location", "Time"), remove = FALSE)
prop_df_long$CellCycle <- factor(prop_df_long$CellCycle, levels = c("G1","S","G2M"))


##### Make figure of stacked bar plot showing proportion of cells in each cell cycle group #####
cell_cycle_loc_time <- ggbarplot(prop_df_long, "Location_Time", "Proportion", add = c("mean_se"), fill = "CellCycle") +
    rotate_x_text(45) +
    scale_fill_manual(values = c("#4393C3", "#92C5DE", "#D1E5F0"))

save_figs(cell_cycle_loc_time,  paste0(outdir,"cell_cycle_stacked_bar_location_time"))


cell_cycle_loc_time_line <- ggbarplot(prop_df_long, "Location_Time", "Proportion", add = c("mean_se"), fill = "CellCycle", facet.by = c("CellLine")) +
    rotate_x_text(45) +
    scale_fill_manual(values = c("#4393C3", "#92C5DE", "#D1E5F0"))

save_figs(cell_cycle_loc_time_line, paste0(outdir,"cell_cycle_stacked_bar_location_time_line"))


##### Make figure for proportion from each individual #####
### All together ###
indiv_prop_df <- prop.table(table(paste0(seurat@meta.data$Location_Time,".", seurat@meta.data$Site_rep) ,seurat@meta.data$Final_Assignment), margin = 1)

indiv_prop_df_long <- as.data.frame(indiv_prop_df)
colnames(indiv_prop_df_long) <- c("Category", "CellLine", "Proportion")
indiv_prop_df_long <- separate(indiv_prop_df_long, sep = "\\.", col = "Category", into = c("Location_Time", "Site_rep"))
indiv_prop_df_long$Location_Time <- gsub("Brisbane_", "Brisbane.", indiv_prop_df_long$Location_Time) %>%
                                gsub("Sydney_", "Sydney.", .) %>%
                                gsub("Melbourne_", "Melbourne.", .)
indiv_prop_df_long$Location_Time <- factor(indiv_prop_df_long$Location_Time, levels = c("Brisbane.Baseline","Brisbane.Village_Day_4", "Melbourne.Baseline", "Melbourne.Village_Day_4", "Sydney.Baseline", "Sydney.Village_Day_4", "Sydney.Thawed_Village_Day_0", "Sydney.Thawed_Village_Day_7"))

indiv_prop_df_long <- separate(indiv_prop_df_long, col = "Location_Time", sep = "\\.", into = c("Location", "Time"), remove = FALSE)


cell_line_loc_time <- ggbarplot(indiv_prop_df_long, "Location_Time", "Proportion", add = c("mean_se"), fill = "CellLine") +
    rotate_x_text(45) +
    scale_fill_manual(values = cell_line_colors)

save_figs(cell_line_loc_time,  paste0(outdir,"cell_line_stacked_bar_location_time"))


##### Make normal QC metric figures #####

### Prepare Seurat object ###
seurat@meta.data$Location_Time <- factor(seurat@meta.data$Location_Time, levels = c("Brisbane_Baseline","Brisbane_Village_Day_4", "Melbourne_Baseline", "Melbourne_Village_Day_4", "Sydney_Baseline", "Sydney_Village_Day_4", "Sydney_Thawed_Village_Day_0", "Sydney_Thawed_Village_Day_7"))
seurat@meta.data$Location <- gsub("\\d","", seurat@meta.data$Site_rep)

seurat@meta.data$Time_Location <- paste0(seurat@meta.data$Time,"_",seurat@meta.data$Location)
seurat@meta.data$Time_Location <- gsub(" ", "_", seurat@meta.data$Time_Location)
seurat@meta.data$Time_Location <- factor(seurat@meta.data$Time_Location, levels = c("Baseline_Brisbane","Baseline_Melbourne","Baseline_Sydney", "Village_Day_4_Brisbane",  "Village_Day_4_Melbourne", "Village_Day_4_Sydney", "Thawed_Village_Day_0_Sydney", "Thawed_Village_Day_7_Sydney"))

seurat@meta.data$Time <- factor(seurat@meta.data$Time , levels = c("Baseline", "Village Day 4", "Thawed Village Day 0", "Thawed Village Day 7"))


### Mitochondrial percent ###
## By Pool ##
mt_pct <- ggplot(seurat@meta.data, aes(Pool, percent.mt, fill = Time)) +
            geom_boxplot() +
            theme_classic() +
            scale_fill_manual(values = time_colors) +
            theme(axis.text.x = element_text(angle = 45,  hjust=1)) +
            ylab("Mitochondrial Percent")


save_figs(mt_pct,  paste0(outdir,"mt_percent_boxplot"))

mt_pct_vio <- ggplot(seurat@meta.data, aes(Pool, percent.mt, fill = Time)) +
            geom_violin() +
            theme_classic() +
            scale_fill_manual(values = time_colors)+
            theme(axis.text.x = element_text(angle = 45,  hjust=1)) +
            ylab("Mitochondrial Percent")


save_figs(mt_pct_vio,  paste0(outdir,"mt_percent_violin"))


## By location and time ##
mt_pct_loc_time <- ggplot(seurat@meta.data, aes(Time_Location, percent.mt, fill = gsub("\\D", "", seurat@meta.data$Site_rep) )) +
            geom_boxplot(outlier.size = 0.5) +
            theme_classic() +
            scale_fill_manual(values = replicate_colors) +
            labs(fill =  "Replicate") +
            theme(axis.text.x = element_text(angle = 45,  hjust=1)) +
            ylab("Mitochondrial Percent")


save_figs(mt_pct_loc_time,  paste0(outdir,"mt_percent_location_time_boxplot"), width = 24)

mt_pct_loc_time_vio <- ggplot(seurat@meta.data, aes(Time_Location, percent.mt, fill  = gsub("\\D", "", seurat@meta.data$Site_rep))) +
            geom_violin() +
            theme_classic() +
            scale_fill_manual(values = replicate_colors) +
            labs(fill =  "Replicate")+
            theme(axis.text.x = element_text(angle = 45,  hjust=1)) +
            ylab("Mitochondrial Percent")


save_figs(mt_pct_loc_time_vio,  paste0(outdir,"mt_percent_location_time_violin"), width = 24)


### N genes ###
## By Pool ##
n_genes <- ggplot(seurat@meta.data, aes(Pool, nFeature_RNA, fill = Time)) +
            geom_boxplot() +
            theme_classic() +
            scale_fill_manual(values = time_colors)+
            theme(axis.text.x = element_text(angle = 45,  hjust=1)) +
            ylab("Number Genes")


save_figs(n_genes,  paste0(outdir,"Ngenes_boxplot"))

n_genes_vio <- ggplot(seurat@meta.data, aes(Pool, nFeature_RNA, fill = Time)) +
            geom_violin() +
            theme_classic() +
            scale_fill_manual(values = time_colors)+
            theme(axis.text.x = element_text(angle = 45,  hjust=1)) +
            ylab("Number Genes")


save_figs(n_genes_vio,  paste0(outdir,"Ngenes_violin"))


## By location and time ##
n_genes_loc_time <- ggplot(seurat@meta.data, aes(Time_Location, nFeature_RNA, fill = gsub("\\D", "", seurat@meta.data$Site_rep) )) +
            geom_boxplot(outlier.size = 0.5) +
            theme_classic() +
            scale_fill_manual(values = replicate_colors) +
            labs(fill =  "Replicate")+
            theme(axis.text.x = element_text(angle = 45,  hjust=1)) +
            ylab("Number Genes")

save_figs(n_genes_loc_time,  paste0(outdir,"Ngenes_location_time_boxplot"), width = 24)

n_genes_loc_time_vio <- ggplot(seurat@meta.data, aes(Time_Location, nFeature_RNA, fill  = gsub("\\D", "", seurat@meta.data$Site_rep))) +
            geom_violin() +
            theme_classic() +
            scale_fill_manual(values = replicate_colors) +
            labs(fill =  "Replicate")+
            theme(axis.text.x = element_text(angle = 45,  hjust=1)) +
            ylab("Number Genes")

save_figs(n_genes_loc_time_vio,  paste0(outdir,"Ngenes_location_time_violin"), width = 24)



### N UMIs ###
n_UMIs <- ggplot(seurat@meta.data, aes(Pool, nCount_RNA, fill = Time)) +
            geom_boxplot() +
            theme_classic() +
            scale_fill_manual(values = time_colors)+
            theme(axis.text.x = element_text(angle = 45,  hjust=1)) +
            ylab("Number UMIs")

save_figs(n_genes,  paste0(outdir,"n_UMIs_RNA_boxplot"))

nCount_RNA_vio <- ggplot(seurat@meta.data, aes(Pool, nCount_RNA, fill = Time)) +
            geom_violin() +
            theme_classic() +
            scale_fill_manual(values = time_colors)+
            theme(axis.text.x = element_text(angle = 45,  hjust=1)) +
            ylab("Number UMIs")

save_figs(n_genes_vio,  paste0(outdir,"n_UMIs_RNA_violin"))


## By location and time ##
nCount_RNA_loc_time <- ggplot(seurat@meta.data, aes(Time_Location, nCount_RNA, fill = gsub("\\D", "", seurat@meta.data$Site_rep) )) +
            geom_boxplot(outlier.size = 0.5) +
            theme_classic() +
            scale_fill_manual(values = replicate_colors) +
            labs(fill =  "Replicate")+
            theme(axis.text.x = element_text(angle = 45,  hjust=1)) +
            ylab("Number UMIs")

save_figs(n_genes_loc_time,  paste0(outdir,"n_UMIs_RNA_location_time_boxplot"), width = 24)

nCount_RNA_loc_time_vio <- ggplot(seurat@meta.data, aes(Time_Location, nCount_RNA, fill  = gsub("\\D", "", seurat@meta.data$Site_rep))) +
            geom_violin() +
            theme_classic() +
            scale_fill_manual(values = replicate_colors) +
            labs(fill =  "Replicate")+
            theme(axis.text.x = element_text(angle = 45,  hjust=1)) +
            ylab("Number UMIs")

save_figs(n_genes_loc_time_vio,  paste0(outdir,"n_UMIs_RNA_location_time_violin"), width = 24)
