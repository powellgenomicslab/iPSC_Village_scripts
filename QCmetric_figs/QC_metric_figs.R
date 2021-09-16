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
library(data.table)


save_figs <- function(plot, basename, width = 17, height = 17, units = "cm"){
    ggsave(plot, filename = paste0(basename,".png"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".pdf"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".eps"), height = height, width = width, units = units)
}

##### Set up colors #####
cell_line_colors <- c("FSA0006" = "#F79E29", "MBE1006" = "#9B2C99", "TOB0421"= "#35369C")
replicate_colors <- c("1" = "#ACD39E", "2" = "#5AAE61", "3" = "#1B7837")
time_colors <- c("Baseline" = "#b9cee4", "Village Day 4" = "#8a92bb", "Thawed Village Day 0" = "#7d57a0", "Thawed Village Day 7" = "#853786")

site_updates <- c("Brisbane" = "Site 1", "Sydney" = "Site 3" ,"Melbourne" = "Site 2")


dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
outdir <- paste0(dir,"/output/QC_metric_figs/")
dir.create(outdir)


##### Read in the data that has been filtered (only filtered on mt %)
seurat <- readRDS(paste0(dir,"output/Seurat_w_Hash_freeze_thaw/CellCycle/seurat_filtered_cell_cycle.rds"))





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


for (location in names(site_updates)){
	prop_df_long$Location <- gsub(location, site_updates[location], prop_df_long$Location)
}
prop_df_long$Village <- ifelse((prop_df_long$Time == "Thawed_Village_Day_0" | prop_df_long$Time == "Baseline"), "Baseline", "Village")
prop_df_long$Cryopreservation <- ifelse((prop_df_long$Time == "Thawed_Village_Day_0" | prop_df_long$Time == "Thawed_Village_Day_7"), "Cryopreserved", "Fresh")
prop_df_long$Cryopreservation <- factor(prop_df_long$Cryopreservation, levels = c("Fresh", "Cryopreserved"))
prop_df_long$Village_Rep <- paste0(prop_df_long$Village, "-", gsub("\\D+","Replicate ", prop_df_long$Site_rep))
prop_df_long$Village_Line <- paste0(prop_df_long$Village, "-",  prop_df_long$CellLine)
prop_df_long$Village_Line <- factor(prop_df_long$Village_Line, levels = c("Baseline-FSA0006", "Baseline-MBE1006", "Baseline-TOB0421", "Village-FSA0006", "Village-MBE1006", "Village-TOB0421"))


##### Make figure of stacked bar plot showing proportion of cells in each cell cycle group #####
# cell_cycle_loc_time <- ggbarplot(prop_df_long, "Location_Time", "Proportion", add = c("mean_se"), fill = "CellCycle") +
#     rotate_x_text(45) +
#     scale_fill_manual(values = c("#4393C3", "#92C5DE", "#D1E5F0"))

# save_figs(cell_cycle_loc_time,  paste0(outdir,"cell_cycle_stacked_bar_location_time"))


# cell_cycle_loc_time_line <- ggbarplot(prop_df_long, "Location_Time", "Proportion", add = c("mean_se"), fill = "CellCycle", facet.by = c("CellLine")) +
#     rotate_x_text(45) +
#     scale_fill_manual(values = c("#4393C3", "#92C5DE", "#D1E5F0"))

# save_figs(cell_cycle_loc_time_line, paste0(outdir,"cell_cycle_stacked_bar_location_time_line"))


## Fresh samples - facet by site and plot each replicate at baseline and village on x axis for each
cell_cycle_fresh <- ggbarplot(prop_df_long[which(prop_df_long$Cryopreservation == "Fresh"),], "Village_Rep", "Proportion", add = c("mean_se"), fill = "CellCycle", facet.by = c("Location"), scales = "free_x", legend = "right") +
    rotate_x_text(45) +
    scale_fill_manual(values = c("#4393C3", "#92C5DE", "#D1E5F0")) +
	theme(axis.title.x=element_blank())

save_figs(cell_cycle_fresh,  paste0(outdir,"cell_cycle_fresh_facet"), width = 13, height = 8.5)


## Fresh samples - facet by site and plot each replicate at baseline and village on x axis for each
cell_cycle_fresh_line <- ggbarplot(prop_df_long[which(prop_df_long$Cryopreservation == "Fresh"),], "Village_Line", "Proportion", add = c("mean_se"), fill = "CellCycle", facet.by = c("Location"), scales = "free_x", legend = "right") +
    rotate_x_text(45) +
    scale_fill_manual(values = c("#4393C3", "#92C5DE", "#D1E5F0")) +
	theme(axis.title.x=element_blank())

save_figs(cell_cycle_fresh_line,  paste0(outdir,"cell_cycle_fresh_facet_by_line"), width = 13, height = 8.5)




## Cryopreserved samples - facet by site and plot each replicate at baseline and village on x axis for each
cell_cycle_cryo <- ggbarplot(prop_df_long[which(prop_df_long$Location == "Site 3"),], "Village_Rep", "Proportion", add = c("mean_se"), fill = "CellCycle", facet.by = c("Cryopreservation"), scales = "free_x", legend = "right") +
    rotate_x_text(45) +
    scale_fill_manual(values = c("#4393C3", "#92C5DE", "#D1E5F0")) +
	theme(axis.title.x=element_blank())

save_figs(cell_cycle_cryo,  paste0(outdir,"cell_cycle_cryo_facet"), width = 10, height = 8.5)




## Cryopreserved samples - facet by site and plot each replicate at baseline and village on x axis for each
cell_cycle_cryo <- ggbarplot(prop_df_long[which(prop_df_long$Location == "Site 3"),], "Village_Line", "Proportion", add = c("mean_se"), fill = "CellCycle", facet.by = c("Cryopreservation"), scales = "free_x", legend = "right") +
    rotate_x_text(45) +
    scale_fill_manual(values = c("#4393C3", "#92C5DE", "#D1E5F0")) +
	theme(axis.title.x=element_blank())

save_figs(cell_cycle_cryo,  paste0(outdir,"cell_cycle_cryo_facet_by_line"), width = 10, height = 8.5)





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

### Make column for location only and update to sites 1, 2 and 3
indiv_prop_df_long <- separate(indiv_prop_df_long, col = "Location_Time", sep = "\\.", into = c("Location", "Time"), remove = FALSE)

for (location in names(site_updates)){
	indiv_prop_df_long$Location <- gsub(location, site_updates[location], indiv_prop_df_long$Location)
}


cell_line_loc_time <- ggbarplot(indiv_prop_df_long, "Location_Time", "Proportion", add = c("mean_se"), fill = "CellLine") +
    rotate_x_text(45) +
    scale_fill_manual(values = cell_line_colors)

save_figs(cell_line_loc_time,  paste0(outdir,"cell_line_stacked_bar_location_time"))



## Proportion of cell lines for just first two timepoints (baseline and village day 4)
indiv_prop_df_long$Time <- gsub("Village_Day_4", "Village", indiv_prop_df_long$Time)
cell_line_loc_time_first2 <- ggbarplot(indiv_prop_df_long[which(indiv_prop_df_long$Time == "Baseline" | indiv_prop_df_long$Time == "Village"),], "Time", "Proportion", add = c("mean_se"), fill = "CellLine", facet.by = "Location", size = 0.25) +
    scale_fill_manual(values = cell_line_colors) +
	rotate_x_text(45) +
	theme(axis.title.x=element_blank())

save_figs(cell_line_loc_time_first2,  paste0(outdir,"cell_line_stacked_bar_location_time_baseline_village_4days"), width = 6, height = 8)


### Compare the FSA0006 level at baseline vs village ###
indiv_prop_df_long <- data.table(indiv_prop_df_long)

t_test_list <- list()
t_test_df <- data.table(unique(indiv_prop_df_long[, c("Location", "CellLine")]), "Paired t-test t" = 0, "Paired t-test P" = 0)
colnames(t_test_df)[1:2] <- c("Site", "Cell Line")

for (site in unique(indiv_prop_df_long$Location)){
	for (line in unique(indiv_prop_df_long$CellLine)){
		t_test_list[[paste0(site, "_", line)]] <- t.test(indiv_prop_df_long[Time == "Baseline" & Location == site & CellLine == line]$Proportion, indiv_prop_df_long[Time == "Village" & Location == site & CellLine == line]$Proportion, paired = TRUE)
		t_test_df[Site == site & `Cell Line` == line, "Paired t-test t"] <- t_test_list[[paste0(site, "_", line)]]$statistic
		t_test_df[Site == site & `Cell Line` == line]$`Paired t-test P` <- t_test_list[[paste0(site, "_", line)]]$p.value
	}
}

fwrite(t_test_df, paste0(outdir,"hPSC_line_proportions.tsv"), sep = "\t")


## Proportion of cell lines for Sydney
indiv_prop_df_long_syd <- indiv_prop_df_long[which(indiv_prop_df_long$Location == "Site 3"),]
indiv_prop_df_long_syd$Village <- ifelse((indiv_prop_df_long_syd$Time == "Thawed_Village_Day_0" | indiv_prop_df_long_syd$Time == "Baseline"), "Baseline", "Village")
indiv_prop_df_long_syd$Cryopreservation <- ifelse((indiv_prop_df_long_syd$Time == "Thawed_Village_Day_0" | indiv_prop_df_long_syd$Time == "Thawed_Village_Day_7"), "Cryo-\npreserved", "Fresh")
indiv_prop_df_long_syd$Cryopreservation <- factor(indiv_prop_df_long_syd$Cryopreservation, levels = c("Fresh", "Cryo-\npreserved"))

cell_line_loc_time_first3 <- ggbarplot(indiv_prop_df_long_syd, "Village", "Proportion", add = c("mean_se"), fill = "CellLine", facet.by = "Cryopreservation", scales = "free_x") +
    scale_fill_manual(values = cell_line_colors) +
	rotate_x_text(45) +
	theme(axis.title.x = element_blank())

save_figs(cell_line_loc_time_first3,  paste0(outdir,"cell_line_stacked_bar_location_time_Sydney"), width = 6, height = 8)



t_test_cryo_list <- list()
t_test_cryo_df <- data.table(unique(indiv_prop_df_long_syd[, c("Cryopreservation", "CellLine")]), "Paired t-test t" = 0, "Paired t-test P" = 0)
colnames(t_test_cryo_df)[1:2] <- c("Site", "Cell Line")

for (site in unique(indiv_prop_df_long_syd$Cryopreservation)){
	for (line in unique(indiv_prop_df_long_syd$CellLine)){
		t_test_cryo_list[[paste0(site, "_", line)]] <- t.test(indiv_prop_df_long_syd[Village == "Baseline" & Cryopreservation == site & CellLine == line]$Proportion, indiv_prop_df_long_syd[Village == "Village" & Cryopreservation == site & CellLine == line]$Proportion, paired = TRUE)
		t_test_cryo_df[Site == site & `Cell Line` == line, "Paired t-test t"] <- t_test_cryo_list[[paste0(site, "_", line)]]$statistic
		t_test_cryo_df[Site == site & `Cell Line` == line]$`Paired t-test P` <- t_test_cryo_list[[paste0(site, "_", line)]]$p.value
	}
}

t_test_cryo_df$Site <- gsub("-\n", "", t_test_cryo_df$Site)

fwrite(t_test_cryo_df, paste0(outdir,"hPSC_line_proportions_cryoperserved.tsv"), sep = "\t")


### Compare Fresh vs Cryo baseline ###
t_test_cryo_base_list <- list()
t_test_cryo_base_df <- data.table(unique(indiv_prop_df_long_syd$CellLine), "Paired t-test t" = 0, "Paired t-test P" = 0)
colnames(t_test_cryo_base_df)[1] <- c("Cell Line")


for (line in unique(indiv_prop_df_long_syd$CellLine)){
	t_test_cryo_base_list[[line]] <- t.test(indiv_prop_df_long_syd[Village == "Baseline" & Cryopreservation == "Fresh" & CellLine == line]$Proportion, indiv_prop_df_long_syd[Village == "Baseline" & Cryopreservation == "Cryo-\npreserved" & CellLine == line]$Proportion, paired = TRUE)
	t_test_cryo_base_df[`Cell Line` == line, "Paired t-test t"] <- t_test_cryo_base_list[[line]]$statistic
	t_test_cryo_base_df[`Cell Line` == line]$`Paired t-test P` <- t_test_cryo_base_list[[line]]$p.value
}





##### Make normal QC metric figures #####
### Prepare Seurat object ###
seurat@meta.data$Location_Time <- factor(seurat@meta.data$Location_Time, levels = c("Brisbane_Baseline","Brisbane_Village_Day_4", "Melbourne_Baseline", "Melbourne_Village_Day_4", "Sydney_Baseline", "Sydney_Village_Day_4", "Sydney_Thawed_Village_Day_0", "Sydney_Thawed_Village_Day_7"))
seurat@meta.data$Location <- gsub("\\d","", seurat@meta.data$Site_rep)
for (location in names(site_updates)){
	seurat@meta.data$Location <- gsub(location, site_updates[location], seurat@meta.data$Location)
}

seurat@meta.data$Village <- ifelse((seurat@meta.data$Time == "Thawed Village Day 0" | seurat@meta.data$Time == "Baseline"), "Baseline", "Village")
seurat@meta.data$Cryopreservation <- ifelse((seurat@meta.data$Time == "Thawed Village Day 0" | seurat@meta.data$Time == "Thawed Village Day 7"), "Cryopreserved", "Fresh")
seurat@meta.data$Cryopreservation <- factor(seurat@meta.data$Cryopreservation, levels = c("Fresh", "Cryopreserved"))


seurat@meta.data$Time_Location <- paste0(seurat@meta.data$Village," ",seurat@meta.data$Location, " ", seurat@meta.data$Cryopreservation)
seurat@meta.data$Time_Location <- factor(seurat@meta.data$Time_Location, levels = c("Baseline Site 1 Fresh","Baseline Site 2 Fresh","Baseline Site 3 Fresh", "Village Site 1 Fresh",  "Village Site 2 Fresh", "Village Site 3 Fresh", "Baseline Site 2 Cryopreserved", "Village Site 2 Cryopreserved"))

seurat@meta.data$Time <- factor(seurat@meta.data$Time , levels = c("Baseline", "Village Day 4", "Thawed Village Day 0", "Thawed Village Day 7"))


### Mitochondrial percent ###
# ## By Pool ##
# mt_pct <- ggplot(seurat@meta.data, aes(Pool, percent.mt, fill = Time)) +
#             geom_boxplot() +
#             theme_classic() +
#             scale_fill_manual(values = time_colors) +
#             theme(axis.text.x = element_text(angle = 45,  hjust=1)) +
#             ylab("Mitochondrial Percent")


# save_figs(mt_pct,  paste0(outdir,"mt_percent_boxplot"))

# mt_pct_vio <- ggplot(seurat@meta.data, aes(Pool, percent.mt, fill = Time)) +
#             geom_violin() +
#             theme_classic() +
#             scale_fill_manual(values = time_colors)+
#             theme(axis.text.x = element_text(angle = 45,  hjust=1)) +
#             ylab("Mitochondrial Percent")


# save_figs(mt_pct_vio,  paste0(outdir,"mt_percent_violin"))


# ## By location and time ##
# mt_pct_loc_time <- ggplot(seurat@meta.data, aes(Time_Location, percent.mt, fill = gsub("\\D", "", seurat@meta.data$Site_rep) )) +
#             geom_boxplot(outlier.size = 0.5) +
#             theme_classic() +
#             scale_fill_manual(values = replicate_colors) +
#             labs(fill =  "Replicate") +
#             theme(axis.text.x = element_text(angle = 45,  hjust=1)) +
#             ylab("Mitochondrial Percent")


# save_figs(mt_pct_loc_time,  paste0(outdir,"mt_percent_location_time_boxplot"), width = 24)

# mt_pct_loc_time_vio <- ggplot(seurat@meta.data, aes(Time_Location, percent.mt, fill  = gsub("\\D", "", seurat@meta.data$Site_rep))) +
#             geom_violin() +
#             theme_classic() +
#             scale_fill_manual(values = replicate_colors) +
#             labs(fill =  "Replicate")+
#             theme(axis.text.x = element_text(angle = 45,  hjust=1)) +
#             ylab("Mitochondrial Percent")


# save_figs(mt_pct_loc_time_vio,  paste0(outdir,"mt_percent_location_time_violin"), width = 24)


## Fresh samples - facet by site and plot each replicate at baseline and village on x axis for each
mt_pct_fresh <- ggplot(seurat@meta.data[which(seurat@meta.data$Cryopreservation == "Fresh"),], aes(Village, percent.mt, fill = gsub("\\D", "", seurat@meta.data[which(seurat@meta.data$Cryopreservation == "Fresh"),]$Site_rep) )) +
            geom_boxplot(outlier.size = 0.25) +
            theme_classic() +
            scale_fill_manual(values = replicate_colors) +
            labs(fill =  "Replicate") +
			theme(axis.title.x = element_blank()) +
            ylab("Mitochondrial Percent") +
			facet_wrap(vars(Location), nrow = 1)

save_figs(mt_pct_fresh,  paste0(outdir,"mt_percent_fresh_facet"), width = 12.5, height = 6)


## Cryopreserved samples - facet by site and plot each replicate at baseline and village on x axis for each
mt_pct_cryo <- ggplot(seurat@meta.data[which(seurat@meta.data$Location == "Site 2"),], aes(Village, percent.mt, fill = gsub("\\D", "", seurat@meta.data[which(seurat@meta.data$Location == "Site 2"),]$Site_rep) )) +
            geom_boxplot(outlier.size = 0.25) +
            theme_classic() +
            scale_fill_manual(values = replicate_colors) +
            labs(fill =  "Replicate") +
			theme(axis.title.x = element_blank()) +
            ylab("Mitochondrial Percent") +
			facet_wrap(vars(Cryopreservation), nrow = 1)

save_figs(mt_pct_cryo,  paste0(outdir,"mt_percent_cryo_facet"), width = 9.5, height = 6)



### N genes ###
# ## By Pool ##
# n_genes <- ggplot(seurat@meta.data, aes(Pool, nFeature_RNA, fill = Time)) +
#             geom_boxplot() +
#             theme_classic() +
#             scale_fill_manual(values = time_colors)+
#             theme(axis.text.x = element_text(angle = 45,  hjust=1)) +
#             ylab("Number Genes")


# save_figs(n_genes,  paste0(outdir,"Ngenes_boxplot"))

# n_genes_vio <- ggplot(seurat@meta.data, aes(Pool, nFeature_RNA, fill = Time)) +
#             geom_violin() +
#             theme_classic() +
#             scale_fill_manual(values = time_colors)+
#             theme(axis.text.x = element_text(angle = 45,  hjust=1)) +
#             ylab("Number Genes")


# save_figs(n_genes_vio,  paste0(outdir,"Ngenes_violin"))


# ## By location and time ##
# n_genes_loc_time <- ggplot(seurat@meta.data, aes(Time_Location, nFeature_RNA, fill = gsub("\\D", "", seurat@meta.data$Site_rep) )) +
#             geom_boxplot(outlier.size = 0.5) +
#             theme_classic() +
#             scale_fill_manual(values = replicate_colors) +
#             labs(fill =  "Replicate")+
#             theme(axis.text.x = element_text(angle = 45,  hjust=1)) +
#             ylab("Number Genes")

# save_figs(n_genes_loc_time,  paste0(outdir,"Ngenes_location_time_boxplot"), width = 24)

# n_genes_loc_time_vio <- ggplot(seurat@meta.data, aes(Time_Location, nFeature_RNA, fill  = gsub("\\D", "", seurat@meta.data$Site_rep))) +
#             geom_violin() +
#             theme_classic() +
#             scale_fill_manual(values = replicate_colors) +
#             labs(fill =  "Replicate")+
#             theme(axis.text.x = element_text(angle = 45,  hjust=1)) +
#             ylab("Number Genes")

# save_figs(n_genes_loc_time_vio,  paste0(outdir,"Ngenes_location_time_violin"), width = 24)



## Fresh samples - facet by site and plot each replicate at baseline and village on x axis for each
n_genes_fresh <- ggplot(seurat@meta.data[which(seurat@meta.data$Cryopreservation == "Fresh"),], aes(Village, nFeature_RNA, fill = gsub("\\D", "", seurat@meta.data[which(seurat@meta.data$Cryopreservation == "Fresh"),]$Site_rep) )) +
            geom_boxplot(outlier.size = 0.25) +
            theme_classic() +
            scale_fill_manual(values = replicate_colors) +
            labs(fill =  "Replicate") +
			theme(axis.title.x = element_blank()) +
            ylab("Number Genes") +
			facet_wrap(vars(Location), nrow = 1)

save_figs(n_genes_fresh,  paste0(outdir,"n_genes_fresh_facet"), width = 13, height = 6)


## Cryopreserved samples - facet by site and plot each replicate at baseline and village on x axis for each
n_genes_cryo <- ggplot(seurat@meta.data[which(seurat@meta.data$Location == "Site 2"),], aes(Village, nFeature_RNA, fill = gsub("\\D", "", seurat@meta.data[which(seurat@meta.data$Location == "Site 2"),]$Site_rep) )) +
            geom_boxplot(outlier.size = 0.25) +
            theme_classic() +
            scale_fill_manual(values = replicate_colors) +
            labs(fill =  "Replicate") +
			theme(axis.title.x = element_blank()) +
            ylab("Number Genes") +
			facet_wrap(vars(Cryopreservation), nrow = 1)

save_figs(n_genes_cryo,  paste0(outdir,"n_genes_cryo_facet"), width = 10, height = 6)



### N UMIs ###
# n_UMIs <- ggplot(seurat@meta.data, aes(Pool, nCount_RNA, fill = Time)) +
#             geom_boxplot() +
#             theme_classic() +
#             scale_fill_manual(values = time_colors)+
#             theme(axis.text.x = element_text(angle = 45,  hjust=1)) +
#             ylab("Number UMIs")

# save_figs(n_genes,  paste0(outdir,"n_UMIs_RNA_boxplot"))

# nCount_RNA_vio <- ggplot(seurat@meta.data, aes(Pool, nCount_RNA, fill = Time)) +
#             geom_violin() +
#             theme_classic() +
#             scale_fill_manual(values = time_colors)+
#             theme(axis.text.x = element_text(angle = 45,  hjust=1)) +
#             ylab("Number UMIs")

# save_figs(n_genes_vio,  paste0(outdir,"n_UMIs_RNA_violin"))


# ## By location and time ##
# nCount_RNA_loc_time <- ggplot(seurat@meta.data, aes(Time_Location, nCount_RNA, fill = gsub("\\D", "", seurat@meta.data$Site_rep) )) +
#             geom_boxplot(outlier.size = 0.5) +
#             theme_classic() +
#             scale_fill_manual(values = replicate_colors) +
#             labs(fill =  "Replicate")+
#             theme(axis.text.x = element_text(angle = 45,  hjust=1)) +
#             ylab("Number UMIs")

# save_figs(n_genes_loc_time,  paste0(outdir,"n_UMIs_RNA_location_time_boxplot"), width = 24)

# nCount_RNA_loc_time_vio <- ggplot(seurat@meta.data, aes(Time_Location, nCount_RNA, fill  = gsub("\\D", "", seurat@meta.data$Site_rep))) +
#             geom_violin() +
#             theme_classic() +
#             scale_fill_manual(values = replicate_colors) +
#             labs(fill =  "Replicate")+
#             theme(axis.text.x = element_text(angle = 45,  hjust=1)) +
#             ylab("Number UMIs")

# save_figs(n_genes_loc_time_vio,  paste0(outdir,"n_UMIs_RNA_location_time_violin"), width = 24)


n_umi_fresh <- ggplot(seurat@meta.data[which(seurat@meta.data$Cryopreservation == "Fresh"),], aes(Village, nCount_RNA, fill = gsub("\\D", "", seurat@meta.data[which(seurat@meta.data$Cryopreservation == "Fresh"),]$Site_rep) )) +
            geom_boxplot(outlier.size = 0.25) +
            theme_classic() +
            scale_fill_manual(values = replicate_colors) +
            labs(fill =  "Replicate") +
			theme(axis.title.x = element_blank()) +
            ylab("Number UMIs") +
			facet_wrap(vars(Location), nrow = 1)

save_figs(n_umi_fresh,  paste0(outdir,"n_umis_fresh_facet"), width = 13.5, height = 6)


## Cryopreserved samples - facet by site and plot each replicate at baseline and village on x axis for each
n_umi_cryo <- ggplot(seurat@meta.data[which(seurat@meta.data$Location == "Site 2"),], aes(Village, nCount_RNA, fill = gsub("\\D", "", seurat@meta.data[which(seurat@meta.data$Location == "Site 2"),]$Site_rep) )) +
            geom_boxplot(outlier.size = 0.25) +
            theme_classic() +
            scale_fill_manual(values = replicate_colors) +
            labs(fill =  "Replicate") +
			theme(axis.title.x = element_blank()) +
            ylab("Number UMIs") +
			facet_wrap(vars(Cryopreservation), nrow = 1)

save_figs(n_umi_cryo,  paste0(outdir,"n_umis_cryo_facet"), width = 10.5, height = 6)
