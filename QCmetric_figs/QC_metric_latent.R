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

site_updates <- c("Brisbane" = "Site 1", "Sydney" = "Site 2" ,"Melbourne" = "Site 3")


dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
outdir <- paste0(dir,"/output/QC_metric_figs/latent/")
dir.create(outdir)


##### Read in Data #####
### seurat object ###
seurat <- readRDS(paste0(dir,"output/Distribution_tests/seurat_integrated_all_times_clustered_1pct_expressing.rds"))

### velocity metadata ###
velo_meta <- fread("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/scVelo/velocyto/scvelo_combat_corrected_all2/metadata.csv", sep = ",")
rownames(velo_meta) <- velo_meta$V1
velo_meta$V1 <- NULL

velo_meta_sub <- velo_meta[,c("Location", "n_unspliced_counts", "latent_time")]
rownames(velo_meta_sub) <- rownames(velo_meta)

seurat <- AddMetaData(seurat, velo_meta_sub)

seurat_noNA <- subset(seurat, subset = latent_time >=0)


df <- data.table(seurat_noNA@meta.data)

for (location in names(site_updates)){
	df$Location <- gsub(location, site_updates[location], df$Location)
}

df$Time <- ifelse(grepl("Baseline", df$Time), "Baseline", 
			ifelse(grepl("Thawed Village Day 0", df$Time), "Baseline", "Village"))

df$Location <- gsub("Site 2cryopreserved", "Site 2 Cryopreserved", df$Location)
df$Location <- factor(df$Location, levels = c("Site 1", "Site 2", "Site 3", "Site 2 Cryopreserved"))


plot <- ggplot(df, aes(latent_time, fill = Final_Assignment)) +
			geom_density(alpha = 0.6, size=1, linetype = "blank") +
			theme_classic() +
			facet_grid(Time ~ Location) +
			scale_fill_manual(values = cell_line_colors) +
			xlab("Pseudotime") +
			ylab("Density")
ggsave(plot, filename = paste0(outdir, "latent_area_plot.png"), width = 9, height = 2)
ggsave(plot, filename = paste0(outdir, "latent_area_plot.pdf"), width = 9, height = 2)



##### Make table of latent time averages for Joseph #####
all_data_4_joseph <- df[,c("Location","Final_Assignment", "Site_rep", "Time", "latent_time")]
table <- all_data_4_joseph %>%
     group_by(Location, Final_Assignment, Site_rep, Time) %>%
     summarise(mean = mean(latent_time), sd = sd(latent_time), N = n(), se = sd(latent_time)/sqrt(n()))

fwrite(all_data_4_joseph, paste0(outdir,"full_dataset_latent_time.csv"), sep = ",")
fwrite(table, paste0(outdir,"summary_table_latent_time.csv"), sep = ",")


##### ANOVA Testing#####
all_data_4_joseph <- fread(paste0(outdir,"full_dataset_latent_time.csv"), sep = ",")
table <- fread(paste0(outdir,"summary_table_latent_time.csv"), sep = ",")

aov1 <- aov(latent_time ~ Time + Final_Assignment, data = all_data_4_joseph[Location == "Site 1"])
aov2 <- aov(latent_time ~ Time + Final_Assignment, data = all_data_4_joseph[Location == "Site 2"])
aov3 <- aov(latent_time ~ Time + Final_Assignment, data = all_data_4_joseph[Location == "Site 3"])
aov2_cryo <- aov(latent_time ~ Time + Final_Assignment, data = all_data_4_joseph[Location == "Site 2 Cryopreserved"])
summary(aov1)
summary(aov2)
summary(aov3)
summary(aov2_cryo)


single_aov1_base <- aov(latent_time ~ Final_Assignment + Site_rep, data = all_data_4_joseph[Location == "Site 1" & Time == "Baseline"])
single_aov1_vill <- aov(latent_time ~ Final_Assignment + Site_rep, data = all_data_4_joseph[Location == "Site 1" & Time == "Village"])
single_aov2_base <- aov(latent_time ~ Final_Assignment + Site_rep, data = all_data_4_joseph[Location == "Site 2" & Time == "Baseline"])
single_aov2_vill <- aov(latent_time ~ Final_Assignment + Site_rep, data = all_data_4_joseph[Location == "Site 2" & Time == "Village"])
single_aov3_base <- aov(latent_time ~ Final_Assignment + Site_rep, data = all_data_4_joseph[Location == "Site 3" & Time == "Baseline"])
single_aov3_vill <- aov(latent_time ~ Final_Assignment + Site_rep, data = all_data_4_joseph[Location == "Site 3" & Time == "Village"])
single_aov2_cryo_base <- aov(latent_time ~ Final_Assignment + Site_rep, data = all_data_4_joseph[Location == "Site 2 Cryopreserved" & Time == "Baseline"])
single_aov2_cryo_vill <- aov(latent_time ~ Final_Assignment + Site_rep, data = all_data_4_joseph[Location == "Site 2 Cryopreserved" & Time == "Village"])
summary(single_aov1_base)
summary(single_aov1_vill)
summary(single_aov2_base)
summary(single_aov2_vill)
summary(single_aov3_base)
summary(single_aov3_vill)
summary(single_aov2_cryo_base)
summary(single_aov2_cryo_vill)



### With average ###
single_aov1_base_mean <- aov(mean ~ Final_Assignment + Site_rep, data = table[Location == "Site 1" & Time == "Baseline"])
single_aov1_vill_mean <- aov(mean ~ Final_Assignment + Site_rep, data = table[Location == "Site 1" & Time == "Village"])
single_aov2_base_mean <- aov(mean ~ Final_Assignment + Site_rep, data = table[Location == "Site 2" & Time == "Baseline"])
single_aov2_vill_mean <- aov(mean ~ Final_Assignment + Site_rep, data = table[Location == "Site 2" & Time == "Village"])
single_aov3_base_mean <- aov(mean ~ Final_Assignment + Site_rep, data = table[Location == "Site 3" & Time == "Baseline"])
single_aov3_vill_mean <- aov(mean ~ Final_Assignment + Site_rep, data = table[Location == "Site 3" & Time == "Village"])
single_aov2_cryo_base_mean <- aov(mean ~ Final_Assignment + Site_rep, data = table[Location == "Site 2 Cryopreserved" & Time == "Baseline"])
single_aov2_cryo_vill_mean <- aov(mean ~ Final_Assignment + Site_rep, data = table[Location == "Site 2 Cryopreserved" & Time == "Village"])
summary(single_aov1_base_mean)
summary(single_aov1_vill_mean)
summary(single_aov2_base_mean)
summary(single_aov2_vill_mean)
summary(single_aov3_base_mean)
summary(single_aov3_vill_mean)
summary(single_aov2_cryo_base_mean)
summary(single_aov2_cryo_vill_mean)

