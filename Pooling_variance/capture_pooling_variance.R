##### Author: Drew Neavin
##### Date: 28 November, 2022
##### Goal: Use the fibroblast data to estimate the pooling variance compared to distribution of cells observed in experiment
##### Fibroblasts 
##### OneK1K
##### This data


library(data.table)
library(tidyverse)
library(ggplot2)
library(Seurat)
library(SeuratDisk)



##### Set up directories #####
datadir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Pooling_variance/"
indir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/CombinedResults/"
prop_dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/preQC/"
amd_dir <- "/directflow/SCCGGroupShare/projects/annsen/analysis/AMD_scRNA/demultiplexing/"
poag_dir <- "/directflow/SCCGGroupShare/projects/annsen/analysis/POAG_scRNA/demuxlet/"
jerber_dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/pool_sd_publicly_available_datasets/Jerber_DN/"
outdir <- paste0(datadir,"figures/")
card_dir <- "/directflow/SCCGGroupShare/projects/nonfar/analysis/cardiac_multiome_directflow/demux_obj/"
dir.create(outdir)



##### Get a list of the village pools #####

onek_pools <- dir(datadir, pattern = "OneK1K_scRNA_Sample")
fibro_pools <- dir(datadir, pattern = "scFibroblast_EQTL_Sample")
village_pools <- c(dir(datadir, pattern = "DRENEA"), dir(datadir, pattern = "Village"))
amd_pools <- dir(amd_dir, pattern = "RP_pool")
poag_pools <- list.files(poag_dir, pattern = "_Demuxlet_ExonOnly.best") %>% gsub("_Demuxlet_ExonOnly.best", "", .)
card_pools <- list.files(card_dir, pattern = "Pool_Day") %>% grep("DemuxALL", ., value = TRUE)



##### Read in results #####
onek1k_results_list <- lapply(onek_pools, function(x){
    fread(paste0(datadir, x, "/combined_assignments_w_combined_assignments.tsv"), sep = "\t")
})
names(onek1k_results_list) <- onek_pools

fibro_results_list <- lapply(fibro_pools, function(x){
    fread(paste0(datadir, x, "/combined_assignments_w_combined_assignments.tsv"), sep = "\t")
})
names(fibro_results_list) <- fibro_pools

### Because of hashing, need to get from seurat objects inseatd
village_seurat <- readRDS("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Seurat_w_Hash_freeze_thaw/Demultiplex_Doublet/Seurat_noDoublets_norm_meta.rds")

amd_results_list <- lapply(amd_pools, function(x){
    fread(paste0(amd_dir, x, "/FinalAssignments.tsv"), sep = "\t")
})
names(amd_results_list) <- amd_pools

poag_results_list <- lapply(poag_pools, function(x){
    fread(paste0(poag_dir, x, "_Demuxlet_ExonOnly.best"), sep = "\t")
})
names(poag_results_list) <- poag_pools


card_results_list <- lapply(card_pools, function(x){
    print(x)
    tmp <- readRDS(paste0(card_dir,x))
    dt <- data.table(tmp@meta.data)
    dt$Pool <- gsub("_DemuxALL.rds", "",x)
    return(dt)
})
names(card_results_list) <- gsub("_DemuxALL.rds", "", card_pools)


cuomo_ipsc <- fread("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/pool_sd_publicly_available_datasets/Cuomo_iPSC/cell_metadata_cols.tsv", sep = "\t")
diff_props <- fread(paste0(prop_dir, "cardiac_diff_prop_lines.tsv"), sep = "\t")
multi_passage_props <- fread(paste0(prop_dir, "multi_passage_prop_lines.tsv"), sep = "\t")


jerber <- fread(paste0(jerber_dir,"cell_numbers.tsv"), sep = "\t")




##### Get Numbers #####
onek1k_results_combined_assignments_list <- lapply(onek1k_results_list, function(x){
    x <- x[,c("MajoritySinglet_DropletType", "MajoritySinglet_Individual_Assignment")]
    colnames(x) <- c("DropletType", "Assignment")
    return(x)
})

fibro_results_combined_assignments_list <- lapply(fibro_results_list, function(x){
    x <- x[,c("MajoritySinglet_DropletType", "MajoritySinglet_Individual_Assignment")]
    colnames(x) <- c("DropletType", "Assignment")
    return(x)
})


amd_results_combined_assignments_list <- lapply(amd_results_list, function(x){
    x <- x[,c("DropletType", "Assignment")]
    return(x)
})


poag_results_combined_assignments_list <- lapply(poag_results_list, function(x){
    data.table(DropletType = ifelse(grepl("DBL", x$BEST), "doublet", ifelse(grepl("SNG", x$BEST), "singlet", "unassigned")), Assignment = ifelse(grepl("DBL", x$BEST), "doublet", ifelse(grepl("SNG", x$BEST), gsub("SNG-", "", x$BEST) %>% gsub("_.+", "", .), "unassigned")))
})


card_results_combined_assignments_list <- lapply(card_results_list, function(x){
    x[,c("Pool", "Assignment")]
})




##### Get the proportions #####
onek1k_table_list <- lapply(names(onek1k_results_combined_assignments_list), function(x){
    tmp <- data.table(data.frame(prop.table(table(onek1k_results_combined_assignments_list[[x]][!(Assignment %in% c("doublet", "unassigned"))]$Assignment))))
    colnames(tmp) <- c("Assignment", "N")
    tmp$Pool <- x
    tmp$Group <- "OneK1K\nPooled"
    return(tmp)
})


fibro_table_list <- lapply(names(fibro_results_combined_assignments_list), function(x){
    tmp <- data.table(data.frame(prop.table(table(fibro_results_combined_assignments_list[[x]][!(Assignment %in% c("doublet", "unassigned"))]$Assignment))))
    colnames(tmp) <- c("Assignment", "N")
    tmp$Pool <- x
    tmp$Group <- "Fibroblasts\nPooled"
    return(tmp)
})


card_table_list <- lapply(names(card_results_combined_assignments_list), function(x){
    tmp <- data.table(data.frame(prop.table(table(card_results_combined_assignments_list[[x]][!(Assignment %in% c("doublet", "unassigned"))]$Assignment))))
    colnames(tmp) <- c("Assignment", "N")
    tmp$Group <- "Cardiomyocyte\nDifferentiation\nPooled"
    tmp$Pool <- x
    return(tmp)
})


amd_table_list <- lapply(names(amd_results_combined_assignments_list), function(x){
    tmp <- data.table(data.frame(prop.table(table(amd_results_combined_assignments_list[[x]][!(Assignment %in% c("doublet", "unassigned"))]$Assignment))))
    colnames(tmp) <- c("Assignment", "N")
    tmp$Pool <- x
    tmp$Group <- "AMD\nPooled"
    return(tmp)
})


poag_table_list <- lapply(names(poag_results_combined_assignments_list), function(x){
    tmp <- data.table(data.frame(prop.table(table(poag_results_combined_assignments_list[[x]][!(Assignment %in% c("doublet", "unassigned"))]$Assignment))))
    colnames(tmp) <- c("Assignment", "N")
    tmp$Pool <- x
    tmp$Group <- "POAG\nPooled"
    return(tmp)
})



cuomo_table_list <- lapply(unique(paste0(cuomo_ipsc$day, cuomo_ipsc$experiment)), function(x){
    tmp <- data.table(data.frame(prop.table(table(cuomo_ipsc[paste0(day, experiment) == x]$donor))))
    colnames(tmp) <- c("Assignment", "N")
    tmp$Pool <- x
    tmp$Group <- "Defendo\nVillage"
    return(tmp)
})


village_table_list <- lapply(unique(paste0(village_seurat@meta.data$Site_rep, village_seurat@meta.data$Time)), function(x){
    dt <- data.table(village_seurat@meta.data)
    tmp <- data.table(data.frame(prop.table(table(dt[paste0(Site_rep, Time) == x]$Final_Assignment))))
    colnames(tmp) <- c("Assignment", "N")
    tmp$Pool <- x
    tmp$Group <- "3-line Village"
    tmp <- tmp[!(Assignment %in% c("doublet", "unassigned"))]
    return(tmp)
})


jerber_table_list <- lapply(unique(paste0(jerber$pool_id, jerber$time_point)), function(x){
    tmp <- jerber[, sum(n_cells), by=.(pool_id,time_point, cell_line)][paste0(pool_id, time_point) == x][,c("cell_line", "V1")]
    tmp[ , prop := V1/sum(V1), ]
    tmp$V1 <- NULL
    colnames(tmp) <- c("Assignment", "N")
    tmp$Pool <- x
    tmp$Group <- "DN\nDifferentiation\nVillage"
    return(tmp)
})


##### Get standard deviation #####
onek1k_sd <- lapply(onek1k_table_list, function(x){
    data.table(Pool = unique(x$Pool), Group = unique(x$Group), sd = sd(x$N))
})


fibro_sd <- lapply(fibro_table_list, function(x){
    data.table(Pool = unique(x$Pool), Group = unique(x$Group), sd = sd(x$N))
})


village_sd <- lapply(village_table_list, function(x){
    data.table(Pool = unique(x$Pool), Group = unique(x$Group), sd = sd(x$N))
})


amd_sd <- lapply(amd_table_list, function(x){
    data.table(Pool = unique(x$Pool), Group = unique(x$Group), sd = sd(x$N))
})


poag_sd <- lapply(poag_table_list, function(x){
    data.table(Pool = unique(x$Pool), Group = unique(x$Group), sd = sd(x$N))
})


cuomo_sd <- lapply(cuomo_table_list, function(x){
    data.table(Pool = unique(x$Pool), Group = unique(x$Group), sd = sd(x$N))
})


jerber_sd <- lapply(jerber_table_list, function(x){
    data.table(Pool = unique(x$Pool), Group = unique(x$Group), sd = sd(x$N))
})


card_sd <- lapply(card_table_list, function(x){
    data.table(Pool = unique(x$Pool), Group = unique(x$Group), sd = sd(x$N))
})


diff_sd <- lapply(unique(diff_props$Day), function(x){
    data.table(Pool = paste0("Cardiac_diff_Day_",x), Group = "Cardiomyocyte\nDifferentiation\nVillage", sd = sd(diff_props[Day == x]$N))
})


multi_passage_sd <- lapply(unique(multi_passage_props$Day), function(x){
    data.table(Pool = paste0("Multi_passage_P",x), Group = "Multi-passage\nVillage", sd = sd(multi_passage_props[Day == x]$N))
})



##### Join datasets together
Proportions_dt <- do.call(rbind, c(onek1k_sd, fibro_sd, diff_sd, multi_passage_sd, amd_sd, poag_sd,cuomo_sd, village_sd, jerber_sd, card_sd))



##### Make a plot for median #####
sd_median <- Proportions_dt[, median(sd), by = .(Group)]


##### Plot distributions of each + this data #####
p_bar <- ggplot(Proportions_dt[Group != "3-line Village"], aes(sd, fill = Group)) +
    geom_histogram(position="identity", alpha = 0.5) +
    theme_classic() +
    facet_wrap(vars(factor(Group, levels = c("OneK1K\nPooled", "Fibroblasts\nPooled", "AMD\nPooled", "POAG\nPooled","Cardiomyocyte\nDifferentiation\nPooled", "Defendo\nVillage", "DN\nDifferentiation\nVillage", "Cardiomyocyte\nDifferentiation\nVillage", "Multi-passage\nVillage"))), ncol = 1, scales = "free_y", strip.position = "right") +
    xlab("Standard Deviation") +
    ylab("Count") +
    geom_vline(data = sd_median, aes(xintercept = V1, color = Group))

ggsave(p_bar, filename = paste0(outdir,"barplot_sd.png"), width = 4, height = 10)
ggsave(p_bar, filename = paste0(outdir,"barplot_se.png"), width = 4, height = 10)




##### Get standard deviation #####
onek1k_se <- lapply(onek1k_table_list, function(x){
    data.table(Pool = unique(x$Pool), Group = unique(x$Group), se = sd(x$N)/sqrt(length(x$N)))
})


fibro_se <- lapply(fibro_table_list, function(x){
    data.table(Pool = unique(x$Pool), Group = unique(x$Group), se = sd(x$N)/sqrt(length(x$N)))
})


village_se <- lapply(village_table_list, function(x){
    data.table(Pool = unique(x$Pool), Group = unique(x$Group), se = sd(x$N)/sqrt(length(x$N)))
})


amd_se <- lapply(amd_table_list, function(x){
    data.table(Pool = unique(x$Pool), Group = unique(x$Group), se = sd(x$N)/sqrt(length(x$N)))
})


poag_se <- lapply(poag_table_list, function(x){
    data.table(Pool = unique(x$Pool), Group = unique(x$Group), se = sd(x$N)/sqrt(length(x$N)))
})


cuomo_se <- lapply(cuomo_table_list, function(x){
    data.table(Pool = unique(x$Pool), Group = unique(x$Group), se = sd(x$N)/sqrt(length(x$N)))
})


jerber_se <- lapply(jerber_table_list, function(x){
    data.table(Pool = unique(x$Pool), Group = unique(x$Group), se = sd(x$N)/sqrt(length(x$N)))
})


card_se <- lapply(card_table_list, function(x){
    data.table(Pool = unique(x$Pool), Group = unique(x$Group), se = sd(x$N)/sqrt(length(x$N)))
})


diff_se <- lapply(unique(diff_props$Day), function(x){
    data.table(Pool = paste0("Cardiac_diff_Day_",x), Group = "Cardiomyocyte\nDifferentiation\nVillage", se = sd(diff_props[Day == x]$N)/sqrt(length(diff_props[Day == x]$N)))
})


multi_passage_se <- lapply(unique(multi_passage_props$Day), function(x){
    data.table(Pool = paste0("Multi_passage_P",x), Group = "Multi-passage\nVillage", se = sd(multi_passage_props[Day == x]$N)/sqrt(length(multi_passage_props[Day == x]$N)))
})



##### Join datasets together
Proportions_dt <- do.call(rbind, c(onek1k_se, fibro_se, diff_se, multi_passage_se, amd_se, poag_se,cuomo_se, village_se, jerber_se, card_se))



##### Make a plot for median #####
se_median <- Proportions_dt[, median(se), by = .(Group)]
se_mean <- Proportions_dt[, mean(se), by = .(Group)]


##### Plot distributions of each + this data #####
p_bar <- ggplot(Proportions_dt[Group != "3-line Village"], aes(se, fill = factor(Group, levels = c("OneK1K\nPooled", "Fibroblasts\nPooled", "AMD\nPooled", "POAG\nPooled","Cardiomyocyte\nDifferentiation\nPooled", "Defendo\nVillage", "DN\nDifferentiation\nVillage", "Cardiomyocyte\nDifferentiation\nVillage", "Multi-passage\nVillage")))) +
    geom_histogram(position="identity", alpha = 0.5) +
    theme_classic() +
    facet_wrap(vars(factor(Group, levels = c("OneK1K\nPooled", "Fibroblasts\nPooled", "AMD\nPooled", "POAG\nPooled","Cardiomyocyte\nDifferentiation\nPooled", "Defendo\nVillage", "DN\nDifferentiation\nVillage", "Cardiomyocyte\nDifferentiation\nVillage", "Multi-passage\nVillage"))), ncol = 1, scales = "free_y", strip.position = "right") +
    xlab("Standard Error") +
    ylab("Count") +
    geom_vline(data = se_mean[Group != "3-line Village"], aes(xintercept = V1, color = factor(Group, levels = c("OneK1K\nPooled", "Fibroblasts\nPooled", "AMD\nPooled", "POAG\nPooled","Cardiomyocyte\nDifferentiation\nPooled", "Defendo\nVillage", "DN\nDifferentiation\nVillage", "Cardiomyocyte\nDifferentiation\nVillage", "Multi-passage\nVillage"))), size = 1, linetype = "dashed") +
    theme(legend.position = "none")

ggsave(p_bar, filename = paste0(outdir,"barplot_se.png"), width = 3, height = 10)
ggsave(p_bar, filename = paste0(outdir,"barplot_se.pdf"), width = 3, height = 10)




fwrite(Proportions_dt[Group != "3-line Village"], paste0(outdir, "pooling_variances.tsv"), sep = "\t")



