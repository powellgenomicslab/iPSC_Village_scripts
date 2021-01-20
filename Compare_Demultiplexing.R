library("tidyr")
library("tidyverse")
library("dplyr")
library("ggplot2")
library("jcolors")

dir <- "/Volumes/ScratchRoot/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/"
outdir <- "/Volumes/ScratchRoot/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/CompareDemultiplexing/"

pools <- dir(paste0(dir,"popscle/"), pattern ="DRENEA")



mad_function <- function(df, column){
  df <- as.data.frame(df)
  print(head(df))
  mad <- mad(df[,column])
  median <- median(df[,column])
  print(paste("The mad of ", column," is:"))
  print(mad)
  print(paste("The median of ", column," is:"))
  print(median)

  df[,paste0(column,"_mad")] <- ifelse((df[,column] > (median - mad) & df[,column] < (median + mad)), "1_mad",
                                                 ifelse((df[,column] > (median - 2*mad) & df[,column] < (median + 2*mad)), "2_mad",
                                                        ifelse((df[,column] > (median - 3*mad) & df[,column] < (median + 3*mad)), "3_mad",
                                                               ifelse((df[,column] > (median - 4*mad) & df[,column] < (median + 4*mad)), "4_mad",
                                                                      ifelse((df[,column] > (median - 5*mad) & df[,column] < (median + 5*mad)), "5_mad", ">5_mad")))))
  return(df)
}


old_name <- c("22_FSA", "29_MBE", "36_TOB00421_i_E8")
new_name <- c("FSA0006","MBE1006","TOB0421")
cell_key <- data.frame(old_name, new_name)

##### Read in demultiplexing results from each of the softwares #####
#### scSplit ####
scSplit <- lapply(pools, FUN = function(x){
  read_delim(paste0(dir,"scSplit/",x,"/Count20/scSplit_result.csv"), delim = "\t")
})
names(scSplit) <- pools

scSplit <- lapply(scSplit, function(x){
  colnames(x) <- paste0("scSplit_",colnames(x))
  return(x)
})

#### souporcell ####
souporcell <- lapply(pools, FUN = function(x){
  read_delim(paste0(dir,"souporcell/",x,"/clusters.tsv"), delim = "\t")
})
names(souporcell) <- pools

souporcell <- lapply(souporcell, function(x){
  colnames(x) <- paste0("souporcell_",colnames(x))
  return(x)
})

#### vireo ####
vireo <- lapply(pools, FUN = function(x){
  read_delim(paste0(dir,"vireo/vireoGenotyped/MAF0.1/",x,"/donor_ids.tsv"), delim = "\t")
})
names(vireo) <- pools

vireo <- lapply(vireo, function(x){
  colnames(x) <- paste0("vireo_",colnames(x))
  return(x)
})

#### vireo no genotype ####
vireo_un <- lapply(pools, function(x){
  read_delim(paste0(dir,"vireo/vireoUngenotyped/MAF0.1/",x,"/donor_ids.tsv"), delim = "\t")
})
names(vireo_un) <- pools

vireo_un <- lapply(vireo_un, function(x){
  colnames(x) <- paste0("vireo_un_",colnames(x))
  return(x)
})

#### demuxlet ####
demuxlet <- lapply(pools, FUN = function(x){
  read_delim(paste0(dir,"popscle/demuxlet/",x,"/Imputed_GPdemuxletOUT.best"), delim = "\t")
})
names(demuxlet) <- pools

demuxlet <- lapply(demuxlet, function(x){
  colnames(x) <- paste0("demuxlet_",colnames(x))
  return(x)
})

#### freemuxlet ####
freemuxlet <- lapply(pools, FUN = function(x){
  read_delim(paste0(dir,"popscle/freemuxlet/",x,"/MAF0.1freemuxletOUT.clust1.samples.gz"), delim = "\t")
})
names(freemuxlet) <- pools

freemuxlet <- lapply(freemuxlet, function(x){
  colnames(x) <- paste0("freemuxlet_",colnames(x))
  return(x)
})

#### Metadata ####
meta <- read_delim(paste0(dir,"Seurat/QC/AllCells_metadata.txt"), delim = "\t")
meta <- meta[,c(1:4,12:14)]

meta <- mad_function(meta, "nCount_RNA")
meta <- mad_function(meta, "nFeature_RNA")
meta <- mad_function(meta, "percent.mt")
meta <- mad_function(meta, "percent.rb")

meta$nCount_RNA_mad <- factor(meta$nCount_RNA_mad, levels = c("1_mad","2_mad","3_mad","4_mad","5_mad",">5_mad"))
meta$nFeature_RNA_mad <- factor(meta$nFeature_RNA_mad, levels = c("1_mad","2_mad","3_mad","4_mad","5_mad",">5_mad"))
meta$percent.mt_mad <- factor(meta$percent.mt_mad, levels = c("1_mad","2_mad","3_mad","4_mad","5_mad",">5_mad"))
meta$percent.rb_mad <- factor(meta$percent.rb_mad, levels = c("1_mad","2_mad","3_mad","4_mad","5_mad",">5_mad"))

##### Join together the dataframes #####
joined_df <- lapply(names(scSplit), FUN = function(x){
  full_join(scSplit[[x]], souporcell[[x]], by = c("scSplit_Barcode" = "souporcell_barcode"))
})
names(joined_df) <- pools

joined_df <- lapply(names(joined_df), function(x){
  full_join(joined_df[[x]], freemuxlet[[x]], by = c("scSplit_Barcode" = "freemuxlet_BARCODE"))
})
names(joined_df) <- pools

joined_df <- lapply(names(joined_df), function(x){
  full_join(joined_df[[x]], vireo[[x]], by = c("scSplit_Barcode" = "vireo_cell"))
})
names(joined_df) <- pools

joined_df <- lapply(names(joined_df), function(x){
  full_join(joined_df[[x]], demuxlet[[x]], by = c("scSplit_Barcode" = "demuxlet_BARCODE"))
})
names(joined_df) <- pools

joined_df <- lapply(names(joined_df), function(x){
  full_join(joined_df[[x]], vireo_un[[x]], c("scSplit_Barcode" = "vireo_un_cell"))
})
names(joined_df) <- pools

##### Change names of the cell lines in the dataframes #####
joined_df <- lapply(joined_df, function(x){
  for (row in 1:nrow(cell_key)){
    print(cell_key[row,"old_name"])
    print(cell_key[row,"new_name"])
    x[] <- lapply(x, function(y){
      gsub(cell_key[row,"old_name"], cell_key[row,"new_name"], y)
    })
  }
  return(x)
  })


##### Rename each of the clusters and then combine into one column for identifying common assignment #####
joined_df <- lapply(joined_df, function(x){
  colnames(x) <- c("Barcode",colnames(x)[2:ncol(x)])
  x$scSplit_Cluster <- gsub("SNG-","scSplit_", x$scSplit_Cluster)
  x$souporcell_assignment <- paste0("souporcell_",x$souporcell_assignment)
  x$freemuxlet_BEST.GUESS <- paste0("freemuxlet_", x$freemuxlet_BEST.GUESS)
  x$vireo_donor_id <- paste0("vireo_", x$vireo_donor_id)
  x$vireo_un_donor_id <- paste0("vireo_un_", x$vireo_un_donor_id)
  x$demuxlet_BEST <- gsub("SNG-","demuxlet_",x$demuxlet_BEST)
  x$combined_assignments <- paste(x$scSplit_Cluster,x$souporcell_assignment,x$freemuxlet_BEST.GUESS,x$vireo_donor_id, x$demuxlet_BEST, x$vireo_un_donor_id, sep = "-")
  return(x)
})

joined_assignment_counts <- lapply(joined_df, function(x){
  df <- as.data.frame(t(table(x$combined_assignments)))
  df <- df[order(df$Freq, decreasing = TRUE),]
  return(df)
})

### Take just the top three to assign cells to those that were shared as singlets across all to reassign to a common assignment ###
joined_assignment_counts_top <- lapply(joined_assignment_counts, function(x){
  x$Var1 <- NULL
  x <- x[1:3,]
  return(x)
})

joined_assignment_counts_top <- lapply(joined_assignment_counts_top, function(x){
  df <- separate(x,col = Var2, into = c("scSplit","souporcell","freemuxlet","vireo","demuxlet","vireo_un"), sep = "-")
  df$common_assignment <- gsub("demuxlet_","",df$demuxlet)
  return(df)
})

joined_assignment_key <- lapply(joined_assignment_counts_top, function(x){
  df <- pivot_longer(x,cols = c("scSplit","souporcell","freemuxlet","vireo","demuxlet","vireo_un"), names_to = "software")
  df$Freq <- NULL
  return(df)
})
names(joined_assignment_key) <- pools


##### Make a dataframe for each pool that has barcode, assignment and software info rbind with other softwares
joined_df_min <- lapply(joined_df, function(x){
  df_scSplit <- x[,c("Barcode","scSplit_Cluster")]
  colnames(df_scSplit) <- c("Barcode","Assignment")
  df_scSplit$Software <- "scSplit"
  df_souporcell <- x[,c("Barcode","souporcell_assignment")]
  colnames(df_souporcell) <- c("Barcode","Assignment")
  df_souporcell$Software <- "souporcell"
  df_freemuxlet <- x[,c("Barcode","freemuxlet_BEST.GUESS")]
  colnames(df_freemuxlet) <- c("Barcode","Assignment")
  df_freemuxlet$Software <- "freemuxlet"
  df_vireo <- x[,c("Barcode", "vireo_donor_id")]
  colnames(df_vireo) <- c("Barcode","Assignment")
  df_vireo$Software <- "vireo"
  df_demuxlet <- x[,c("Barcode", "demuxlet_BEST")]
  colnames(df_demuxlet) <- c("Barcode", "Assignment")
  df_demuxlet$Software <- "demuxlet"
  df_vireo_un <- x[,c("Barcode","vireo_un_donor_id")]
  colnames(df_vireo_un) <- c("Barcode","Assignment")
  df_vireo_un$Software <- "vireo_un"
  df <- rbind(df_scSplit, df_souporcell, df_freemuxlet, df_vireo, df_demuxlet, df_vireo_un)
  return(df)
})
names(joined_df_min) <- pools


##### Add the combined donor information by left_joining
joined_df_min <- lapply(names(joined_df_min), function(x){
  df <- left_join(joined_df_min[[x]], joined_assignment_key[[x]], by = c("Assignment" = "value", "Software" = "software"))
  df$Pool <- x
  return(df)
})


joined_df_min <- lapply(joined_df_min, function(x){
  x$common_assignment[is.na(x$common_assignment)] <- "doublet/unassigned"
  return(x)
})
names(joined_df_min) <- pools

##### Combine the pool dataframes together to do a facet_grid plot
joined_df4facet <- do.call(rbind,joined_df_min)
joined_df4facet

plot <- ggplot(joined_df4facet, aes(common_assignment, fill = Software)) +
  geom_bar(position="dodge") +
  facet_wrap(~Pool, nrow = 2) +
  theme_classic() +
  scale_fill_jcolors(palette = "pal4") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot, filename = paste0(outdir,"6Pools_software_bar_plot.png"))
         

##### Add the assignments to the joined_df #####
joined_df_min_wide <- lapply(joined_df_min, function(x){
  pivot_wider(x, names_from = Software, values_from = c(Assignment,common_assignment))
})

joined_df <- lapply(names(joined_df), function(x){
  left_join(joined_df[[x]], joined_df_min_wide[[x]], by = c("Barcode"))
})
names(joined_df) <- pools

### Add pool information to Barcodes ###
joined_df <- lapply(joined_df, function(x){
  x$Barcode <- paste0(x$Pool, "_", x$Barcode)
  x$demuxlet_N.SNP[is.na(x$demuxlet_N.SNP)] <- 0 ### Note, demuxlet will not return cells that have no SNPs
  return(x)
})

### Rbind the pool dataframes together

joined_df_combined <- do.call(rbind, joined_df)


##### Add in unassigned to the common assignments #####
joined_df_combined$common_assignment_scSplit <- ifelse(is.na(joined_df_combined$scSplit_Cluster), "quality_removed", joined_df_combined$common_assignment_scSplit)
joined_df_combined$common_assignment_souporcell <- ifelse(joined_df_combined$souporcell_status == "unassigned", "unassigned", joined_df_combined$common_assignment_souporcell)
joined_df_combined$common_assignment_freemuxlet <- ifelse(joined_df_combined$freemuxlet_DROPLET.TYPE == "AMB", "unassigned", joined_df_combined$common_assignment_freemuxlet)
joined_df_combined$common_assignment_vireo <- ifelse(joined_df_combined$vireo_donor_id == "vireo_unassigned", "unassigned", joined_df_combined$common_assignment_vireo)
joined_df_combined$common_assignment_demuxlet <- ifelse(grepl("AMB",joined_df_combined$Assignment_demuxlet), "unassigned", 
                                                          ifelse(is.na(joined_df_combined$Assignment_demuxlet),"quality_removed",joined_df_combined$common_assignment_demuxlet))
joined_df_combined$common_assignment_vireo_un <- ifelse(joined_df_combined$vireo_un_donor_id == "vireo_un_unassigned", "unassigned", joined_df_combined$common_assignment_vireo_un)


joined_df_combined$common_assignment_scSplit <- ifelse(joined_df_combined$common_assignment_scSplit == "doublet/unassigned", "doublet", joined_df_combined$common_assignment_scSplit)
joined_df_combined$common_assignment_souporcell <- ifelse(joined_df_combined$common_assignment_souporcell == "doublet/unassigned", "doublet", joined_df_combined$common_assignment_souporcell)
joined_df_combined$common_assignment_freemuxlet <- ifelse(joined_df_combined$common_assignment_freemuxlet == "doublet/unassigned", "doublet", joined_df_combined$common_assignment_freemuxlet)
joined_df_combined$common_assignment_vireo <- ifelse(joined_df_combined$common_assignment_vireo == "doublet/unassigned", "doublet", joined_df_combined$common_assignment_vireo)
joined_df_combined$common_assignment_demuxlet <- ifelse(joined_df_combined$common_assignment_demuxlet == "doublet/unassigned", "doublet", joined_df_combined$common_assignment_demuxlet)
joined_df_combined$common_assignment_vireo_un <- ifelse(joined_df_combined$common_assignment_vireo_un == "doublet/unassigned", "doublet", joined_df_combined$common_assignment_vireo_un)

### Combine dataframe with metadata ###
joined_df_combined$Barcode <- gsub("-1","",joined_df_combined$Barcode)
joined_df_combined_meta <- left_join(joined_df_combined, meta, by = "Barcode")

joined_df_combined_meta4facet <- pivot_longer(joined_df_combined_meta,cols = c(vireo_n_vars, vireo_un_n_vars, demuxlet_N.SNP, freemuxlet_NUM.SNPS), names_to = "software", values_to = "n_vars")
joined_df_combined_meta4facet <- full_join(joined_df_combined_meta4facet,joined_df_combined_meta4facet, by = colnames(joined_df_combined_meta4facet)[1:81])
joined_df_combined_meta4facet$n_vars.x <- as.numeric(as.character(joined_df_combined_meta4facet$n_vars.x))
joined_df_combined_meta4facet$n_vars.y <- as.numeric(as.character(joined_df_combined_meta4facet$n_vars.y))

##### Plots #####
color_vars_discrete <- c("ScrubletDoublet")
color_vars_common_assign <- c("common_assignment_vireo_un","common_assignment_vireo","common_assignment_demuxlet", "common_assignment_freemuxlet","common_assignment_souporcell","common_assignment_scSplit")
color_vars_discrete_seq <- c("nCount_RNA_mad", "nFeature_RNA_mad", "percent.mt_mad","percent.rb_mad")
color_vars_continuous <- c("percent.mt", "percent.rb","nFeature_RNA", "nCount_RNA")

common_assignment_colors <- c(jcolors("default"),"sienna")
names(common_assignment_colors) <- c("FSA006","MBE1006","doublet","TOB421","unassigned","quality_removed")

for (var in color_vars_discrete){
  plot <- ggplot(joined_df_combined_meta4facet, aes(log(n_vars.x + 1), log(n_vars.y +1))) +
    geom_point(alpha = 0.5, aes_string(color = var), size = 0.7) +
    theme_classic() +
    facet_grid(cols = vars(software.x), rows = vars(software.y), scales = "free") +
    scale_color_jcolors("pal9") +
    theme(text = element_text(size=25))
  ggsave(plot, filename = paste0(outdir,"compare_N_vars_",var,".png"),width = 12, height = 9)
}

for (var in color_vars_common_assign){
  plot <- ggplot(joined_df_combined_meta4facet, aes(log(n_vars.x + 1), log(n_vars.y +1))) +
    geom_point(alpha = 0.5, aes_string(color = var), size = 0.7) +
    theme_classic() +
    facet_grid(cols = vars(software.x), rows = vars(software.y), scales = "free") +
    scale_color_manual(values = common_assignment_colors) +
    theme(text = element_text(size=25))
  ggsave(plot, filename = paste0(outdir,"compare_N_vars_",var,".png"),width = 12, height = 9)
}

for (var in color_vars_continuous){
  plot <- ggplot(joined_df_combined_meta4facet, aes(log(n_vars.x+1), log(n_vars.y+1))) +
    geom_point(alpha = 0.5, aes_string(color = var), size = 0.7) +
    theme_classic() +
    facet_grid(cols = vars(software.x), rows = vars(software.y), scales = "free") +
    scale_color_jcolors_contin(palette = "pal12") +
    theme(text = element_text(size=25))
  ggsave(plot, filename = paste0(outdir,"compare_N_vars_",var,".png"),width = 12, height = 9)
}

colors <- jcolors("pal12")[c(1,4,6,8,11,13)]
for (var in color_vars_discrete_seq){
  plot <- ggplot(joined_df_combined_meta4facet, aes(log(n_vars.x+1), log(n_vars.y+1))) +
    geom_point(alpha = 0.5, aes_string(color = var), size = 0.7) +
    theme_classic() +
    facet_grid(cols = vars(software.x), rows = vars(software.y), scales = "free") +
    scale_color_manual(values = colors) +
    theme(text = element_text(size=25))
  ggsave(plot, filename = paste0(outdir,"compare_N_vars_",var,".png"),width = 12, height = 9)
}



### Plot probability of doublets vs singlets ###
plot <- ggplot(joined_df_combined_meta4facet, aes(souporcell_log_prob_singleton, souporcell_log_prob_doublet)) +
  geom_point() +
  theme_classic() +
  facet_grid(cols = vars(software.x), rows = vars(software.y), scales = "free") +
  scale_color_manual(values = colors) +
  theme(text = element_text(size=25))
ggsave(plot, filename = paste0(outdir,"compare_N_probs_souporcell_probabilities.png"),width = 12, height = 9)



### Plot number of doublets vs QC ###
joined_df_combined_meta$FSA006 <- NA
joined_df_combined_meta$MBE1006 <- NA
joined_df_combined_meta$TOB421 <- NA
joined_df_combined_meta$doublet <- NA
joined_df_combined_meta$unassigned <- NA
joined_df_combined_meta$quality_removed <- NA

for (row in 1:nrow(joined_df_combined_meta)){
  joined_df_combined_meta$FSA006[row] <- length(which(joined_df_combined_meta[row,c(grep("common_assignment",colnames(joined_df_combined_meta)))] == "FSA006"))
  joined_df_combined_meta$MBE1006[row] <- length(which(joined_df_combined_meta[row,c(grep("common_assignment",colnames(joined_df_combined_meta)))] == "MBE1006"))
  joined_df_combined_meta$TOB421[row] <- length(which(joined_df_combined_meta[row,c(grep("common_assignment",colnames(joined_df_combined_meta)))] == "TOB421"))
  joined_df_combined_meta$doublet[row] <- length(which(joined_df_combined_meta[row,c(grep("common_assignment",colnames(joined_df_combined_meta)))] == "doublet"))
  joined_df_combined_meta$unassigned[row] <- length(which(joined_df_combined_meta[row,c(grep("common_assignment",colnames(joined_df_combined_meta)))] == "unassigned"))
  joined_df_combined_meta$quality_removed[row] <- length(which(joined_df_combined_meta[row,c(grep("common_assignment",colnames(joined_df_combined_meta)))] == "quality_removed"))
}


joined_df_combined_meta_long <- pivot_longer(joined_df_combined_meta, cols = c(FSA006, MBE1006, TOB421, doublet, unassigned, quality_removed), names_to = "cell_assignment", values_to = "n_softwares")
joined_df_combined_meta_long$cell_assignment <- factor(joined_df_combined_meta_long$cell_assignment, levels = c("FSA006", "MBE1006","TOB421","doublet","unassigned","quality_removed"))



### Bar plots for each QC metric by number of souftware identified groups ###
for (var in color_vars_discrete){
  plot <- ggplot(joined_df_combined_meta_long, aes(n_softwares)) +
    geom_bar(aes_string(fill = var),position = "fill") +
    theme_classic() +
    scale_fill_jcolors(palette = "pal9") +
    theme(text = element_text(size=25)) +
    facet_wrap(facets =vars(cell_assignment))
  ggsave(plot, filename = paste0(outdir,"compare_N_vars_facet_",var,".png"),width = 12, height = 9)
}






### Bar plots for each QC metric by number of souftware identified doublets ###
for (var in color_vars_common_assign){
  plot <- ggplot(joined_df_combined_meta_long, aes(n_softwares)) +
    geom_bar(aes_string(fill = var),position = "fill") +
    theme_classic() +
    scale_fill_manual(values = common_assignment_colors) +
    theme(text = element_text(size=25)) +
    facet_wrap(facets =vars(cell_assignment))
  ggsave(plot, filename = paste0(outdir,"compare_N_vars_facet_",var,".png"),width = 12, height = 9)
}


colors <- jcolors("pal12")[c(1,4,6,8,11,13)]
for (var in color_vars_discrete_seq){
  plot <- ggplot(joined_df_combined_meta_long, aes(n_softwares)) +
    geom_bar(aes_string(fill = var),position = "fill") +
    theme_classic() +
    xlim(-1, 7) +
    scale_fill_manual(values = colors) +
    theme(text = element_text(size=25)) +
    facet_wrap(facets =vars(cell_assignment))
  ggsave(plot, filename = paste0(outdir,"compare_N_vars_facet_",var,".png"),width = 12, height = 9)
}



##### Make table to use for metadata #####
joined_df_combined_meta$FinalAssignment <- ifelse(joined_df_combined_meta$FSA006 >= 4, "FSA006", 
                                                  ifelse(joined_df_combined_meta$MBE1006 >= 4, "MBE1006", 
                                                         ifelse(joined_df_combined_meta$TOB421 >= 4, "TOB421", 
                                                                ifelse(joined_df_combined_meta$doublet >= 4, "doublet", 
                                                                       ifelse(joined_df_combined_meta$unassigned >= 4, "unassigned",
                                                                              ifelse(joined_df_combined_meta$quality_removed >= 4, "quality_removed","combination_fail_remove"))))))
write_delim(joined_df_combined_meta, paste0(outdir,"doublet_metadata.txt"), delim = "\t")



