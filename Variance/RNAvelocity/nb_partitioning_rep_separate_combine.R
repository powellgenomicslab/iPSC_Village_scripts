library(tidyverse)
library(data.table)
library(Seurat)



dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
datadir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/RNAvelocity/nb_partitioning_rep_separate/gene_separated/"
outdir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/RNAvelocity/nb_partitioning_rep_separate/combined/"
dir.create(outdir)




save_figs <- function(plot, basename, width = 17, height = 17, units = "cm"){
    ggsave(plot, filename = paste0(basename,".png"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".pdf"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".eps"), height = height, width = width, units = units)
}


##### Set up colors #####
variable_colors <- c(Village = "#A2B0D0", Replicate = "#64A66B", Line = "#68319B") 
line_colors <- c(FSA0006 = "#F79E29", MBE1006 = "#9B2C99", TOB0421 = "#35369C")

site_updates <- c("Brisbane" = "Site 1", "Sydney_cryopreserved" = "Site 2 Cryopreserved", "Sydney" = "Site 2" ,"Melbourne" = "Site 3")





if (!file.exists(paste0(outdir,"nb_variance_partitionin_df.tsv"))){
	##### Get list of files #####
	dir_list <- (datadir)
	dir_list <- dir_list[!(dir_list %in% "logs")]
	fits_files <- list.files(paste0(datadir, dir_list), pattern = "_variable_variance_partition_coefficients.rds")
	names(fits_files) <- gsub("_variable_variance_partition_coefficients.rds", "", fits_files)


	fits_files_list <- list()
	for (rep in dir_list){
		fits_files_list[[rep]] <- fits_files[which(grepl(rep, fits_files))]
	}


	##### Read in files #####
	fits_list <- lapply(names(fits_files_list), function(x){
		lapply(fits_files_list[[x]], function(y){
			readRDS(paste0(datadir,x,"/",y))
		})
	})
	names(fits_list) <- names(fits_files_list)
	

	fits <- lapply(names(fits_list), function(x){
		temp <- lapply(names(fits_list[[x]]), function(y){
			fits_list[[x]][[y]]$Replicate <- gsub("_ENSG.+","",y)
			return(fits_list[[x]][[y]])
		})
		binded <- do.call(rbind, temp)
		return(binded)
	})

	##### Combine fits Results #####
	fits_df <- do.call(rbind, fits)

	write_delim(fits_df, paste0(outdir,"nb_variance_partitionin_df.tsv"), delim = "\t")

} else {
	fits_df <- read_delim(paste0(outdir,"nb_variance_partitionin_df.tsv"), delim = "\t")
}

##### Remove faulty model genes and locations #####
### Read in faulty log files ###
faulty <- fread(paste0(datadir,"faulty_runs.txt"), header = FALSE)


### Read in the seurat data objects
seurat_list <- list()
for (replicate in unique(fits_df$Replicate)){
	seurat_list[[replicate]] <- readRDS(paste0("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/RNAvelocity/nb_partitioning_rep_separate/data/", replicate, "_seurat_1pct.rds"))
}

### Make dataframe of location x number for faulty models ###
faulty_df <- data.table(Replicate = gsub("_nb\\.o.+", "", faulty$V1), Number = gsub(".+nb\\.o\\d+\\.", "", faulty$V1))


## Add gene id to the faulty dataframe
faulty_df$ENSG <- NA
for (row in 1:nrow(faulty_df)){
	print(row)
	faulty_df$ENSG[row] <- rownames(seurat_list[[faulty_df$Replicate[row]]])[as.numeric(faulty_df$Number[row])] 
}


### Remove the faulty genes ###
fits_df <- setDT(fits_df)[!faulty_df, on = c("Replicate", "Gene" = "ENSG")]



##### Make long dataframe for plotting #####
fits_df_long <- pivot_longer(fits_df, cols = c("Line", "Village"), names_to = "Covariate", values_to = "Variance_Explained")
fits_df_long$Variance_Explained <- round(fits_df_long$Variance_Explained,6)
fits_df_long$Location <- gsub("\\d_Q\\d","",fits_df_long$Replicate)

fits_df_long <- left_join(fits_df_long, GeneConversion, by = c("Gene" = "ENSG_ID"))
fits_df_long$Location_Replicate_Quintile <- fits_df_long$Replicate



##### Update Site names #####
for (location in names(site_updates)){
	print(location)
	fits_df_long$Location <- gsub(location, site_updates[location], fits_df_long$Location)
	fits_df_long$Location_Replicate_Quintile <- gsub(location, paste0(site_updates[location], " Replicate "), fits_df_long$Location_Replicate_Quintile)
	fits_df_long$Replicate <- gsub(location, "Replicate ", fits_df_long$Replicate)
}

fits_df_long$Replicate <- gsub("_.+", "", fits_df_long$Replicate)
fits_df_long$Location_Replicate_Quintile <- gsub("_", " - ", fits_df_long$Location_Replicate_Quintile)
fits_df_long$Quintile <- gsub(".+- Q", "Q", fits_df_long$Location_Replicate_Quintile)



### Variance explained correlation between sites ###
fits_df_long$Replicate <- gsub(" - Q\\d")
fits_df_long <- data.table(fits_df_long)

fwrite(fits_df_long, sep = "\t", paste0(outdir, "var_explained_long.tsv"))



fits_df_long <- fread(sep = "\t", paste0(outdir, "var_explained_long.tsv"))

#### Check genes of interest ####
pluri_genes <- read_delim(paste0(dir,"data/pluripotency_genes.tsv"), delim = "\t", col_names = "Gene_ID")


##### Add gene IDs for easy identification downstream #####
GeneConversion1 <- read_delim("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/DRENEA_1/outs/filtered_feature_bc_matrix/features.tsv.gz", col_names = F, delim = "\t")
GeneConversion2 <- read_delim("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/Village_B_1_week/outs/filtered_feature_bc_matrix/features.tsv.gz", col_names = F, delim = "\t")

GeneConversion <- unique(rbind(GeneConversion1, GeneConversion2))
GeneConversion <- GeneConversion[!duplicated(GeneConversion$X1),]
GeneConversion$X3 <- NULL
colnames(GeneConversion) <- c("ENSG_ID", "Gene_ID")

pluri_genes <- left_join(pluri_genes, GeneConversion, by = "Gene_ID")


fits_df_long_pluri_genes <- left_join(pluri_genes, fits_df_long, by = c("ENSG_ID" = "Gene", "Gene_ID"))
fits_df_long_pluri_genes <- na.omit(fits_df_l ong_pluri_genes)

fits_df_long_pluri_genes_avg <- fits_df_long_pluri_genes %>% group_by(Gene_ID, ENSG_ID, Covariate, Location) %>% summarize(Variance_Explained = mean(Variance_Explained))

fits_df_long_pluri_genes_main_avg$Location <- factor(fits_df_long_pluri_genes_main_avg$Location, levels = c("Site 1", "Site 2", "Site 2 Cryopreserved", "Site 3"))

pPluri_Genes_Cont <- ggplot() +
						geom_bar(data = fits_df_long_pluri_genes_avg, aes(Location, Variance_Explained*100, fill = Covariate), position = "dodge", stat = "identity") +
						geom_point(data = fits_df_long_pluri_genes, aes(Location, Variance_Explained*100, fill = Covariate), size = 0.5, position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0)) +
						theme_classic() +
						facet_wrap(Gene_ID ~ Quintile) +
						# scale_color_manual(values = variable_colors) +
						scale_fill_manual(values = variable_colors) +
						theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1)) +
						ylab("Percent Gene Expression\nVariance Explained") +
						theme(axis.title.x=element_blank())
save_figs(pPluri_Genes_Cont, paste0(outdir, "Pluripotent_Gene_Variable_Contributions"), width = 50, height = 10)


##### Only use main markers for main figure #####
fits_df_long_pluri_genes_main <- fits_df_long_pluri_genes[which(fits_df_long_pluri_genes$Gene_ID %in% c("NANOG", "POU5F1", "SOX2", "MYC")),]

fits_df_long_pluri_genes_main_avg <- fits_df_long_pluri_genes_main %>% group_by(Gene_ID, ENSG_ID, Covariate, Location, Quintile) %>% summarize(Variance_Explained = mean(Variance_Explained))

fits_df_long_pluri_genes_main_avg$Location <- factor(fits_df_long_pluri_genes_main_avg$Location, levels = c("Site 1", "Site 2", "Site 2 Cryopreserved", "Site 3"))


pPluri_Genes_Cont <- ggplot() +
						geom_bar(data = fits_df_long_pluri_genes_main_avg, aes(Quintile, Variance_Explained*100, fill = Covariate), position = "dodge", stat = "identity", width  = 0.75) +
						geom_point(data = fits_df_long_pluri_genes_main, aes(Quintile, Variance_Explained*100, fill = Covariate), size = 0.5, position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0)) +
						theme_classic() +
						facet_grid(Gene_ID ~ Location, scales = "free") +
						# coord_flip() +
						scale_fill_manual(values = variable_colors) +
						ylab("Percent Gene Expression\nVariance Explained") +
						theme(axis.title.y=element_blank(),
						legend.position = "none") 
save_figs(pPluri_Genes_Cont, paste0(outdir, "Pluripotent_Gene_Variable_Contributions_main"), width = 20, height = 20)





### Main
fits_df_long_pluri_genes_main_avg_cryo <- fits_df_long_pluri_genes_main_avg[grepl("Sydney", fits_df_long_pluri_genes_main_avg$Location),]
fits_df_long_pluri_genes_main_avg_cryo$Location <- factor((gsub("Sydney_cryopreserved", "Cryopreserved", fits_df_long_pluri_genes_main_avg_cryo$Location) %>% gsub("Sydney", "Fresh",.)), levels = c("Cryopreserved", "Fresh"))

fits_df_long_pluri_genes_main_cryo <- fits_df_long_pluri_genes_main[grepl("Sydney", fits_df_long_pluri_genes_main$Location),]
fits_df_long_pluri_genes_main_cryo$Location <- factor((gsub("Sydney_cryopreserved", "Cryopreserved", fits_df_long_pluri_genes_main_cryo$Location) %>% gsub("Sydney", "Fresh",.)), levels = c("Cryopreserved", "Fresh"))


pPluri_Genes_Cont_cryo_main <- ggplot() +
						geom_bar(data = fits_df_long_pluri_genes_main_avg_cryo, aes(Location, Variance_Explained*100, fill = Covariate), position = "dodge", stat = "identity", width  = 0.75) +
						geom_point(data = fits_df_long_pluri_genes_main_cryo, aes(Location, Variance_Explained*100, fill = Covariate), size = 0.5, position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0)) +
						theme_classic() +
						facet_grid(Gene_ID ~ .) +
						coord_flip() +
						scale_fill_manual(values = variable_colors) +
						ylab("Percent Gene Expression\nVariance Explained") +
						theme(axis.title.y=element_blank(),
						legend.position = "none") 
save_figs(pPluri_Genes_Cont_cryo_main, paste0(outdir, "Pluripotent_Gene_Variable_Contributions_main_cryo"), width = 20, height = 20)


