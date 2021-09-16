library(tidyverse)



dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
datadir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/nb_partitioning_rep_separate/gene_separated/"
outdir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/nb_partitioning_rep_separate/combined/"
dir.create(outdir)




save_figs <- function(plot, basename, width = 17, height = 17, units = "cm"){
    ggsave(plot, filename = paste0(basename,".png"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".pdf"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".eps"), height = height, width = width, units = units)
}



##### Set up colors #####
variable_colors <- c(Village = "#A2B0D0", Replicate = "#64A66B", Line = "#68319B") 
line_colors <- c(FSA0006 = "#F79E29", MBE1006 = "#9B2C99", TOB0421 = "#35369C")

site_updates <- c("Brisbane" = "Site 1", "Melbourne" = "Site 2" ,"Sydney" = "Site 3")




if (!file.exists(paste0(outdir,"nb_variance_partitionin_df.tsv"))){
	##### Get list of files #####
	fits_files <- list.files(datadir, pattern = "_variable_variance_partition_coefficients.rds")
	names(fits_files) <- gsub("_variable_variance_partition_coefficients.rds", "", fits_files)


	fits_files_list <- list()
	for (rep in c("Brisbane1", "Brisbane2", "Brisbane3", "Melbourne1", "Melbourne2", "Melbourne3", "Sydney1", "Sydney2", "Sydney3", "Sydney_cryopreserved1", "Sydney_cryopreserved2", "Sydney_cryopreserved3")){
		fits_files_list[[rep]] <- fits_files[which(grepl(rep, fits_files))]
	}


	##### Read in files #####
	fits_list <- lapply(fits_files_list, function(x){
		lapply(x, function(y){
			readRDS(paste0(datadir,y))
		})
	})
	

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





##### Make long dataframe for plotting #####
fits_df_long <- pivot_longer(fits_df, cols = c("Line", "Village"), names_to = "Covariate", values_to = "Variance_Explained")
fits_df_long$Variance_Explained <- round(fits_df_long$Variance_Explained,6)
fits_df_long$Location <- gsub("\\d","",fits_df_long$Replicate)

for (location in names(site_updates)){
	fits_df_long$Location <- gsub(location, site_updates[location], fits_df_long$Location)
}


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


fits_df_long_pluri_genes <- left_join(pluri_genes, fits_df_long, by = c("ENSG_ID" = "Gene"))
fits_df_long_pluri_genes <- na.omit(fits_df_long_pluri_genes)

fits_df_long_pluri_genes_avg <- fits_df_long_pluri_genes %>% group_by(Gene_ID, ENSG_ID, Covariate, Location) %>% summarize(Variance_Explained = mean(Variance_Explained))


pPluri_Genes_Cont <- ggplot() +
						geom_bar(data = fits_df_long_pluri_genes_avg[which(fits_df_long_pluri_genes_avg$Location != "Site 3_cryopreserved"),], aes(Location, Variance_Explained*100, fill = Covariate), position = "dodge", stat = "identity") +
						geom_point(data = fits_df_long_pluri_genes[which(fits_df_long_pluri_genes$Location != "Site 3_cryopreserved"),], aes(Location, Variance_Explained*100, fill = Covariate), size = 0.5, position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0)) +
						theme_classic() +
						facet_wrap(Gene_ID ~ ., nrow = 3) +
						# scale_color_manual(values = variable_colors) +
						scale_fill_manual(values = variable_colors) +
						theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1)) +
						ylab("Percent Gene Expression Variance Explained") +
						theme(axis.title.x=element_blank(),
							legend.position="bottom")
save_figs(pPluri_Genes_Cont, paste0(outdir, "Pluripotent_Gene_Variable_Contributions"), width = 20, height = 12)


##### Only use main markers for main figure #####
fits_df_long_pluri_genes_main <- fits_df_long_pluri_genes[which(fits_df_long_pluri_genes$Gene_ID %in% c("NANOG", "POU5F1", "SOX2", "MYC")),]

fits_df_long_pluri_genes_main_avg <- fits_df_long_pluri_genes_main %>% group_by(Gene_ID, ENSG_ID, Covariate, Location) %>% summarize(Variance_Explained = mean(Variance_Explained))

fits_df_long_pluri_genes_main_avg$Location <- factor(fits_df_long_pluri_genes_main_avg$Location, levels = c("Site 3", "Site 2", "Site 1", "Site 3_cryopreserved"))


pPluri_Genes_Cont <- ggplot() +
						geom_bar(data = fits_df_long_pluri_genes_main_avg[which(fits_df_long_pluri_genes_main_avg$Location != "Site 3_cryopreserved"),], aes(Location, Variance_Explained*100, fill = Covariate), position = "dodge", stat = "identity", width  = 0.75) +
						geom_point(data = fits_df_long_pluri_genes_main[which(fits_df_long_pluri_genes_main$Location != "Site 3_cryopreserved"),], aes(Location, Variance_Explained*100, fill = Covariate), size = 0.5, position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0)) +
						theme_classic() +
						facet_grid(Gene_ID ~ .) +
						coord_flip() +
						scale_fill_manual(values = variable_colors) +
						ylab("Percent Gene Expression\nVariance Explained") +
						theme(axis.title.y=element_blank(),
						legend.position = "none") 
save_figs(pPluri_Genes_Cont, paste0(outdir, "Pluripotent_Gene_Variable_Contributions_main"), width = 8, height = 8)



##### Figures for cryopreserved figures #####
fits_df_long_pluri_genes_cryo <- fits_df_long_pluri_genes[grepl("Site 3", fits_df_long_pluri_genes$Location),]
fits_df_long_pluri_genes_cryo$Location <- factor((gsub("Site 3_cryopreserved", "Cryopreserved", fits_df_long_pluri_genes_cryo$Location) %>% gsub("Sydney", "Fresh",.)), levels = c("Fresh", "Cryopreserved"))


fits_df_long_pluri_genes_main_avg_cryo <- fits_df_long_pluri_genes_avg[grepl("Site 3", fits_df_long_pluri_genes_avg$Location),]
fits_df_long_pluri_genes_main_avg_cryo$Location <- factor((gsub("Site 3_cryopreserved", "Cryopreserved", fits_df_long_pluri_genes_main_avg_cryo$Location) %>% gsub("Site 3", "Fresh",.)), levels = c("Fresh", "Cryopreserved"))



pPluri_Genes_Cont_cryo <- ggplot() +
						geom_bar(data = fits_df_long_pluri_genes_main_avg_cryo, aes(Location, Variance_Explained*100, fill = Covariate), position = "dodge", stat = "identity") +
						geom_point(data = fits_df_long_pluri_genes_cryo, aes(Location, Variance_Explained*100, fill = Covariate), size = 0.5, position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0)) +
						theme_classic() +
						facet_wrap(Gene_ID ~ ., nrow = 3) +
						# scale_color_manual(values = variable_colors) +
						scale_fill_manual(values = variable_colors) +
						theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
							axis.title.x=element_blank(),
							legend.position="bottom") +
						ylab("Percent Gene Expression\nVariance Explained")
save_figs(pPluri_Genes_Cont_cryo, paste0(outdir, "Pluripotent_Gene_Variable_Contributions_cryo"), width = 18, height = 12)


### Main
fits_df_long_pluri_genes_main_avg_cryo <- fits_df_long_pluri_genes_main_avg[grepl("Site 2", fits_df_long_pluri_genes_main_avg$Location),]
fits_df_long_pluri_genes_main_avg_cryo$Location <- factor((gsub("Site 3_cryopreserved", "Cryopreserved", fits_df_long_pluri_genes_main_avg_cryo$Location) %>% gsub("Site 2", "Fresh",.)), levels = c("Cryopreserved", "Fresh"))

fits_df_long_pluri_genes_main_cryo <- fits_df_long_pluri_genes_main[grepl("Site 2", fits_df_long_pluri_genes_main$Location),]
fits_df_long_pluri_genes_main_cryo$Location <- factor((gsub("Site 3_cryopreserved", "Cryopreserved", fits_df_long_pluri_genes_main_cryo$Location) %>% gsub("Site 2", "Fresh",.)), levels = c("Cryopreserved", "Fresh"))


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
save_figs(pPluri_Genes_Cont_cryo_main, paste0(outdir, "Pluripotent_Gene_Variable_Contributions_main_cryo"), width = 10, height = 8)


