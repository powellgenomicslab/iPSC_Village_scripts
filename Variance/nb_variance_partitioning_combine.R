library(tidyverse)



dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
datadir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/nb_variance_partitioning/gene_separated/"
outdir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/nb_variance_partitioning/combined/"
dir.create(outdir)




save_figs <- function(plot, basename, width = 17, height = 17, units = "cm"){
    ggsave(plot, filename = paste0(basename,".png"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".pdf"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".eps"), height = height, width = width, units = units)
}



##### Set up colors #####
variable_colors <- c(Village = "#A2B0D0", Replicate = "#64A66B", Line = "#68319B") 
line_colors <- c(FSA0006 = "#F79E29", MBE1006 = "#9B2C99", TOB0421 = "#35369C")





if (!file.exists(paste0(outdir,"nb_variance_partitionin_df.tsv"))){
	##### Get list of files #####
	fits_files <- list.files(datadir, pattern = "_variable_variance_partition_coefficients.rds")
	names(fits_files) <- gsub("_variable_variance_partition_coefficients.rds", "", fits_files)


	##### Read in files #####
	fits_list <- lapply(fits_files, function(x){
		readRDS(paste0(datadir,x))
	})

	fits <- lapply(names(fits_list), function(x){
		temp <- lapply(names(fits_list[[x]]), function(y){
			fits_list[[x]][[y]]$Location <- y
			fits_list[[x]][[y]]$Gene <- x
			return(fits_list[[x]][[y]])
		})
		binded <- do.call(rbind, temp)
		return(binded)
	})

	##### Combine fits Results #####
	fits_df <- do.call(rbind, fits)

	write_delim(fits_df, paste0(outdir,"nb_variance_partitionin_df.tsv"), delim = "\t")

} else {
	fits_df <- read_delim(paste0(outdir,"nb_variance_partitionin_df.tsv"), sep = "\t")
}


##### Make long dataframe for plotting #####
fits_df_long <- pivot_longer(fits_df, cols = c("Line", "Village", "Replicate"), names_to = "Covariate", values_to = "Variance_Explained")
fits_df_long$Variance_Explained <- round(fits_df_long$Variance_Explained,6)



##### Make some figures!!! #####
pTotal_Cont <- ggplot(fits_df_long, aes(Variance_Explained, color = Covariate)) +
	geom_density() +
	facet_grid(vars(Location)) +
	theme_classic()  +
	# scale_y_continuous(trans = "log10") +
	scale_x_continuous(trans = "log10") +
	scale_color_manual(values = variable_colors) +
	geom_vline(xintercept = 0.01, linetype="dashed")

save_figs(pTotal_Cont, paste0(outdir, "Total_Contribution_Histogram_cov"))


pVar_Explained_box <- ggplot(fits_df_long, aes(x = Covariate, y = Variance_Explained, color = Covariate)) +
	geom_boxplot(outlier.size = 0.5) +
	facet_grid(vars(Location)) +
	scale_color_manual(values = variable_colors) +
	theme_classic()

save_figs(pVar_Explained_box, paste0(outdir, "Total_Contribution_Boxplot_cov"))

	
fits_df_long[which(fits_df_long$Variance_Explained > 0.5 & fits_df_long$Covariate == "Village"),]
table(fits_df_long[which(fits_df_long$Variance_Explained > 0.4 & fits_df_long$Covariate == "Village"),]$Gene)
table(fits_df_long[which(fits_df_long$Variance_Explained > 0.5 & fits_df_long$Covariate == "Line"),]$Gene)
head(fits_df_long[which(fits_df_long$Covariate == "Village"),][rev(order(na.omit(fits_df_long[which(fits_df_long$Covariate == "Village"),"Variance_Explained"]))),])



### Variance explained correlation between sites ###
fits_df_joined <- full_join(fits_df_long, fits_df_long, by = c("Gene", "Covariate"))
fits_df_joined <- fits_df_joined[which(fits_df_joined$Location.x != fits_df_joined$Location.y),]

Rsquared <- unique(fits_df_joined[,c("Location.x", "Location.y", "Covariate")])
Rsquared$Rsquared <- NA
Rsquared$pearson <- NA

for (row in 1:nrow(Rsquared)){
	Rsquared$Rsquared[row] <- round(summary(lm(Variance_Explained.y ~ Variance_Explained.x, fits_df_joined[which(fits_df_joined$Location.x == Rsquared$Location.x[row] & 
																						fits_df_joined$Location.y == Rsquared$Location.y[row] & 
																						fits_df_joined$Covariate== Rsquared$Covariate[row]),]))$r.squared,2)
	Rsquared$pearson[row] <- round(cor(fits_df_joined[which(fits_df_joined$Location.x == Rsquared$Location.x[row] & 
																						fits_df_joined$Location.y == Rsquared$Location.y[row] & 
																						fits_df_joined$Covariate== Rsquared$Covariate[row]),]$Variance_Explained.y, 
							fits_df_joined[which(fits_df_joined$Location.x == Rsquared$Location.x[row] & 
																						fits_df_joined$Location.y == Rsquared$Location.y[row] & 
																						fits_df_joined$Covariate== Rsquared$Covariate[row]),]$Variance_Explained.x, use="complete.obs"), 2)
}




pCorr_Location_Point <- ggplot(fits_df_joined, aes(Variance_Explained.x, Variance_Explained.y, color = Covariate)) +
	geom_point(size = 0.5, alpha = 0.25) +
	facet_grid(rows = vars(Location.y), cols = vars(Location.x), switch="both") +
	theme_classic() +
	scale_color_manual(values = variable_colors) +
	xlab("Contribution to Gene Variation") +
	ylab("Contribution to Gene Variation") +
	ylim(0,1.2) +
	geom_text(x = 0, y = 1, aes(label = paste0("R^2 = ", Rsquared)), data = Rsquared[which(Rsquared$Covariate == "Village"),], hjust = 0)+
	geom_text(x = 0, y = 1.2, aes(label = paste0("R^2 = ", Rsquared)), data = Rsquared[which(Rsquared$Covariate == "Line"),], hjust = 0) +
	geom_text(x = 0, y = 1.1, aes(label = paste0("R^2 = ", Rsquared)), data = Rsquared[which(Rsquared$Covariate == "Replicate"),], hjust = 0)
save_figs(pCorr_Location_Point, paste0(outdir, "Correlation_Point_All_Covariates"), width = 20, height = 18)



pCorr_Location_Point_Line <- ggplot(fits_df_joined[which(fits_df_joined$Covariate == "Line"),], aes(Variance_Explained.x, Variance_Explained.y, color = Covariate)) +
	geom_point(size = 0.5, alpha = 0.25) +
	facet_grid(rows = vars(Location.y), cols = vars(Location.x)) +
	scale_color_manual(values = variable_colors[["Line"]]) +
	theme_classic() +
	ylim(0,1.1) +
	geom_text(x = 0, y = 1.1, aes(label = paste0("R^2 = ", Rsquared)), data = Rsquared[which(Rsquared$Covariate == "Line"),], hjust = 0) 
save_figs(pCorr_Location_Point_Line, paste0(outdir, "Correlation_Point_Line"), width = 20, height = 18)


pCorr_Location_Point_Village <- ggplot(fits_df_joined[which(fits_df_joined$Covariate == "Village"),], aes(Variance_Explained.x, Variance_Explained.y, color = Covariate)) +
	geom_point(size = 0.5, alpha = 0.25) +
	facet_grid(rows = vars(Location.y), cols = vars(Location.x)) +
	scale_color_manual(values = variable_colors[["Village"]]) +
	theme_classic() +
	geom_text(x = 0, y = 0.65, aes(label = paste0("R^2 = ", Rsquared)), data = Rsquared[which(Rsquared$Covariate == "Village"),], hjust = 0)
save_figs(pCorr_Location_Point_Village, paste0(outdir, "Correlation_Point_Village"), width = 20, height = 18)


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


pPluri_Genes_Cont <- ggplot(fits_df_long_pluri_genes, aes(Location, Variance_Explained, fill = Covariate)) +
						geom_bar(position = "dodge", stat = "identity") +
						theme_classic() +
						facet_wrap(Gene_ID ~ ., nrow = 1) +
						scale_fill_manual(values = variable_colors) +
						theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1)) +
						ylab("Gene Expression Variance Explained")
save_figs(pPluri_Genes_Cont, paste0(outdir, "Pluripotent_Gene_Variable_Contributions"), width = 50, height = 10)


