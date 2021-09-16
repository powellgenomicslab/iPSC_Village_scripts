library(tidyverse)
library(ggplot2)
library(ggrepel)



dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
datadir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/nb_variance_partitioning_w_cryo/gene_separated/"
outdir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/nb_variance_partitioning_w_cryo/combined/"
dir.create(outdir)




save_figs <- function(plot, basename, width = 17, height = 17, units = "cm"){
    ggsave(plot, filename = paste0(basename,".png"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".pdf"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".eps"), height = height, width = width, units = units)
}



##### Set up colors #####
variable_colors <- c(Village = "#A2B0D0", Replicate = "#64A66B", Line = "#68319B") 
line_colors <- c(FSA0006 = "#F79E29", MBE1006 = "#9B2C99", TOB0421 = "#35369C")




##### Add gene IDs for easy identification downstream #####
GeneConversion1 <- read_delim("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/DRENEA_1/outs/filtered_feature_bc_matrix/features.tsv.gz", col_names = F, delim = "\t")
GeneConversion2 <- read_delim("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/Village_B_1_week/outs/filtered_feature_bc_matrix/features.tsv.gz", col_names = F, delim = "\t")

GeneConversion <- unique(rbind(GeneConversion1, GeneConversion2))
GeneConversion <- GeneConversion[!duplicated(GeneConversion$X1),]
GeneConversion$X3 <- NULL
colnames(GeneConversion) <- c("ENSG_ID", "Gene_ID")




if (!file.exists(paste0(outdir,"nb_variance_partitionin_df.tsv"))){
	##### Get list of files #####
	fits_files <- list.files(datadir, pattern = "_variable_variance_partition_coefficients.rds")
	names(fits_files) <- gsub("_variable_variance_partition_coefficients.rds", "", fits_files)


	fits_files_list <- list()
	for (location in c("Brisbane", "Melbourne" ,"Sydney_E", "Sydney_Cryopreserved")){
		fits_files_list[[gsub("_E", "", location)]] <- fits_files[which(grepl(location, fits_files))]
	}


	##### Read in files #####
	fits_list <- lapply(fits_files_list, function(x){
		lapply(x, function(y){
			readRDS(paste0(datadir,y))
		})
	})
	

	fits <- lapply(names(fits_list), function(x){
		temp <- lapply(names(fits_list[[x]]), function(y){
			fits_list[[x]][[y]]$Location <- y
			return(fits_list[[x]][[y]])
		})
		binded <- do.call(rbind, temp)
		return(binded)
	})

	##### Combine fits Results #####
	fits_df <- do.call(rbind, fits)

	fits_df$Location <- gsub("_ENSG.+", "", fits_df$Location)

	write_delim(fits_df, paste0(outdir,"nb_variance_partitionin_df.tsv"), delim = "\t")

} else {
	fits_df <- read_delim(paste0(outdir,"nb_variance_partitionin_df.tsv"), sep = "\t")
}


##### Make long dataframe for plotting #####
fits_df_long <- pivot_longer(fits_df, cols = c("Line", "Village", "Replicate"), names_to = "Covariate", values_to = "Variance_Explained")
fits_df_long$Variance_Explained <- round(fits_df_long$Variance_Explained,6)

fits_df_long <- left_join(fits_df_long, GeneConversion, by = c("Gene" = "ENSG_ID"))


##### Make some figures!!! #####
### FRESH ###
pTotal_Cont <- ggplot(fits_df_long[which(fits_df_long$Location != "Sydney_Cryopreserved"),], aes(Variance_Explained*100, color = Covariate)) +
	geom_density() +
	facet_grid(vars(Location)) +
	theme_classic()  +
	# scale_y_continuous(trans = "log10") +
	scale_x_continuous(trans = "log10") +
	scale_color_manual(values = variable_colors) +
	geom_vline(xintercept = 1, linetype="dashed")

save_figs(pTotal_Cont, paste0(outdir, "Total_Contribution_Histogram_cov"))


table(subset(fits_df_long[which(fits_df_long$Location != "Sydney_Cryopreserved"),], (Variance_Explained > 0.90))$Gene_ID)
subset(fits_df_long[which(fits_df_long$Location != "Sydney_Cryopreserved"),], (Variance_Explained > 0.90))
fits_df_long[grepl("XIST", fits_df_long$Gene_ID),]
fits_df_long[grepl("MT-ATP6", fits_df_long$Gene_ID),]

table(subset(fits_df_long[which(fits_df_long$Location != "Sydney_Cryopreserved" & fits_df_long$Covariate == "Village"),], (Variance_Explained > 0.46))$Gene_ID)



pVar_Explained_box <- ggplot(fits_df_long[which(fits_df_long$Location != "Sydney_Cryopreserved"),], aes(x = Covariate, y = Variance_Explained*100, color = Covariate)) +
	geom_boxplot(outlier.size = 0.5) +
	facet_grid(vars(Location)) +
	scale_color_manual(values = variable_colors) +
	theme_classic() +
	ylab("Percent Gene Expression\nVariance Explained") +
	geom_text_repel(data = subset(fits_df_long[which(fits_df_long$Location != "Sydney_Cryopreserved"),], 
							((Gene_ID %in% c("CHCHD2", "EIF1AY") & Covariate == "Line") | 
							(Gene_ID %in% c("SARS2") & Covariate == "Replicate") | 
							(Gene_ID %in% c("MT-ATP6", "MT-ND1") & Covariate == "Village"))), 
							aes(label = Gene_ID), box.padding =1, direction = "x")

save_figs(pVar_Explained_box, paste0(outdir, "Total_Contribution_Boxplot_cov"), width = 18)

	

fits_df_long[which(fits_df_long$Variance_Explained > 0.5 & fits_df_long$Covariate == "Village"),]
table(fits_df_long[which(fits_df_long$Variance_Explained > 0.4 & fits_df_long$Covariate == "Village"),]$Gene)
table(fits_df_long[which(fits_df_long$Variance_Explained > 0.5 & fits_df_long$Covariate == "Line"),]$Gene)
head(fits_df_long[which(fits_df_long$Covariate == "Village"),][rev(order(na.omit(fits_df_long[which(fits_df_long$Covariate == "Village"),"Variance_Explained"]))),])




### Variance explained correlation between sites ###
fits_df_joined <- full_join(fits_df_long, fits_df_long, by = c("Gene", "Covariate", "Gene_ID"))
fits_df_joined <- fits_df_joined[which(fits_df_joined$Location.x != fits_df_joined$Location.y),]


### identify the combinatinons to use
pairs_df <- as.data.frame(t(combn(unique(fits_df_joined$Location.x), 2, simplify = TRUE)))
colnames(pairs_df) <- c("Location.x", "Location.y")

fits_df_joined <- left_join(pairs_df, fits_df_joined)
fits_df_joined <- na.omit(fits_df_joined)


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

Rsquared <- na.omit(Rsquared)



pCorr_Location_Point_Line <- ggplot(fits_df_joined[which(fits_df_joined$Covariate == "Line" & fits_df_joined$Location.x != "Sydney_Cryopreserved" & fits_df_joined$Location.y != "Sydney_Cryopreserved"),], aes(Variance_Explained.x*100, Variance_Explained.y*100, color = Covariate)) +
	geom_point(size = 0.5, alpha = 0.25) +
	facet_grid(rows = vars(Location.y), cols = vars(Location.x)) +
	scale_color_manual(values = variable_colors[["Line"]]) +
	theme_classic() +
	ylim(0,110)+
	geom_text(x = 0, y = 110, aes(label = paste0("R^2 = ", Rsquared)), data = Rsquared[which(Rsquared$Covariate == "Line" & Rsquared$Location.x != "Sydney_Cryopreserved" & Rsquared$Location.y != "Sydney_Cryopreserved"),], hjust = 0) +
	xlab("Percent Gene Expression Variance Explained") +
	ylab("Percent Gene Expression Variance Explained") +
	geom_text_repel(data = subset(fits_df_joined[which(fits_df_joined$Location.x != "Sydney_Cryopreserved" & fits_df_joined$Location.y != "Sydney_Cryopreserved"),], 
							((Gene_ID %in% c("CHCHD2", "EIF1AY") & Covariate == "Line"))), 
							aes(label = Gene_ID), box.padding =0.5)
save_figs(pCorr_Location_Point_Line, paste0(outdir, "Correlation_Point_Line"), width = 15, height = 12)


pCorr_Location_Point_Village <- ggplot(fits_df_joined[which(fits_df_joined$Covariate == "Village" & fits_df_joined$Location.x != "Sydney_Cryopreserved" & fits_df_joined$Location.y != "Sydney_Cryopreserved"),], aes(Variance_Explained.x*100, Variance_Explained.y*100, color = Covariate)) +
	geom_point(size = 0.5, alpha = 0.25) +
	facet_grid(rows = vars(Location.y), cols = vars(Location.x)) +
	scale_color_manual(values = variable_colors[["Village"]]) +
	theme_classic() +
	ylim(0,70)+
	geom_text(x = 0, y = 70, aes(label = paste0("R^2 = ", Rsquared)), data = Rsquared[which(Rsquared$Covariate == "Village" & Rsquared$Location.x != "Sydney_Cryopreserved" & Rsquared$Location.y != "Sydney_Cryopreserved"),], hjust = 0) +
	xlab("Percent Gene Expression Variance Explained") +
	ylab("Percent Gene Expression Variance Explained") +
	geom_text_repel(data = subset(fits_df_joined[which(fits_df_joined$Location.x != "Sydney_Cryopreserved" & fits_df_joined$Location.y != "Sydney_Cryopreserved"),], 
							(Gene_ID %in% c("MT-ATP6", "MT-ND1") & Covariate == "Village")), 
							aes(label = Gene_ID), box.padding =1)
save_figs(pCorr_Location_Point_Village, paste0(outdir, "Correlation_Point_Village"), width = 15, height = 12)






### CRYOPRESERVED ###
fits_df_long_cryo <- fits_df_long[grepl("Sydney", fits_df_long$Location),]
fits_df_long_cryo$Location <- gsub("Sydney_Cryopreserved", "Cryopreserved", fits_df_long_cryo$Location) %>% gsub("Sydney", "Fresh", .)


pTotal_Cont_Cryo <- ggplot(fits_df_long_cryo aes(Variance_Explained*100, color = Covariate)) +
	geom_density() +
	facet_grid(vars(Location)) +
	theme_classic()  +
	# scale_y_continuous(trans = "log10") +
	scale_x_continuous(trans = "log10") +
	scale_color_manual(values = variable_colors) +
	geom_vline(xintercept = 1, linetype="dashed")

save_figs(pTotal_Cont_Cryo, paste0(outdir, "Total_Contribution_Histogram_cov_cryo"))



pVar_Explained_box_Cryo <- ggplot(fits_df_long_cryo, aes(x = Covariate, y = Variance_Explained*100, color = Covariate)) +
	geom_boxplot(outlier.size = 0.5) +
	facet_grid(vars(Location)) +
	scale_color_manual(values = variable_colors) +
	theme_classic() +
	ylab("Percent Gene Expression\nVariance Explained") +
	geom_text_repel(data = subset(fits_df_long_cryo, 
							((Gene_ID %in% c("CHCHD2", "ZNF71", "TBL1Y") & Covariate == "Line") | 
							(Gene_ID %in% c("MT-CO3") & Covariate == "Replicate") | 
							(Gene_ID %in% c("MT-CO3", "MT-ATP6") & Covariate == "Village"))), 
							aes(label = Gene_ID), box.padding =0.7, direction = "x")

save_figs(pVar_Explained_box_Cryo, paste0(outdir, "Total_Contribution_Boxplot_cov_cryo"), width = 18)

	
table(subset(fits_df_long_cryo, (Variance_Explained > 0.90))$Gene_ID)
subset(fits_df_long_cryo, (Variance_Explained > 0.90))
fits_df_long_cryo[which(fits_df_long_cryo$Gene_ID == "ZNF71"),]
subset(fits_df_long_cryo[which(fits_df_long_cryo$Covariate == "Replicate"),], Variance_Explained > 0.2)
table(subset(fits_df_long_cryo[which(fits_df_long_cryo$Covariate == "Village"),], Variance_Explained > 0.5)$Gene_ID)



### Variance explained correlation between sites ###
fits_df_joined_cryo <- full_join(fits_df_long_cryo, fits_df_long_cryo, by = c("Gene", "Covariate", "Gene_ID"))
fits_df_joined_cryo <- fits_df_joined_cryo[which(fits_df_joined_cryo$Location.x != fits_df_joined_cryo$Location.y),]


### identify the combinatinons to use
pairs_df <- as.data.frame(t(combn(unique(fits_df_joined_cryo$Location.x), 2, simplify = TRUE)))
colnames(pairs_df) <- c("Location.x", "Location.y")

fits_df_joined_cryo <- left_join(pairs_df, fits_df_joined_cryo)
fits_df_joined_cryo <- na.omit(fits_df_joined_cryo)


Rsquared_cryo <- unique(fits_df_joined_cryo[,c("Location.x", "Location.y", "Covariate")])
Rsquared_cryo$Rsquared <- NA
Rsquared_cryo$pearson <- NA


for (row in 1:nrow(Rsquared_cryo)){
	Rsquared_cryo$Rsquared[row] <- round(summary(lm(Variance_Explained.y ~ Variance_Explained.x, fits_df_joined_cryo[which(fits_df_joined_cryo$Location.x == Rsquared_cryo$Location.x[row] & 
																						fits_df_joined_cryo$Location.y == Rsquared_cryo$Location.y[row] & 
																						fits_df_joined_cryo$Covariate== Rsquared_cryo$Covariate[row]),]))$r.squared,2)
	Rsquared_cryo$pearson[row] <- round(cor(fits_df_joined_cryo[which(fits_df_joined_cryo$Location.x == Rsquared_cryo$Location.x[row] & 
																						fits_df_joined_cryo$Location.y == Rsquared_cryo$Location.y[row] & 
																						fits_df_joined_cryo$Covariate== Rsquared_cryo$Covariate[row]),]$Variance_Explained.y, 
					fits_df_joined_cryo[which(fits_df_joined_cryo$Location.x == Rsquared_cryo$Location.x[row] & 
																						fits_df_joined_cryo$Location.y == Rsquared_cryo$Location.y[row] & 
																						fits_df_joined_cryo$Covariate== Rsquared_cryo$Covariate[row]),]$Variance_Explained.x, use="complete.obs"), 2)
}

Rsquared_cryo <- na.omit(Rsquared_cryo)



pCorr_Location_Point_Line_cryo <- ggplot(fits_df_joined_cryo[which(fits_df_joined_cryo$Covariate == "Line"),], aes(Variance_Explained.x*100, Variance_Explained.y*100, color = Covariate)) +
	geom_point(size = 0.5, alpha = 0.25) +
	scale_color_manual(values = variable_colors[["Line"]]) +
	theme_classic() +
	ylim(0,110)+
	geom_text(x = 0, y = 110, aes(label = paste0("R^2 = ", Rsquared)), data = Rsquared_cryo[which(Rsquared_cryo$Covariate == "Line"),], hjust = 0) +
	xlab("Percent Gene Expression\nVariance Explained - Fresh") +
	ylab("Percent Gene Expression\nVariance Explained - Cryopreserved") +
	geom_text_repel(data = subset(fits_df_joined_cryo[which(fits_df_joined_cryo$Covariate == "Line"),], 
							(Gene_ID %in% c("CHCHD2", "ZNF71", "TBL1Y") & Covariate == "Line")), 
							aes(label = Gene_ID), box.padding = 1)

save_figs(pCorr_Location_Point_Line_cryo, paste0(outdir, "Correlation_Point_Line_Cryo"), width = 10, height = 8)


pCorr_Location_Point_Village_cryo <- ggplot(fits_df_joined_cryo[which(fits_df_joined_cryo$Covariate == "Village"),], aes(Variance_Explained.x*100, Variance_Explained.y*100, color = Covariate)) +
	geom_point(size = 0.5, alpha = 0.25) +
	scale_color_manual(values = variable_colors[["Village"]]) +
	theme_classic() +
	ylim(0,70)+
	geom_text(x = 0, y = 70, aes(label = paste0("R^2 = ", Rsquared)), data = Rsquared_cryo[which(Rsquared_cryo$Covariate == "Village"),], hjust = 0) +
	xlab("Percent Gene Expression\nVariance Explained - Fresh") +
	ylab("Percent Gene Expression\nVariance Explained - Cryopreserved") +
	geom_text_repel(data = subset(fits_df_joined_cryo[which(fits_df_joined_cryo$Covariate == "Village"),], 
							(Gene_ID %in% c("MT-CO3", "MT-ATP6") & Covariate == "Village")), 
							aes(label = Gene_ID), box.padding = 1)

save_figs(pCorr_Location_Point_Village_cryo, paste0(outdir, "Correlation_Point_Village_Cryo"), width = 10, height = 8)

