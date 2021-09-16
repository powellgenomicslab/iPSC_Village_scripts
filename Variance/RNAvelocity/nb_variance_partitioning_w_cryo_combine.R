library(tidyverse)
library(ggplot2)
library(ggrepel)
library(data.table)
library(Seurat)
library(viridis)
library(colorspace)
library(RColorBrewer)
library(ggsignif)


dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
seurat_dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/RNAvelocity/nb_partitioning/data/"
datadir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/RNAvelocity/nb_partitioning/gene_separated/"
outdir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/RNAvelocity/nb_partitioning/combined/"
dir.create(outdir)




save_figs <- function(plot, basename, width = 17, height = 17, units = "cm"){
    ggsave(plot, filename = paste0(basename,".png"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".pdf"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".eps"), height = height, width = width, units = units)
}



##### Set up colors #####
variable_colors <- c(Village = "#A2B0D0", Replicate = "#64A66B", Line = "#68319B") 
new_variable_colors <- c(Line = "#E5716F", Replicate = "#81BD9B", Village = "#7698AD")
line_colors <- c(FSA0006 = "#F79E29", MBE1006 = "#9B2C99", TOB0421 = "#35369C")

site_updates <- c("Brisbane" = "Site 1", "Melbourne" = "Site 2" ,"Sydney" = "Site 3")


##### Add gene IDs for easy identification downstream #####
GeneConversion1 <- read_delim("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/DRENEA_1/outs/filtered_feature_bc_matrix/features.tsv.gz", col_names = F, delim = "\t")
GeneConversion2 <- read_delim("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/Village_B_1_week/outs/filtered_feature_bc_matrix/features.tsv.gz", col_names = F, delim = "\t")

GeneConversion <- unique(rbind(GeneConversion1, GeneConversion2))
GeneConversion <- GeneConversion[!duplicated(GeneConversion$X1),]
GeneConversion$X3 <- NULL
colnames(GeneConversion) <- c("ENSG_ID", "Gene_ID")




if (!file.exists(paste0(outdir,"nb_variance_partitionin_df.tsv"))){
	##### Get list of files #####
	locations <- dir(seurat_dir)[!(dir(seurat_dir) %in% "logs")] %>% gsub("_SCT_seurat_1pct.rds", "", .)

	fits_files <- lapply(locations, function(x){
		tmp <- list.files(datadir, pattern = "_variable_variance_partition_coefficients.rds")
		names(tmp) <- gsub("_variable_variance_partition_coefficients.rds", "", tmp)
		return(tmp)
	})
	names(fits_files) <- locations


	##### Read in files #####
	fits_list <- lapply(fits_files, function(x){
		lapply(x, function(y){
			readRDS(paste0(datadir,y))
		})
	})
	

	fits <- lapply(names(fits_list), function(x){
		print(x)
		binded <- do.call(rbind, fits_list[[x]])
		return(binded)
	})

	##### Combine fits Results #####
	fits_df <- do.call(rbind, fits)

	fits_df$Location <- gsub("_ENSG.+", "",rownames(fits_df))

	write_delim(fits_df, paste0(outdir,"nb_variance_partitionin_df.tsv"), delim = "\t")

} else {
	fits_df <- read_delim(paste0(outdir,"nb_variance_partitionin_df.tsv"), delim = "\t")
}

##### Remove faulty model genes and locations #####
### Read in faulty log files ###
faulty <- fread(paste0(datadir,"faulty_runs.txt"), header = FALSE)


### Read in the seurat data objects
seurat_list <- list()
for (location in unique(fits_df$Location)){
	seurat_list[[location]] <- readRDS(paste0("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/RNAvelocity/nb_partitioning/data/",location,"_SCT_seurat_1pct.rds"))
}

### Make dataframe of location x number for faulty models ###
faulty_df <- data.table(Location = gsub("_nb\\.o.+", "", faulty$V1), Number = gsub(".+nb\\.o\\d+\\.", "", faulty$V1))

## Add gene id to the faulty dataframe
faulty_df$ENSG <- NA
for (row in 1:nrow(faulty_df)){
	print(row)
	faulty_df$ENSG[row] <- rownames(seurat_list[[faulty_df$Location[row]]])[as.numeric(faulty_df$Number[row])] 
}


### Remove the faulty genes ###
fits_df <- setDT(fits_df)[!faulty_df, on = c("Location", "Gene" = "ENSG")]



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

for (location in names(site_updates)){
	fits_df_long$Location <- gsub(location, site_updates[location], fits_df_long$Location)
}

fits_df_long$Location <- gsub("Cryo", "\nCryo", fits_df_long$Location)

fits_df_long$Quartile <- gsub(".+_Q", "Q", fits_df_long$Location)
fits_df_long$Location <- gsub("_Q\\d", "", fits_df_long$Location)
fits_df_long$Location <- gsub("_", "", fits_df_long$Location)

pVar_Explained_box <- ggplot(fits_df_long, aes(x = Covariate, y = Variance_Explained*100, color = Covariate)) +
	geom_boxplot(outlier.size = 0.25) +
	facet_grid(Location ~ Quartile) +
	scale_color_manual(values = variable_colors) +
	theme_classic() +
	ylab("Percent Gene Expression Variance Explained")  +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
save_figs(pVar_Explained_box, paste0(outdir, "Total_Contribution_Boxplot_cov"), width = 20, height = 12)



fits_df_long[which(fits_df_long$Variance_Explained > 0.5 & fits_df_long$Covariate == "Village"),]
table(fits_df_long[which(fits_df_long$Variance_Explained > 0.4 & fits_df_long$Covariate == "Village"),]$Gene)
table(fits_df_long[which(fits_df_long$Variance_Explained > 0.5 & fits_df_long$Covariate == "Line"),]$Gene)
head(fits_df_long[which(fits_df_long$Covariate == "Village"),][rev(order(na.omit(fits_df_long[which(fits_df_long$Covariate == "Village"),"Variance_Explained"]))),])




### Variance explained correlation between sites ###
fits_df_long$Location_Quartile <- paste0(fits_df_long$Location, "-", fits_df_long$Quartile)
fits_df_long <- data.table(fits_df_long)

fwrite(fits_df_long, sep = "\t", paste0(outdir, "var_explained_long.tsv"))

### Subset
fits_df_long_line <- fits_df_long[Covariate == "Line"]
fits_df_long_rep <- fits_df_long[Covariate == "Replicate"]
fits_df_long_village <- fits_df_long[Covariate == "Village"]


### identify the combinatinons to use
pairs_df <- as.data.table(t(combn(sort(unique(fits_df_long$Location_Quartile)), 2, simplify = TRUE)))
colnames(pairs_df) <- c("Location_Quartile.x", "Location_Quartile.y")


### Self Join ###
fits_df_long_line_joined <- merge(pairs_df, fits_df_long_line, by.x = "Location_Quartile.x", by.y = "Location_Quartile", allow.cartesian=TRUE)
fits_df_long_line_joined <- merge(fits_df_long_line_joined, fits_df_long_line, by.x = c("Location_Quartile.y", "Gene", "Covariate", "Gene_ID"), by.y = c("Location_Quartile", "Gene", "Covariate", "Gene_ID"), allow.cartesian=TRUE)
fits_df_long_line_joined <- na.omit(fits_df_long_line_joined)

fits_df_long_rep_joined <- merge(pairs_df, fits_df_long_rep, by.x = "Location_Quartile.x", by.y = "Location_Quartile", allow.cartesian=TRUE)
fits_df_long_rep_joined <- merge(fits_df_long_rep_joined, fits_df_long_rep, by.x = c("Location_Quartile.y", "Gene", "Covariate", "Gene_ID"), by.y = c("Location_Quartile", "Gene", "Covariate", "Gene_ID"), allow.cartesian=TRUE)
fits_df_long_rep_joined <- na.omit(fits_df_long_rep_joined)

fits_df_long_village_joined <- merge(pairs_df, fits_df_long_village, by.x = "Location_Quartile.x", by.y = "Location_Quartile", allow.cartesian=TRUE)
fits_df_long_village_joined <- merge(fits_df_long_village_joined, fits_df_long_village, by.x = c("Location_Quartile.y", "Gene", "Covariate", "Gene_ID"), by.y = c("Location_Quartile", "Gene", "Covariate", "Gene_ID"), allow.cartesian=TRUE)
fits_df_long_village_joined <- na.omit(fits_df_long_village_joined)


### Remove the ones that are merges between the same ###
fits_df_long_line_joined <- fits_df_long_line_joined[which(fits_df_long_line_joined$Location_Quartile.x != fits_df_long_line_joined$Location_Quartile.y),]
fits_df_long_rep_joined <- fits_df_long_rep_joined[which(fits_df_long_rep_joined$Location_Quartile.x != fits_df_long_rep_joined$Location_Quartile.y),]
fits_df_long_village_joined <- fits_df_long_village_joined[which(fits_df_long_village_joined$Location_Quartile.x != fits_df_long_village_joined$Location_Quartile.y),]


fits_df_long_joined_list <- list(fits_df_long_line_joined, fits_df_long_rep_joined, fits_df_long_village_joined)
names(fits_df_long_joined_list) <- c("Line", "Replicate", "Village")


# ################
# lapply(fits_df_long_joined_list, function(x){
# 	x$Location <- gsub("Site 3", "Site 4", x$Location) %>% gsub("Site 2", "Site 3", .) %>% gsub("Site 4", "Site 3", .) 
# 	x$Location_Quartile.x <- gsub("Site 3", "Site 4", x$Location_Quartile.x) %>% gsub("Site 2", "Site 3", .) %>% gsub("Site 4", "Site 3", .) 
# 	x$Location_Quartile.y <- gsub("Site 3", "Site 4", x$Location_Quartile.y) %>% gsub("Site 2", "Site 3", .) %>% gsub("Site 4", "Site 3", .) 
# 	return(x)
# })

# ################



saveRDS(fits_df_long_joined_list, paste0(outdir, "variance_explained_long_list.rds"))
fits_df_long_joined_list <- readRDS(paste0(outdir, "variance_explained_long_list.rds"))

Rsquared_list <- lapply(names(fits_df_long_joined_list), function(x){
	tmp <- pairs_df
	tmp$Covariate <- x
	tmp$Rsquared <- NA
	return(data.table(tmp))
})
names(Rsquared_list) <- names(fits_df_long_joined_list)



Rsquared_list <- lapply(names(Rsquared_list), function(x){
	print(x)
	for (row in 1:nrow(Rsquared_list[[x]])){
		print(row)
		Rsquared_list[[x]]$Rsquared[row] <- round(summary(lm(Variance_Explained.y ~ Variance_Explained.x, fits_df_long_joined_list[[x]][Location_Quartile.x == Rsquared_list[[x]]$Location_Quartile.x[row] & Location_Quartile.y == Rsquared_list[[x]]$Location_Quartile.y[row],]))$r.squared,2)
		# Rsquared_list[[x]]$pearson[row] <- round(cor(fits_df_long_joined_list[[x]][which(fits_df_long_joined_list[[x]]$Location_Quartile.x == Rsquared_list[[x]]$Location_Quartile.x[row] & 
		# 																					fits_df_long_joined_list[[x]]$Location_Quartile.y == Rsquared_list[[x]]$Location_Quartile.y[row]),]$Variance_Explained.y, 
		# 						fits_df_long_joined_list[[x]][which(fits_df_long_joined_list[[x]]$Location_Quartile.x == Rsquared_list[[x]]$Location_Quartile.x[row] & 
		# 																					fits_df_long_joined_list[[x]]$Location_Quartile.y == Rsquared_list[[x]]$Location_Quartile.y[row]),]$Variance_Explained.x, use="complete.obs"), 2)
	}

	Rsquared_list[[x]] <- na.omit(Rsquared_list[[x]])
	return(Rsquared_list[[x]])
})
names(Rsquared_list) <- names(fits_df_long_joined_list)

saveRDS(Rsquared_list, paste0(outdir, "R_squared_list.rds"))




fits_df_long_joined_list <- readRDS(paste0(outdir, "variance_explained_long_list.rds"))
Rsquared_list <- readRDS(paste0(outdir, "R_squared_list.rds"))


Rsquared_list <- lapply(Rsquared_list, function(x){
	x$Location_Quartile.x <- gsub("Site 3", "Site 4", x$Location_Quartile.x) %>% gsub("Site 2", "Site 3", .) %>% gsub("Site 4", "Site 2", .) %>% gsub("Cryo", "\nCryo", .)
	x$Location_Quartile.y <- gsub("Site 3", "Site 4", x$Location_Quartile.y) %>% gsub("Site 2", "Site 3", .) %>% gsub("Site 4", "Site 2", .) %>% gsub("Cryo", "\nCryo", .)
	return(x)
})


pCorr_tile <- ggplot(Rsquared_list[["Line"]], aes(Location_Quartile.x, Location_Quartile.y, fill = Rsquared)) +
				geom_tile() +
				theme_classic() +
				theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
				scale_fill_continuous_sequential(palette = "Blues")
ggsave(pCorr_tile, filename = paste0(outdir,"Rsquared_heatmap.png"), height = 5.85)



##### Make lists that are only same site #####
Rsquared_list_subset <- lapply(Rsquared_list, function(x){
	tmp <- list()
	for (location in unique(gsub("-.+", "", unique(c(x$Location_Quartile.x, x$Location_Quartile.y))))){
		tmp[[location]] <- x[grepl(paste0(location, "-"), x$Location_Quartile.x) & grepl(paste0(location, "-"), Location_Quartile.y)]
		tmp[[location]]$Location <- location
		tmp[[location]]$Quintile.x <- gsub(".+-", "", tmp[[location]]$Location_Quartile.x)
		tmp[[location]]$Quintile.y <- gsub(".+-", "", tmp[[location]]$Location_Quartile.y)
	}
	tmp2 <- do.call(rbind, tmp)
	return(tmp2)
})


pCorr_tile <- ggplot(Rsquared_list_subset[["Line"]], aes(Quintile.x, Quintile.y, fill = Rsquared)) +
				geom_tile() +
				theme_classic() +
				facet_wrap(vars(Location), nrow = 2) +
				theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
				scale_fill_continuous_sequential(palette = "Purples") +
				xlab("Quintile") +
				ylab("Quintile") +
				labs(fill = "Line\nR^2")
ggsave(pCorr_tile, filename = paste0(outdir,"Rsquared_heatmap_location_facet.png"), height = 3.25, width = 3.5)
ggsave(pCorr_tile, filename = paste0(outdir,"Rsquared_heatmap_location_facet.pdf"), height = 3.25, width = 3.5)


pCorr_Village<- ggplot(Rsquared_list_subset[["Village"]], aes(Quintile.x, Quintile.y, fill = Rsquared)) +
				geom_tile() +
				theme_classic() +
				facet_wrap(vars(Location), nrow = 2) +
				theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
				scale_fill_continuous_sequential(palette = "Blues") +				
				xlab("Quintile") +
				ylab("Quintile") +
				labs(fill = "Village\nR^2")
ggsave(pCorr_Village, filename = paste0(outdir,"Rsquared_heatmap_village_facet.png"), height = 3.25, width = 3.5)
ggsave(pCorr_Village, filename = paste0(outdir,"Rsquared_heatmap_village_facet.pdf"), height = 3.25, width = 3.5)





##### Make lists that are only same quintile #####
fits_df_long <- fread(sep = "\t", paste0(outdir, "var_explained_long.tsv"))


Rsquared_list_subset_quint <- lapply(Rsquared_list, function(x){
	tmp <- list()
	for (quintile in unique(fits_df_long$Quartile)){
		tmp[[quintile]] <- x[grepl(paste0("-", quintile), x$Location_Quartile.x) & grepl(paste0("-", quintile), x$Location_Quartile.y)]
		tmp[[quintile]]$Quintile <- quintile
		tmp[[quintile]]$Location.x <- gsub("-.+", "", tmp[[quintile]]$Location_Quartile.x)
		tmp[[quintile]]$Location.y <- gsub("-.+", "", tmp[[quintile]]$Location_Quartile.y)
	}
	tmp2 <- do.call(rbind, tmp)
	return(tmp2)
})

bluecols <- colorRampPalette(brewer.pal(9, "Blues"))
purplecols <- colorRampPalette(brewer.pal(9, "Purples"))

pCorr_Line_quant <- ggplot(Rsquared_list_subset_quint[["Line"]], aes(Location.x, Location.y, fill = Rsquared)) +
				geom_tile() +
				theme_classic() +
				facet_wrap(vars(Quintile), nrow = 1) +
				theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
				scale_fill_gradientn(colours = purplecols(100), limits=c(0, 1))
				# scale_fill_continuous_sequential(palette = "Purples")
ggsave(pCorr_Line_quant, filename = paste0(outdir,"Rsquared_heatmap_line_facet_quint.png"), height = 2.5, width = 6.5)

pCorr_Village_quant <- ggplot(Rsquared_list_subset_quint[["Village"]], aes(Location.x, Location.y, fill = Rsquared)) +
				geom_tile() +
				theme_classic() +
				facet_wrap(vars(Quintile), nrow = 1) +
				theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
				scale_fill_gradientn(colours = bluecols(100), limits=c(0, 1))
				# scale_fill_continuous_sequential(palette = "Blues", begin = 0, end = 1)
ggsave(pCorr_Village_quant, filename = paste0(outdir,"Rsquared_heatmap_village_facet_quint.png"), height = 2.5, width = 6.5)



pCorr_Location_Point_Line <- ggplot(fits_df_long_joined_list[["Line"]], aes(Variance_Explained.x*100, Variance_Explained.y*100, color = Covariate)) +
	geom_point(size = 0.5, alpha = 0.5) +
	facet_grid(rows = vars(Location.y), cols = vars(Location.x)) +
	scale_color_manual(values = variable_colors[["Line"]]) +
	theme_classic() +
	ylim(0,110)+
	# geom_text(size = 2.5, x = 0, y = 110, aes(label = paste0("R^2 = ", Rsquared)), data = Rsquared_list[["Line"]], hjust = 0) +
	xlab("Percent Gene Expression Variance Explained") +
	ylab("Percent Gene Expression Variance Explained") 
save_figs(pCorr_Location_Point_Line, paste0(outdir, "Correlation_Point_Line"), width = 50, height = 5*7.25)


pCorr_Location_Point_Village <- ggplot(fits_df_joined[which(fits_df_joined$Covariate == "Village" & fits_df_joined$Location.x != "Site 2_Cryopreserved" & fits_df_joined$Location.y != "Site 2_Cryopreserved"),], aes(Variance_Explained.x*100, Variance_Explained.y*100, color = Covariate)) +
	geom_point(size = 0.5, alpha = 0.5) +
	facet_grid(rows = vars(Location.y), cols = vars(Location.x)) +
	scale_color_manual(values = variable_colors[["Village"]]) +
	theme_classic() +
	ylim(0,70)+
	geom_text(size = 2.5, x = 0, y = 70, aes(label = paste0("R^2 = ", Rsquared)), data = Rsquared[which(Rsquared$Covariate == "Village" & Rsquared$Location.x != "Site 2_Cryopreserved" & Rsquared$Location.y != "Site 2_Cryopreserved"),], hjust = 0) +
	xlab("Percent Gene Expression\nVariance Explained") +
	ylab("Percent Gene Expression\nVariance Explained") +
	geom_text_repel(size = 2.5, data = subset(fits_df_joined[which(fits_df_joined$Location.x != "Site 2_Cryopreserved" & fits_df_joined$Location.y != "Site 2_Cryopreserved"),], 
							(Gene_ID %in% c("MT-ATP6", "MT-ND1") & Covariate == "Village")), 
							aes(label = Gene_ID), box.padding =1)
save_figs(pCorr_Location_Point_Village, paste0(outdir, "Correlation_Point_Village"), width = 10, height = 7.25)




##### Make a point figure instead of heatmap for Correlation between sites for each quintile
Rsquared_list_subset <- do.call(rbind, Rsquared_list_subset_quint)

pR2_point <- ggplot(Rsquared_list_subset[Covariate != "Replicate"], aes(Covariate, Rsquared, color = Quintile)) +
				geom_jitter(width = 0.1, size = 0.75) +
				theme_classic() +
				ylab("R^2") +
				scale_color_viridis(discrete = TRUE) +
				geom_signif(comparisons = list(c("Line", "Village")), 
								map_signif_level=TRUE, y = 1,
								test = "t.test", test.args=list(alternative = "two.sided", var.equal = FALSE, paired=TRUE))

ggsave(pR2_point, filename = paste0(outdir, "R2_scatter.png"), width = 2.5, height = 3)
ggsave(pR2_point, filename = paste0(outdir, "R2_scatter.pdf"), width = 2.5, height = 3)


t.test(Rsquared_list_subset[Covariate == "Village"]$Rsquared, Rsquared_list_subset[Covariate == "Line"]$Rsquared, alternative = "two.sided", var.equal = FALSE, paired=TRUE)

