library(tidyverse)
library(data.table)
library(ggrepel)



dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
datadir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/nb_partitioning_village_separate/gene_separated/"
outdir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/nb_partitioning_village_separate/combined/"
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


##### Add gene IDs for easy identification downstream #####
GeneConversion1 <- read_delim("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/DRENEA_1/outs/filtered_feature_bc_matrix/features.tsv.gz", col_names = F, delim = "\t")
GeneConversion2 <- read_delim("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/Village_B_1_week/outs/filtered_feature_bc_matrix/features.tsv.gz", col_names = F, delim = "\t")

GeneConversion <- unique(rbind(GeneConversion1, GeneConversion2))
GeneConversion <- GeneConversion[!duplicated(GeneConversion$X1),]
GeneConversion$X3 <- NULL
colnames(GeneConversion) <- c("ENSG_ID", "Gene_ID")



if (!file.exists(paste0(outdir,"nb_variance_partitioning_df.tsv"))){
	##### Get list of files #####
	fits_files <- list.files(datadir, pattern = "_variable_variance_partition_coefficients.rds")
	names(fits_files) <- gsub("_variable_variance_partition_coefficients.rds", "", fits_files)


	fits_files_list <- list()
	fits_files_list[["Baseline"]][["Brisbane"]] <- fits_files[which(grepl("Baseline", fits_files) & grepl("Brisbane", fits_files))]
	fits_files_list[["Baseline"]][["Melbourne"]] <- fits_files[which(grepl("Baseline", fits_files) & grepl("Melbourne", fits_files))]
	fits_files_list[["Baseline"]][["Sydney_Fresh"]] <- fits_files[which(grepl("Baseline", fits_files) & grepl("Sydney_Fresh", fits_files))]
	fits_files_list[["Baseline"]][["Sydney_Cryopreserved"]] <- fits_files[which(grepl("Baseline", fits_files) & grepl("Sydney_Cryopreserved", fits_files))]
	fits_files_list[["Village"]][["Brisbane"]] <- fits_files[which(grepl("Village", fits_files) & grepl("Brisbane", fits_files))]
	fits_files_list[["Village"]][["Melbourne"]] <- fits_files[which(grepl("Village", fits_files) & grepl("Melbourne", fits_files))]
	fits_files_list[["Village"]][["Sydney_Fresh"]] <- fits_files[which(grepl("Village", fits_files) & grepl("Sydney_Fresh", fits_files))]
	fits_files_list[["Village"]][["Sydney_Cryopreserved"]] <- fits_files[which(grepl("Village", fits_files) & grepl("Sydney_Cryopreserved", fits_files))]

	##### Read in files #####
	fits_list <- lapply(fits_files_list, function(x){
		lapply(x, function(y){
			lapply(y, function(z){
				readRDS(paste0(datadir,z))
			})
		})
	})

	fits <- lapply(names(fits_list), function(x){
		temp <- lapply(names(fits_list[[x]]), function(y){
			temp2 <- lapply(names(fits_list[[x]][[y]]), function(z){
				fits_list[[x]][[y]][[z]]$Village <- x
				fits_list[[x]][[y]][[z]]$Location <- y
				return(fits_list[[x]][[y]][[z]])
			})
			bind1 <- do.call(rbind, temp2)
		})
		binded <- do.call(rbind, temp)
		return(binded)
	})

	##### Combine fits Results #####
	fits_df <- do.call(rbind, fits)

	write_delim(fits_df, paste0(outdir,"nb_variance_partitioning_df.tsv"), delim = "\t")

} else {
	fits_df <- fread(paste0(outdir,"nb_variance_partitioning_df.tsv"), sep = "\t")
}

fits_df$Location <- gsub("_Fresh", "", fits_df$Location)


##### Remove faulty model genes and locations #####
### Read in faulty log files ###
faulty <- fread(paste0(datadir,"faulty_runs.txt"), header = FALSE)


### Read in the seurat data objects
seurat_list <- list()
for (location in unique(fits_df$Location)){
	seurat_list[[location]] <- readRDS(paste0("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/nb_partitioning_cryo/data/",location,"_SCT_seurat_1pct.rds"))
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


for (location in names(site_updates)){
	fits_df$Location <- gsub(location, site_updates[location], fits_df$Location) %>% gsub("_Cryopreserved", " Cryopreserved", .) %>% gsub("_Fresh", "", .)
}




##### Make long dataframe for plotting #####
fits_df_long <- pivot_longer(fits_df, cols = c("Line", "Village", "Replicate"), names_to = "Covariate", values_to = "Variance_Explained")
fits_df_long$Variance_Explained <- round(fits_df_long$Variance_Explained,6)

fits_df_long <- left_join(fits_df_long, GeneConversion, by = c("Gene" = "ENSG_ID"))


##### Make long dataframe for plotting #####
fits_df_long <- pivot_longer(fits_df, cols = c("Line", "Replicate"), names_to = "Covariate", values_to = "Variance_Explained")
fits_df_long$Variance_Explained <- round(fits_df_long$Variance_Explained,6)

fits_df_long <- left_join(fits_df_long, GeneConversion, by = c("Gene" = "ENSG_ID"))



##### Make some figures!!! #####
pTotal_Cont <- ggplot(fits_df_long, aes(Variance_Explained, color = Covariate)) +
	geom_density() +
	facet_grid(Location ~ Village) +
	theme_classic()  +
	# scale_y_continuous(trans = "log10") +
	scale_x_continuous(trans = "log10") +
	scale_color_manual(values = variable_colors) +
	geom_vline(xintercept = 0.01, linetype="dashed")

save_figs(pTotal_Cont, paste0(outdir, "Total_Contribution_Histogram_cov"))


pVar_Explained_box <- ggplot(fits_df_long, aes(x = Covariate, y = Variance_Explained, color = Covariate)) +
	geom_boxplot(outlier.size = 0.5) +
	facet_grid(Location ~ Village) +
	scale_color_manual(values = variable_colors) +
	theme_classic()

save_figs(pVar_Explained_box, paste0(outdir, "Total_Contribution_Boxplot_cov"))




### Variance explained correlation between sites ###
fits_df_long$Location_Village <- paste0(fits_df_long$Location, ".", fits_df_long$Village)
fits_df_joined <- full_join(fits_df_long, fits_df_long, by = c("Gene", "Covariate", "Gene_ID"))

fits_df_joined <- fits_df_joined[which(fits_df_joined$Location_Village.x != fits_df_joined$Location_Village.y),]
fits_df_joined <- fits_df_joined[which(fits_df_joined$Location.x == fits_df_joined$Location.y),]
fits_df_joined <- fits_df_joined[which(fits_df_joined$Village.x == "Baseline"),]

fits_df_joined$Location.x <- factor(fits_df_joined$Location.x, levels = c("Site 1", "Site 2", "Site 3", "Site 3 Cryopreserved"))


Rsquared <- unique(fits_df_joined[,c("Location_Village.x", "Location_Village.y","Location.x", "Location.y", "Village.x", "Village.y", "Covariate")])
Rsquared$Rsquared <- NA
Rsquared$pearson <- NA

for (row in 1:nrow(Rsquared)){
	Rsquared$Rsquared[row] <- round(summary(lm(Variance_Explained.y ~ Variance_Explained.x, fits_df_joined[which(fits_df_joined$Location.x == Rsquared$Location.x[row] & 
																						fits_df_joined$Location.y == Rsquared$Location.y[row] & 
																						fits_df_joined$Village.x == Rsquared$Village.x[row] & 
																						fits_df_joined$Village.y == Rsquared$Village.y[row] & 
																						fits_df_joined$Covariate== Rsquared$Covariate[row]),]))$r.squared,2)
	Rsquared$pearson[row] <- round(cor(fits_df_joined[which(fits_df_joined$Location.x == Rsquared$Location.x[row] & 
																						fits_df_joined$Location.y == Rsquared$Location.y[row] & 
																						fits_df_joined$Village.x == Rsquared$Village.x[row] & 
																						fits_df_joined$Village.y == Rsquared$Village.y[row] & 
																						fits_df_joined$Covariate== Rsquared$Covariate[row]),]$Variance_Explained.y, 
							fits_df_joined[which(fits_df_joined$Location.x == Rsquared$Location.x[row] & 
																						fits_df_joined$Location.y == Rsquared$Location.y[row] & 
																						fits_df_joined$Village.x == Rsquared$Village.x[row] & 
																						fits_df_joined$Village.y == Rsquared$Village.y[row] & 
																						fits_df_joined$Covariate== Rsquared$Covariate[row]),]$Variance_Explained.x, use="complete.obs"), 2)
}


as.data.frame(Rsquared[which(Rsquared$Covariate == "Line"),])

Rsquared <- data.table(Rsquared)
fits_df_joined <- data.table(fits_df_joined)



pCorr_Location_Point <- ggplot(fits_df_joined, aes(Variance_Explained.x, Variance_Explained.y, color = Covariate)) +
	geom_point(size = 0.5, alpha = 0.5) +
	facet_grid(rows = vars(Location_Village.y), cols = vars(Location_Village.x), switch="both") +
	theme_classic() +
	scale_color_manual(values = variable_colors) +
	xlab("Contribution to Gene Variation") +
	ylab("Contribution to Gene Variation") +
	ylim(0,1.2) +
	geom_text(x = 0, y = 1.2, aes(label = paste0("R = ", pearson)), data = Rsquared[which(Rsquared$Covariate == "Line"),], hjust = 0) +
	geom_text(x = 0, y = 1.1, aes(label = paste0("R = ", pearson)), data = Rsquared[which(Rsquared$Covariate == "Replicate"),], hjust = 0)
save_figs(pCorr_Location_Point, paste0(outdir, "Correlation_Point_All_Covariates"), width = 20, height = 18)



pCorr_Location_Point_Line <- ggplot(fits_df_joined[Covariate == "Line" & Location.x != "Site 3 Cryopreserved"], aes(Variance_Explained.x, Variance_Explained.y, color = Covariate)) +
	geom_point(size = 0.5, alpha = 0.5) +
	facet_wrap(vars(Location.x), nrow = 1) +
	scale_color_manual(values = variable_colors[["Line"]]) +
	theme_classic() +
	ylim(0,1.05) +
	xlab("Baseline Variance Explained") +
	ylab("Village Variance\nExplained") +
	geom_text(x = 0, y = 1, aes(label = paste0("R = ", pearson)), data = Rsquared[Covariate == "Line" & Location.x != "Site 3 Cryopreserved"], hjust = 0, color = "black", size = 3) +
	theme(panel.spacing.x = unit(5, "mm")) +
	geom_text_repel(size = 3, data = subset(fits_df_joined[Covariate == "Line" & Location.x != "Site 3 Cryopreserved"],
							(Gene_ID %in% c("CHCHD2", "EIF1AY") & Covariate == "Line")), 
							aes(label = Gene_ID), box.padding = 1, color = "black")
save_figs(pCorr_Location_Point_Line, paste0(outdir, "Correlation_Point_Line_Time"), width = 14, height = 5)



##### Make Fig 1 subplot #####
pCorr_Location_Point_Line_site1 <- ggplot(fits_df_joined[Covariate == "Line" & Location.x == "Site 1" & Location.y == "Site 1"], aes(Variance_Explained.x*100, Variance_Explained.y*100, color = Covariate)) +
	geom_point(size = 0.5, alpha = 0.5) +
	scale_color_manual(values = variable_colors[["Line"]]) +
	theme_classic() +
	ylim(0,110)+
	theme(legend.position = "none")+
	geom_text(size = 4.25, x = 0, y = 105, aes(label = paste0("R = ", pearson)), data = Rsquared[Covariate == "Line" & Location.x == "Site 1" & Location.y == "Site 1",], hjust = 0, color = "black") +
	xlab("Baseline Percent Gene Expression\nVariance Explained") +
	ylab("Village Percent Gene Expression\nVariance Explained")
save_figs(pCorr_Location_Point_Line_site1, paste0(outdir, "Correlation_Point_Line_Time_site_1"), width = 8, height = 8)




fits_df_joined_cryo <- fits_df_joined[grepl("Site 3", Location.x)]
fits_df_joined_cryo$Time <- ifelse(fits_df_joined_cryo$Location.x == "Site 3 Cryopreserved", "Cryopreserved", "Fresh")

Rsquared_cryo <- Rsquared[grepl("Site 3", Location.x)]
Rsquared_cryo$Time <- ifelse(Rsquared_cryo$Location.x == "Site 3 Cryopreserved", "Cryopreserved", "Fresh")
Rsquared_cryo$Time <- factor(Rsquared_cryo$Time, levels = c("Fresh", "Cryopreserved"))


pCorr_Location_Point_Line_Cryo <- ggplot(fits_df_joined_cryo[Covariate == "Line"], aes(Variance_Explained.x, Variance_Explained.y, color = Covariate)) +
	geom_point(size = 0.5, alpha = 0.5) +
	facet_wrap(vars(Time), nrow = 1) +
	scale_color_manual(values = variable_colors[["Line"]]) +
	theme_classic() +
	ylim(0,1.05) +
	xlab("Baseline Variance Explained") +
	ylab("Village Variance\nExplained") +
	geom_text(x = 0, y = 1, aes(label = paste0("R = ", pearson)), data = Rsquared_cryo[Covariate == "Line"], hjust = 0, color = "black", size = 3) +
	theme(panel.spacing.x = unit(5, "mm")) +
	geom_text_repel(size = 3, data = subset(fits_df_joined_cryo[Covariate == "Line"], 
							(Gene_ID %in% c("CHCHD2", "EIF1AY") & Covariate == "Line")), 
							aes(label = Gene_ID), box.padding = 1, color = "black")
save_figs(pCorr_Location_Point_Line_Cryo, paste0(outdir, "Correlation_Point_Line_Time_cryopreserved"), width = 11, height = 5)





##### Filter to 10% and rerun #####
##### Read in Data ######
location_SCT <- readRDS(paste0(dir,"output/Variance_contributions_lm_scale_data/seurat_SCT_wRB_sydney_regress_all_covs.rds"))
cryo_SCT <- readRDS(paste0(dir,"output/Variance_contributions_lm_scale_data_Sydney/seurat_SCT_wRB_sydney_regress_all_covs.rds"))


### Update names to combine objects ###
names(location_SCT) <- c("Baseline", "Village")
location_SCT <- lapply(location_SCT, function(x){
	names(x) <- c("Brisbane", "Melbourne", "Sydney_Fresh")
	return(x)
})


cryo_SCT_updated <- list()
cryo_SCT_updated[["Baseline"]][["Sydney_Fresh"]] <- cryo_SCT[["Baseline"]]
cryo_SCT_updated[["Baseline"]][["Sydney_Cryopreserved"]] <- cryo_SCT[["Thawed Village Day 0"]]
cryo_SCT_updated[["Village"]][["Sydney_Fresh"]] <- cryo_SCT[["Village Day 4"]]
cryo_SCT_updated[["Village"]][["Sydney_Cryopreserved"]] <- cryo_SCT[["Thawed Village Day 7"]]


SCT_combined <- location_SCT
SCT_combined[["Baseline"]][["Sydney_Cryopreserved"]] <- cryo_SCT_updated[["Baseline"]][["Sydney_Cryopreserved"]]
SCT_combined[["Village"]][["Sydney_Cryopreserved"]] <- cryo_SCT_updated[["Village"]][["Sydney_Cryopreserved"]]


SCT_combined <- lapply(SCT_combined, function(x){
	lapply(x, function(y){
		subset(y, features = rownames(y)[which(rowSums(y[["SCT"]]@counts > 0)/ncol(y[["SCT"]]@counts) >= 0.01)])
	})
})