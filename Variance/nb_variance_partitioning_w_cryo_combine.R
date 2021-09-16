library(tidyverse)
library(ggplot2)
library(ggrepel)
library(data.table)
library(Seurat)
library(gamlss)



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
	fits_df <- read_delim(paste0(outdir,"nb_variance_partitionin_df.tsv"), delim = "\t")
}


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


pVar_Explained_box <- ggplot(fits_df_long[which(fits_df_long$Location != "Site 3_Cryopreserved"),], aes(x = Covariate, y = Variance_Explained*100, color = Covariate)) +
	geom_boxplot(outlier.size = 0.5) +
	facet_grid(vars(Location)) +
	scale_color_manual(values = variable_colors) +
	theme_classic() +
	ylab("Percent Gene Expression\nVariance Explained") +
	geom_text_repel(size = 2, data = subset(fits_df_long[which(fits_df_long$Location != "Site 3_Cryopreserved"),], 
							((Gene_ID %in% c("CHCHD2", "EIF1AY") & Covariate == "Line") | 
							(Gene_ID %in% c("MT-ATP6", "MT-ND1") & Covariate == "Village"))), 
							aes(label = Gene_ID), box.padding =1, direction = "x")

save_figs(pVar_Explained_box, paste0(outdir, "Total_Contribution_Boxplot_cov"), width = 10, height = 8)

	

fits_df_long[which(fits_df_long$Variance_Explained > 0.5 & fits_df_long$Covariate == "Village"),]
table(fits_df_long[which(fits_df_long$Variance_Explained > 0.4 & fits_df_long$Covariate == "Village"),]$Gene)
table(fits_df_long[which(fits_df_long$Variance_Explained > 0.5 & fits_df_long$Covariate == "Line"),]$Gene)
head(fits_df_long[which(fits_df_long$Covariate == "Village"),][rev(order(na.omit(fits_df_long[which(fits_df_long$Covariate == "Village"),"Variance_Explained"]))),])


### Get average of each covariate contribution to gene expression variation
mean(fits_df_long[which(fits_df_long$Covariate == "Line" & fits_df_long$Location != "Site 3_Cryopreserved"),]$Variance_Explained)
mean(fits_df_long[which(fits_df_long$Covariate == "Village" & fits_df_long$Location != "Site 3_Cryopreserved"),]$Variance_Explained)
mean(fits_df_long[which(fits_df_long$Covariate == "Replicate" & fits_df_long$Location != "Site 3_Cryopreserved"),]$Variance_Explained)



### Variance explained correlation between sites ###
fits_df_joined <- full_join(fits_df_long, fits_df_long, by = c("Gene", "Covariate", "Gene_ID"))
fits_df_joined <- fits_df_joined[which(fits_df_joined$Location.x != fits_df_joined$Location.y),]


### identify the combinatinons to use
pairs_df <- as.data.frame(t(combn(unique(fits_df_joined$Location.x), 2, simplify = TRUE)))
colnames(pairs_df) <- c("Location.x", "Location.y")

fits_df_joined <- left_join(pairs_df, fits_df_joined)
fits_df_joined <- na.omit(fits_df_joined)


Rsquared <- unique(fits_df_joined[,c("Location.x", "Location.y", "Covariate")])

combinations <- data.frame(t(combn(unique(c(fits_df_joined[,c("Location.x")], fits_df_joined[,c("Location.y")])), 2, simplify = TRUE)))
colnames(combinations) <- c("Location.x", "Location.y")

Rsquared <- left_join(combinations, Rsquared)

Rsquared$Rsquared <- NA
Rsquared$`P value` <- NA
Rsquared$pearson <- NA
Rsquared$spearman <- NA
Rsquared$kendall <- NA

for (row in 1:nrow(Rsquared)){
	Rsquared$Rsquared[row] <- round(summary(lm(Variance_Explained.y ~ 0 + Variance_Explained.x, fits_df_joined[which(fits_df_joined$Location.x == Rsquared$Location.x[row] & 
																						fits_df_joined$Location.y == Rsquared$Location.y[row] & 
																						fits_df_joined$Covariate == Rsquared$Covariate[row]),]))$r.squared,2)
	# Rsquared$`P value`[row] <- signif(summary(lm(Variance_Explained.y ~ Variance_Explained.x, fits_df_joined[which(fits_df_joined$Location.x == Rsquared$Location.x[row] & 
	# 																					fits_df_joined$Location.y == Rsquared$Location.y[row] & 
	# 																					fits_df_joined$Covariate == Rsquared$Covariate[row]),]))$coefficients[2,4], digits = 3)
	Rsquared$`P value`[row] <- signif(summary(lm(Variance_Explained.y ~ 0 + Variance_Explained.x, fits_df_joined[which(fits_df_joined$Location.x == Rsquared$Location.x[row] & 
																						fits_df_joined$Location.y == Rsquared$Location.y[row] & 
																						fits_df_joined$Covariate == Rsquared$Covariate[row]),]))$coefficients[1,4], digits = 3)
	Rsquared$pearson[row] <- round(cor(fits_df_joined[which(fits_df_joined$Location.x == Rsquared$Location.x[row] & 
																						fits_df_joined$Location.y == Rsquared$Location.y[row] & 
																						fits_df_joined$Covariate == Rsquared$Covariate[row]),]$Variance_Explained.y,
																						fits_df_joined[which(fits_df_joined$Location.x == Rsquared$Location.x[row] & 
																						fits_df_joined$Location.y == Rsquared$Location.y[row] & 
																						fits_df_joined$Covariate == Rsquared$Covariate[row]),]$Variance_Explained.x, use="complete.obs"), digits = 2)
	Rsquared$spearman[row] <- round(cor(fits_df_joined[which(fits_df_joined$Location.x == Rsquared$Location.x[row] & 
																						fits_df_joined$Location.y == Rsquared$Location.y[row] & 
																						fits_df_joined$Covariate == Rsquared$Covariate[row]),]$Variance_Explained.y,
																						fits_df_joined[which(fits_df_joined$Location.x == Rsquared$Location.x[row] & 
																						fits_df_joined$Location.y == Rsquared$Location.y[row] & 
																						fits_df_joined$Covariate == Rsquared$Covariate[row]),]$Variance_Explained.x, use="complete.obs", method = "spearman"), digits = 2)
	Rsquared$kendall[row] <- round(cor(fits_df_joined[which(fits_df_joined$Location.x == Rsquared$Location.x[row] & 
																						fits_df_joined$Location.y == Rsquared$Location.y[row] & 
																						fits_df_joined$Covariate == Rsquared$Covariate[row]),]$Variance_Explained.y,
																						fits_df_joined[which(fits_df_joined$Location.x == Rsquared$Location.x[row] & 
																						fits_df_joined$Location.y == Rsquared$Location.y[row] & 
																						fits_df_joined$Covariate == Rsquared$Covariate[row]),]$Variance_Explained.x, use="complete.obs", method = "kendall"), digits = 2)
}

Rsquared <- na.omit(Rsquared)
Rsquared <- data.table(Rsquared)

Rsquared_fresh <- Rsquared[Location.y != "Site 3_Cryopreserved"]
Rsquared_cryo <- Rsquared[(Location.y == "Site 3" | Location.y == "Site 3_Cryopreserved") & Location.x == "Site 3"]

Rsquared_fresh[Covariate == "Village"]
Rsquared_fresh[Covariate == "Line"]

t.test(Rsquared_fresh[Covariate == "Line"]$Rsquared, Rsquared_fresh[Covariate == "Village"]$Rsquared, paired = TRUE)

fits_df_joined <- data.table(fits_df_joined)


pCorr_Location_Point_Line <- ggplot(fits_df_joined[Covariate == "Line" & Location.x != "Site 3_Cryopreserved" & Location.y != "Site 3_Cryopreserved"], aes(Variance_Explained.x*100, Variance_Explained.y*100, color = Covariate)) +
	geom_point(size = 0.5, alpha = 0.5) +
	facet_grid(rows = vars(Location.y), cols = vars(Location.x)) +
	scale_color_manual(values = variable_colors[["Line"]]) +
	theme_classic() +
	ylim(0,110)+
	geom_text(size = 2.5, x = 0, y = 105, aes(label = paste0("R = ", pearson)), data = Rsquared[Covariate == "Line" & Location.x != "Site 3_Cryopreserved" & Location.y != "Site 3_Cryopreserved",], hjust = 0, color = "black") +
	xlab("Percent Gene Expression\nVariance Explained") +
	ylab("Percent Gene Expression\nVariance Explained") +
	geom_text_repel(size = 2.5, data = subset(fits_df_joined[Location.x != "Site 3_Cryopreserved" & Location.y != "Site 3_Cryopreserved",], 
							((Gene_ID %in% c("CHCHD2", "EIF1AY") & Covariate == "Line"))), 
							aes(label = Gene_ID), box.padding =0.5, color = "black")
save_figs(pCorr_Location_Point_Line, paste0(outdir, "Correlation_Point_Line"), width = 10, height = 7.25)


pCorr_Location_Point_Village <- ggplot(fits_df_joined[Covariate == "Village" & Location.x != "Site 3_Cryopreserved" & Location.y != "Site 3_Cryopreserved"], aes(Variance_Explained.x*100, Variance_Explained.y*100, color = Covariate)) +
	geom_point(size = 0.5, alpha = 0.5) +
	facet_grid(rows = vars(Location.y), cols = vars(Location.x)) +
	scale_color_manual(values = variable_colors[["Village"]]) +
	theme_classic() +
	ylim(0,75)+
	geom_text(size = 2.5, x = 0, y = 70, aes(label = paste0("R = ", pearson)), data = Rsquared[Covariate == "Village" & Location.x != "Site 3_Cryopreserved" & Location.y != "Site 3_Cryopreserved"], hjust = 0, color = "black") +
	xlab("Percent Gene Expression\nVariance Explained") +
	ylab("Percent Gene Expression\nVariance Explained") 
	# geom_text_repel(size = 2.5, data = subset(fits_df_joined[Location.x != "Site 3_Cryopreserved" & Location.y != "Site 3_Cryopreserved"], 
	# 						(Gene_ID %in% c("MT-ATP6", "MT-ND1") & Covariate == "Village")), 
	# 						aes(label = Gene_ID), box.padding =1, color = "black")
save_figs(pCorr_Location_Point_Village, paste0(outdir, "Correlation_Point_Village"), width = 10, height = 7.25)



##### Test for heteroskadicity (if variance is consistent across x or y) ####
fits_df_joined <- data.table(fits_df_joined)

locations_df <- unique(fits_df_joined[(Covariate == "Village" & Location.x != "Site 3_Cryopreserved" & Location.y != "Site 3_Cryopreserved"),.(Location.x, Location.y)])

heteroskadicity_df_list <- list() 

for (row in 1:nrow(locations_df)){
	heteroskadicity_df_list[["Village"]][[paste0(locations_df[row,Location.x], " - ", locations_df[row,Location.y])]] <- data.table("Residuals" = resid(lm(Variance_Explained.x ~ 0 +Variance_Explained.y, data = fits_df_joined[Covariate == "Village" & Location.x == locations_df[row,Location.x] & Location.y == locations_df[row,Location.y],])), 
								y = fits_df_joined[Covariate == "Village" & Location.x == locations_df[row,Location.x] & Location.y == locations_df[row,Location.y], "Variance_Explained.y"])
	heteroskadicity_df_list[["Line"]][[paste0(locations_df[row,Location.x], locations_df[row,Location.y])]] <- data.table("Residuals" = resid(lm(Variance_Explained.x ~ 0 + Variance_Explained.y, data = fits_df_joined[Covariate == "Line" & Location.x == locations_df[row,Location.x] & Location.y == locations_df[row,Location.y],])), 
								y = fits_df_joined[Covariate == "Line" & Location.x == locations_df[row,Location.x] & Location.y == locations_df[row,Location.y], "Variance_Explained.y"])
}


lapply(names(heteroskadicity_df_list), function(x){
	lapply(names(heteroskadicity_df_list[[x]]), function(y){
		plot <- ggplot(heteroskadicity_df_list[[x]][[y]], aes(y.Variance_Explained.y*100, Residuals)) +
			geom_point() +
			theme_classic()
		ggsave(plot, filename = gsub(" ", "_", paste0(outdir,x, "_",y,"_residuals.png")))
	})
})



##### Use fitDist from gamlss package to compare different distributions for the data and select the best - test first #####
realAll_fit_list <- list()

for (row in 1:nrow(locations_df)){
	realAll_fit_list[["Village"]][[paste0(locations_df[row,Location.x], " - ", locations_df[row,Location.y])]] <- fitDist(fits_df_joined[Covariate == "Village" & Location.x == locations_df[row,Location.x] & Location.y == locations_df[row,Location.y],Variance_Explained.y], k = 2, type = "realAll")
	realAll_fit_list[["Line"]][[paste0(locations_df[row,Location.x], " - ", locations_df[row,Location.y])]] <- fitDist(fits_df_joined[Covariate == "Line" & Location.x == locations_df[row,Location.x] & Location.y == locations_df[row,Location.y],Variance_Explained.y], k = 2, type = "realAll")
}





RbGeneList <- read.delim(file = "/directflow/SCCGGroupShare/projects/DrewNeavin/References/RibosomalGeneList_GeneID_ENSG.txt")
MtGeneList <- read.delim(file = "/directflow/SCCGGroupShare/projects/DrewNeavin/References/MtGeneList_GeneID_ENSG.txt")
fits_df_joined_noMT <- fits_df_joined[-c(which(fits_df_joined$Gene_ID %in% MtGeneList$GeneID)),]
# fits_df_joined_noMT <- fits_df_joined_noMT[-c(which(fits_df_joined_noMT$Gene_ID %in% RbGeneList$GeneID)),]

Rsquared_noMT <- unique(fits_df_joined_noMT[,c("Location.x", "Location.y", "Covariate")])

Rsquared_noMT <- left_join(combinations, Rsquared_noMT)
Rsquared_noMT <- data.table(Rsquared_noMT)

Rsquared_noMT <- Rsquared_noMT[Location.x != Location.y]
Rsquared_noMT <- na.omit(Rsquared_noMT)

Rsquared_noMT$Rsquared <- NA
Rsquared_noMT$`P value` <- NA



for (row in 1:nrow(Rsquared_noMT)){
	Rsquared_noMT$Rsquared[row] <- round(summary(lm(Variance_Explained.y ~ Variance_Explained.x, fits_df_joined_noMT[which(fits_df_joined_noMT$Location.x == Rsquared_noMT$Location.x[row] & 
																						fits_df_joined_noMT$Location.y == Rsquared_noMT$Location.y[row] & 
																						fits_df_joined_noMT$Covariate == Rsquared_noMT$Covariate[row]),]))$r.squared,2)
	Rsquared_noMT$`P value`[row] <- signif(summary(lm(Variance_Explained.y ~ Variance_Explained.x, fits_df_joined_noMT[which(fits_df_joined_noMT$Location.x == Rsquared_noMT$Location.x[row] & 
																						fits_df_joined_noMT$Location.y == Rsquared_noMT$Location.y[row] & 
																						fits_df_joined_noMT$Covariate == Rsquared_noMT$Covariate[row]),]))$coefficients[2,4], digits = 3)
	Rsquared_noMT$pearson[row] <- round(cor(fits_df_joined_noMT[which(fits_df_joined_noMT$Location.x == Rsquared_noMT$Location.x[row] & 
																						fits_df_joined_noMT$Location.y == Rsquared_noMT$Location.y[row] & 
																						fits_df_joined_noMT$Covariate == Rsquared_noMT$Covariate[row]),]$Variance_Explained.y,
																						fits_df_joined_noMT[which(fits_df_joined_noMT$Location.x == Rsquared_noMT$Location.x[row] & 
																						fits_df_joined_noMT$Location.y == Rsquared_noMT$Location.y[row] & 
																						fits_df_joined_noMT$Covariate == Rsquared_noMT$Covariate[row]),]$Variance_Explained.x, use="complete.obs"), digits = 2)

}

Rsquared_noMT <- data.table(Rsquared_noMT)
fits_df_joined_noMT <- data.table(fits_df_joined_noMT)


pCorr_Location_Point_Village_noMT <- ggplot(fits_df_joined_noMT[which(fits_df_joined_noMT$Covariate == "Village" & fits_df_joined_noMT$Location.x != "Site 3_Cryopreserved" & fits_df_joined_noMT$Location.y != "Site 3_Cryopreserved"),], aes(Variance_Explained.x*100, Variance_Explained.y*100, color = Covariate)) +
	geom_point(size = 0.5, alpha = 0.5) +
	facet_grid(rows = vars(Location.y), cols = vars(Location.x), scales = "free") +
	scale_color_manual(values = variable_colors[["Village"]]) +
	theme_classic() +
	ylim(0,75)+
	geom_text(size = 2.5, x = 0, y = 70, aes(label = paste0("R = ", pearson)), data = Rsquared_noMT[which(Rsquared_noMT$Covariate == "Village" & Rsquared_noMT$Location.x != "Site 3_Cryopreserved" & Rsquared_noMT$Location.y != "Site 3_Cryopreserved"),], hjust = 0, color = "black") +
	xlab("Percent Gene Expression\nVariance Explained") +
	ylab("Percent Gene Expression\nVariance Explained") 
save_figs(pCorr_Location_Point_Village_noMT, paste0(outdir, "Correlation_Point_Village_noMT"), width = 10, height = 7.25)


pCorr_Location_Point_Line_noMT <- ggplot(fits_df_joined_noMT[which(fits_df_joined_noMT$Covariate == "Line" & fits_df_joined_noMT$Location.x != "Site 3_Cryopreserved" & fits_df_joined_noMT$Location.y != "Site 3_Cryopreserved"),], aes(Variance_Explained.x*100, Variance_Explained.y*100, color = Covariate)) +
	geom_point(size = 0.5, alpha = 0.5) +
	facet_grid(rows = vars(Location.y), cols = vars(Location.x)) +
	scale_color_manual(values = variable_colors[["Line"]]) +
	theme_classic() +
	ylim(0,110)+
	geom_text(size = 2.5, x = 0, y = 110, aes(label = paste0("R = ", pearson)), data = Rsquared_noMT[which(Rsquared_noMT$Covariate == "Line" & Rsquared_noMT$Location.x != "Site 3_Cryopreserved" & Rsquared_noMT$Location.y != "Site 3_Cryopreserved"),], hjust = 0, color = "black") +
	xlab("Percent Gene Expression\nVariance Explained") +
	ylab("Percent Gene Expression\nVariance Explained") +
	geom_text_repel(size = 2.5, data = subset(fits_df_joined_noMT[which(fits_df_joined_noMT$Location.x != "Site 3_Cryopreserved" & fits_df_joined_noMT$Location.y != "Site 3_Cryopreserved"),], 
							((Gene_ID %in% c("CHCHD2", "EIF1AY") & Covariate == "Line"))), 
							aes(label = Gene_ID), box.padding =0.5, color = "black")
save_figs(pCorr_Location_Point_Line_noMT, paste0(outdir, "Correlation_Point_Line_noMT"), width = 10, height = 7.25)




### CRYOPRESERVED ###
fits_df_long_cryo <- fits_df_long[grepl("Site 3", fits_df_long$Location),]
fits_df_long_cryo$Location <- gsub("Site 3_Cryopreserved", "Cryopreserved", fits_df_long_cryo$Location) %>% gsub("Site 3", "Fresh", .)

fits_df_long_cryo$Location <- factor(fits_df_long_cryo$Location, levels = c("Fresh", "Cryopreserved"))

pTotal_Cont_Cryo <- ggplot(fits_df_long_cryo, aes(Variance_Explained*100, color = Covariate)) +
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
	geom_text_repel(size = 2.5, data = subset(fits_df_long_cryo, 
							((Gene_ID %in% c("CHCHD2", "EIF1AY") & Covariate == "Line") | 
							(Gene_ID %in% c("MT-CO3") & Covariate == "Replicate") | 
							(Gene_ID %in% c("MT-CO3", "MT-ATP6") & Covariate == "Village"))), 
							aes(label = Gene_ID), box.padding =0.7, direction = "x")

save_figs(pVar_Explained_box_Cryo, paste0(outdir, "Total_Contribution_Boxplot_cov_cryo"), width = 10, height = 9)

	
table(subset(fits_df_long_cryo, (Variance_Explained > 0.90))$Gene_ID)
subset(fits_df_long_cryo, (Variance_Explained > 0.90))
fits_df_long_cryo[which(fits_df_long_cryo$Gene_ID == "ZNF71"),]
subset(fits_df_long_cryo[which(fits_df_long_cryo$Covariate == "Replicate"),], Variance_Explained > 0.2)
table(subset(fits_df_long_cryo[which(fits_df_long_cryo$Covariate == "Village"),], Variance_Explained > 0.5)$Gene_ID)



### Variance explained correlation between sites ###
fits_df_joined_cryo <- full_join(fits_df_long_cryo, fits_df_long_cryo, by = c("Gene", "Covariate", "Gene_ID"))
fits_df_joined_cryo <- fits_df_joined_cryo[which(fits_df_joined_cryo$Location.x != fits_df_joined_cryo$Location.y),]


### identify the combinatinons to use
pairs_df <- as.data.table(t(combn(unique(fits_df_joined_cryo$Location.x), 2, simplify = TRUE)))
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
	geom_point(size = 0.5, alpha = 0.5) +
	scale_color_manual(values = variable_colors[["Line"]]) +
	theme_classic() +
	ylim(0,110)+
	geom_text(size = 3, x = 0, y = 110, aes(label = paste0("R = ", pearson)), data = Rsquared_cryo[which(Rsquared_cryo$Covariate == "Line"),], hjust = 0, color = "black") +
	xlab("Percent Gene Expression\nVariance Explained - Fresh") +
	ylab("Percent Gene Expression\nVariance Explained - Cryopreserved") +
	geom_text_repel(size = 3, data = subset(fits_df_joined_cryo[which(fits_df_joined_cryo$Covariate == "Line"),], 
							(Gene_ID %in% c("CHCHD2", "EIF1AY") & Covariate == "Line")), 
							aes(label = Gene_ID), box.padding = 1, color = "black")

save_figs(pCorr_Location_Point_Line_cryo, paste0(outdir, "Correlation_Point_Line_Cryo"), width = 10, height = 7.75)


pCorr_Location_Point_Village_cryo <- ggplot(fits_df_joined_cryo[which(fits_df_joined_cryo$Covariate == "Village"),], aes(Variance_Explained.x*100, Variance_Explained.y*100, color = Covariate)) +
	geom_point(size = 0.5, alpha = 0.5) +
	scale_color_manual(values = variable_colors[["Village"]]) +
	theme_classic() +
	ylim(0,70)+
	geom_text(size = 3, x = 0, y = 70, aes(label = paste0("R = ", pearson)), data = Rsquared_cryo[which(Rsquared_cryo$Covariate == "Village"),], hjust = 0, color = "black") +
	xlab("Percent Gene Expression\nVariance Explained - Fresh") +
	ylab("Percent Gene Expression\nVariance Explained - Cryopreserved") 
	# geom_text_repel(size = 3, data = subset(fits_df_joined_cryo[which(fits_df_joined_cryo$Covariate == "Village"),], 
	# 						(Gene_ID %in% c("MT-CO3", "MT-ATP6") & Covariate == "Village")), 
	# 						aes(label = Gene_ID), box.padding = 1, color = "black")

save_figs(pCorr_Location_Point_Village_cryo, paste0(outdir, "Correlation_Point_Village_Cryo"), width = 10, height = 7.75)



Rsquared_cryo_noMT <- Rsquared_noMT[grepl("Site 3", Rsquared_noMT$Location.x) & grepl("Site 3", Rsquared_noMT$Location.y)]

Rsquared_cryo_noMT$Location.x <- gsub("Site 3", "Fresh", Rsquared_cryo_noMT$Location.x)
Rsquared_cryo_noMT$Location.y <- gsub("Site 3_Cryopreserved", "Cryopreserved", Rsquared_cryo_noMT$Location.y)


# RbGeneList <- read.delim(file = "/directflow/SCCGGroupShare/projects/DrewNeavin/References/RibosomalGeneList_GeneID_ENSG.txt")
MtGeneList <- read.delim(file = "/directflow/SCCGGroupShare/projects/DrewNeavin/References/MtGeneList_GeneID_ENSG.txt")
fits_df_joined_cryo_noMT <- fits_df_joined_cryo[-c(which(fits_df_joined_cryo$Gene_ID %in% MtGeneList$GeneID)),]
# fits_df_joined_cryo_noMT <- fits_df_joined_cryo_noMT[-c(which(fits_df_joined_cryo_noMT$Gene_ID %in% RbGeneList$GeneID)),]



pCorr_Location_Point_Village_cryo_noMT <- ggplot(fits_df_joined_cryo_noMT[which(fits_df_joined_cryo_noMT$Covariate == "Village"),], aes(Variance_Explained.x*100, Variance_Explained.y*100, color = Covariate)) +
	geom_point(size = 0.5, alpha = 0.5) +
	scale_color_manual(values = variable_colors[["Village"]]) +
	theme_classic() +
	ylim(0,70)+
	geom_text(size = 3, x = 0, y = 70, aes(label = paste0("R = ", pearson)), data = Rsquared_cryo_noMT[which(Rsquared_cryo_noMT$Covariate == "Village"),], hjust = 0, color = "black") +
	xlab("Percent Gene Expression\nVariance Explained - Fresh") +
	ylab("Percent Gene Expression\nVariance Explained - Cryopreserved")

save_figs(pCorr_Location_Point_Village_cryo_noMT, paste0(outdir, "Correlation_Point_Village_Cryo_noMT"), width = 10, height = 7.75)


pCorr_Location_Point_Line_cryo_noMT <- ggplot(fits_df_joined_cryo_noMT[which(fits_df_joined_cryo_noMT$Covariate == "Line"),], aes(Variance_Explained.x*100, Variance_Explained.y*100, color = Covariate)) +
	geom_point(size = 0.5, alpha = 0.5) +
	scale_color_manual(values = variable_colors[["Line"]]) +
	theme_classic() +
	ylim(0,110)+
	geom_text(size = 3, x = 0, y = 110, aes(label = paste0("R = ", pearson)), data = Rsquared_cryo_noMT[which(Rsquared_cryo_noMT$Covariate == "Line"),], hjust = 0, color = "black") +
	xlab("Percent Gene Expression\nVariance Explained - Fresh") +
	ylab("Percent Gene Expression\nVariance Explained - Cryopreserved") +
	geom_text_repel(size = 3, data = subset(fits_df_joined_cryo_noMT[which(fits_df_joined_cryo_noMT$Covariate == "Line"),], 
							(Gene_ID %in% c("CHCHD2", "EIF1AY") & Covariate == "Line")), 
							aes(label = Gene_ID), box.padding = 1, color = "black")

save_figs(pCorr_Location_Point_Line_cryo_noMT, paste0(outdir, "Correlation_Point_Line_Cryo_noMT"), width = 10, height = 7.75)
