library(tidyverse)
library(Seurat)
library(matrixStats)
library(ggpubr)
library(data.table)


##### Set Up Directories #####
dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
datadir <- paste0(dir,"output/All_data_integrated_remove_bad/")
outdir <- paste0(dir,"output/Expression_Correlations/")

dir.create(outdir)



##### Set Colors for Figures #####
line_colors <- c(FSA0006 = "#F79E29", MBE1006 = "#9B2C99", TOB0421 = "#35369C")

site_updates <- c("Brisbane" = "Site 1", "Sydney" = "Site 3" ,"Melbourne" = "Site 2")


##### Figure Function #####
save_figs <- function(plot, basename, width = 17, height = 17, units = "cm"){
    ggsave(plot, filename = paste0(basename,".png"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".pdf"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".eps"), height = height, width = width, units = units)
}




##### Read in Data #####
seurat <- readRDS(paste0(datadir,"seurat_integrated_all_times_clustered.rds"))

colnames(seurat@meta.data)
unique(seurat@meta.data$Location)
unique(seurat@meta.data$Time)
unique(seurat@meta.data$Final_Assignment)
unique(seurat@meta.data$MULTI_ID)

seurat@meta.data$Location <- gsub("_.+", "", seurat@meta.data$Location)


seurat_list <- lapply(unique(seurat@meta.data$Location), function(x){
	subset(seurat, subset = Location == x)
})
names(seurat_list) <- unique(seurat@meta.data$Location)


seurat_list_list <- lapply(names(seurat_list), function(x){
	temp <- list()
	for (day in unique(seurat_list[[x]]@meta.data$Time)){
		temp[[x]][[day]] <- subset(seurat_list[[x]], subset = Time == day)
	}
	return(temp[[x]])
})
names(seurat_list_list) <- unique(seurat@meta.data$Location)


seurat_list_list_list <- lapply(names(seurat_list_list), function(x){
	temp2 <- lapply(names(seurat_list_list[[x]]), function(y){
		temp <- list()
		for (assign in unique(seurat_list_list[[x]][[y]]$Final_Assignment)){
			temp[[x]][[y]][[assign]] <- subset(seurat_list_list[[x]][[y]], subset = Final_Assignment == assign)
		}
		return(temp[[x]][[y]])
	})
	names(temp2) <- unique(seurat_list[[x]]@meta.data$Time)
	return(temp2)
})
names(seurat_list_list_list) <- unique(seurat@meta.data$Location)

seurat_list_list_list_list <- lapply(names(seurat_list_list), function(x){
	temp2 <- lapply(names(seurat_list_list_list[[x]]), function(y){
		temp3 <- lapply(names(seurat_list_list_list[[x]][[y]]), function(z){
			temp <- list()
			for (rep in unique(seurat_list_list_list[[x]][[y]][[z]]$MULTI_ID)){
				temp[[x]][[y]][[z]][[rep]] <- subset(seurat_list_list_list[[x]][[y]][[z]], subset = MULTI_ID == rep)
			}
			return(temp[[x]][[y]][[z]])
		})
		print(temp3)
		names(temp3) <- unique(seurat_list_list[[x]][[y]]$Final_Assignment)
		return(temp3)
	})
	names(temp2) <- unique(seurat_list[[x]]@meta.data$Time)
	return(temp2)
})
names(seurat_list_list_list_list) <- unique(seurat@meta.data$Location)


### 10% to start
seurat_list_list_list_list_sub <- lapply(seurat_list_list_list_list, function(x){
	lapply(x, function(y){
		lapply(y, function(z){
			lapply(z, function(rep){
				subset(rep, features = rownames(seurat)[which((rowSums(seurat[["SCT"]]@counts > 0)/ncol(seurat[["SCT"]]@counts)) > 0.1)])
			})
		})
	})
})
names(seurat_list_list_list_list_sub) <- unique(seurat_sub@meta.data$Location)

saveRDS(seurat_list_list_list_list_sub, paste0(outdir, "subset_seurat_list.rds"))
seurat_list_list_list_list_sub <- readRDS(paste0(outdir, "subset_seurat_list.rds"))


summary_list <- lapply(names(seurat_list_list_list_list), function(x){
	temp3 <- lapply(names(seurat_list_list_list_list[[x]]), function(y){
		temp2 <- lapply(names(seurat_list_list_list_list[[x]][[y]]), function(z){
			temp <- lapply(names(seurat_list_list_list_list[[x]][[y]][[z]]), function(rep){
				data.frame(Gene = rownames(seurat_list_list_list_list[[x]][[y]][[z]][[rep]]), Mean = rowMeans(seurat_list_list_list_list[[x]][[y]][[z]][[rep]][["SCT"]]@scale.data), Replicate = rep, Line = z, Time = y, Location = x)
			})
			# print(head(temp))
			do.call(rbind, temp)
		})
		do.call(rbind, temp2)
	})
	do.call(rbind, temp3)
})

summary <- do.call(rbind, summary_list)
summary$Replicate <- gsub("Brisbane", "Replicate", summary$Replicate) %>% gsub("Melbourne", "Replicate", .) %>% gsub("Sydney", "Replicate", .)
summary_wide <- pivot_wider(summary, names_from = Replicate, values_from = Mean)
summary_wide$SD <- rowSds(as.matrix(summary_wide[,c("Replicate1", "Replicate2", "Replicate3")]))
summary_wide$Mean <- rowMeans(as.matrix(summary_wide[,c("Replicate1", "Replicate2", "Replicate3")]))


summary_wide_wide <- inner_join(summary_wide, summary_wide, by = c("Gene", "Line", "Location"))
summary_wide_wide <- summary_wide_wide[which(summary_wide_wide$Time.x != summary_wide_wide$Time.y),]
summary_wide_wide <- summary_wide_wide[which(summary_wide_wide$Time.x == "Baseline"),]



pExpression_Correlation <- list()
for (location in unique(summary_wide_wide$Location)){
	pExpression_Correlation[[location]] <- ggplot(summary_wide_wide[which(summary_wide_wide$Location == location),], aes(Mean.x, Mean.y)) +
		geom_point() +
		theme_classic()
	
	ggsave(pExpression_Correlation[[location]], filename = paste0(outdir,location,"_Expression_Correlation.png"))
}







##### Redo with separately assessed by SCT #####
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



SCT_combined_list <- lapply(names(SCT_combined), function(x){
	temp2 <- lapply(names(SCT_combined[[x]]), function(y){
		temp <- list()
		for (assign in unique(SCT_combined[[x]][[y]]$Final_Assignment)){
			temp[[x]][[y]][[assign]] <- subset(SCT_combined[[x]][[y]], subset = Final_Assignment == assign)
		}
		return(temp[[x]][[y]])
	})
	names(temp2) <- unique(names(SCT_combined[[x]]))
	return(temp2)
})
names(SCT_combined_list) <- names(SCT_combined)


SCT_combined_list_list <- lapply(names(SCT_combined_list), function(x){
	temp2 <- lapply(names(SCT_combined_list[[x]]), function(y){
		temp3 <- lapply(names(SCT_combined_list[[x]][[y]]), function(z){
			temp <- list()
			for (rep in unique(SCT_combined_list[[x]][[y]][[z]]$MULTI_ID)){
				temp[[x]][[y]][[z]][[rep]] <- subset(SCT_combined_list[[x]][[y]][[z]], subset = MULTI_ID == rep)
			}
			return(temp[[x]][[y]][[z]])
		})
		print(temp3)
		names(temp3) <- unique(SCT_combined[[x]][[y]]$Final_Assignment)
		return(temp3)
	})
	names(temp2) <- names(SCT_combined[[x]])
	return(temp2)
})
names(SCT_combined_list_list) <- names(SCT_combined)


saveRDS(SCT_combined_list_list, paste0(outdir, "subset_seurat_list_sep.rds"))
SCT_combined_list_list <- readRDS(paste0(outdir, "subset_seurat_list_sep.rds"))


summary_SCT_list <- lapply(names(SCT_combined_list_list), function(x){
	temp3 <- lapply(names(SCT_combined_list_list[[x]]), function(y){
		temp2 <- lapply(names(SCT_combined_list_list[[x]][[y]]), function(z){
			temp <- lapply(names(SCT_combined_list_list[[x]][[y]][[z]]), function(rep){
				data.frame(Gene = rownames(SCT_combined_list_list[[x]][[y]][[z]][[rep]]), Mean = rowMeans(SCT_combined_list_list[[x]][[y]][[z]][[rep]][["SCT"]]@scale.data), N = ncol(SCT_combined_list_list[[x]][[y]][[z]][[rep]][["SCT"]]@scale.data), Replicate = rep, Line = z, Location = y, Time = x)
			})
			# print(head(temp))
			do.call(rbind, temp)
		})
		do.call(rbind, temp2)
	})
	do.call(rbind, temp3)
})



summary_SCT <- do.call(rbind, summary_SCT_list)
summary_SCT$Replicate <- gsub("Brisbane", "Replicate", summary_SCT$Replicate) %>% gsub("Melbourne", "Replicate", .) %>% gsub("Sydney", "Replicate", .)
summary_SCT_wide <- pivot_wider(summary_SCT, names_from = Replicate, values_from = c(Mean, N))
summary_SCT_wide$SD <- rowSds(as.matrix(summary_SCT_wide[,c("Mean_Replicate1", "Mean_Replicate2", "Mean_Replicate3")]))
summary_SCT_wide$Mean <- (summary_SCT_wide$Mean_Replicate1 * summary_SCT_wide$N_Replicate1 + summary_SCT_wide$Mean_Replicate2 * summary_SCT_wide$N_Replicate2 + summary_SCT_wide$Mean_Replicate3 * summary_SCT_wide$N_Replicate3)/(summary_SCT_wide$N_Replicate1 + summary_SCT_wide$N_Replicate2 + summary_SCT_wide$N_Replicate3)


summary_SCT_wide_wide <- inner_join(summary_SCT_wide, summary_SCT_wide, by = c("Gene", "Line", "Location"))
summary_SCT_wide_wide <- summary_SCT_wide_wide[which(summary_SCT_wide_wide$Time.x != summary_SCT_wide_wide$Time.y),]
summary_SCT_wide_wide <- summary_SCT_wide_wide[which(summary_SCT_wide_wide$Time.x == "Baseline"),]

as.data.frame(summary_SCT_wide_wide[which(summary_SCT_wide_wide$Line == "TOB0421" & summary_SCT_wide_wide$Mean.y < -2.5), ]) ### ENSG00000129824
as.data.frame(summary_SCT_wide_wide[which(summary_SCT_wide_wide$Gene == "ENSG00000129824"),])
as.data.frame(summary_SCT_wide_wide[which(summary_SCT_wide_wide$Line == "TOB0421" & summary_SCT_wide_wide$Mean.y >4), ]) ### ENSG00000106153 CHCHD2
as.data.frame(summary_SCT_wide_wide[which(summary_SCT_wide_wide$Gene == "ENSG00000106153"),])


pExpression_Correlation_SCT <- list()
for (location in unique(summary_SCT_wide_wide$Location)){
	pExpression_Correlation_SCT[[location]] <- ggplot(summary_SCT_wide_wide[which(summary_SCT_wide_wide$Location == location),], aes(Mean.x, Mean.y)) +
		geom_point() +
		theme_classic() +
		facet_wrap(vars(Line), nrow = 1, scales = "free")
	
	ggsave(pExpression_Correlation_SCT[[location]], filename = paste0(outdir,location,"_Expression_Correlation_norm_sep.png"), width = 15)
}


Rsquared <- unique(summary_SCT_wide_wide[,c("Line","Location","Time.x","Time.y")])
Rsquared$Rsquared <- NA
Rsquared$pearson
Rsquared$spearman

for (row in 1:nrow(Rsquared)){
	Rsquared$Rsquared[row] <- summary(lm(Mean.y ~ Mean.x, summary_SCT_wide_wide[which(summary_SCT_wide_wide$Line == Rsquared$Line[row] & 
																						summary_SCT_wide_wide$Location == Rsquared$Location[row] & 
																						summary_SCT_wide_wide$Time.y == Rsquared$Time.y[row] &
																						summary_SCT_wide_wide$Time.x == Rsquared$Time.x[row]),]))$r.squared
	Rsquared$pearson[row] <- cor(summary_SCT_wide_wide[which(summary_SCT_wide_wide$Line == Rsquared$Line[row] & 
																						summary_SCT_wide_wide$Location == Rsquared$Location[row] & 
																						summary_SCT_wide_wide$Time.y == Rsquared$Time.y[row] &
																						summary_SCT_wide_wide$Time.x == Rsquared$Time.x[row]),]$Mean.y, 
							summary_SCT_wide_wide[which(summary_SCT_wide_wide$Line == Rsquared$Line[row] & 
																						summary_SCT_wide_wide$Location == Rsquared$Location[row] & 
																						summary_SCT_wide_wide$Time.y == Rsquared$Time.y[row] &
																						summary_SCT_wide_wide$Time.x == Rsquared$Time.x[row]),]$Mean.x)
	Rsquared$spearman[row] <- cor(summary_SCT_wide_wide[which(summary_SCT_wide_wide$Line == Rsquared$Line[row] & 
																						summary_SCT_wide_wide$Location == Rsquared$Location[row] & 
																						summary_SCT_wide_wide$Time.y == Rsquared$Time.y[row] &
																						summary_SCT_wide_wide$Time.x == Rsquared$Time.x[row]),]$Mean.y, 
							summary_SCT_wide_wide[which(summary_SCT_wide_wide$Line == Rsquared$Line[row] & 
																						summary_SCT_wide_wide$Location == Rsquared$Location[row] & 
																						summary_SCT_wide_wide$Time.y == Rsquared$Time.y[row] &
																						summary_SCT_wide_wide$Time.x == Rsquared$Time.x[row]),]$Mean.x, method = "spearman")
}

Rsquared







# ##### Filter out the genes not expressed in 1% of cells #####
# SCT_combined_1pct <- lapply(SCT_combined, function(x){
# 	lapply(x, function(y){
# 		subset(y, features = rownames(y)[which((rowSums(y[["SCT"]]@counts > 0)/ncol(y[["SCT"]]@counts)) > 0.01)])
# 	})
# })



# summary_SCT_list_1pct <- lapply(names(SCT_combined_1pct), function(village){
# 	temp <- lapply(names(SCT_combined_1pct[[village]]), function(location){
# 		summary_SCT[which(summary_SCT$Location == location & summary_SCT$Time == village & summary_SCT$Gene %in% rownames(SCT_combined_1pct[[village]][[location]])), ]
# 	})
# 	do.call(rbind, temp)
# })


# summary_SCT_1pct <- do.call(rbind, summary_SCT_list_1pct)
# summary_SCT_1pct$Replicate <- gsub("Brisbane", "Replicate", summary_SCT_1pct$Replicate) %>% gsub("Melbourne", "Replicate", .) %>% gsub("Sydney", "Replicate", .)
# summary_SCT_1pct_wide <- pivot_wider(summary_SCT_1pct, names_from = Replicate, values_from = c(Mean, N))
# summary_SCT_1pct_wide$SD <- rowSds(as.matrix(summary_SCT_1pct_wide[,c("Mean_Replicate1", "Mean_Replicate2", "Mean_Replicate3")]))
# summary_SCT_1pct_wide$Mean <- (summary_SCT_1pct_wide$Mean_Replicate1 * summary_SCT_1pct_wide$N_Replicate1 + summary_SCT_1pct_wide$Mean_Replicate2 * summary_SCT_1pct_wide$N_Replicate2 + summary_SCT_1pct_wide$Mean_Replicate3 * summary_SCT_1pct_wide$N_Replicate3)/(summary_SCT_1pct_wide$N_Replicate1 + summary_SCT_1pct_wide$N_Replicate2 + summary_SCT_1pct_wide$N_Replicate3)


# summary_SCT_1pct_wide_wide <- inner_join(summary_SCT_1pct_wide, summary_SCT_1pct_wide, by = c("Gene", "Line", "Location"))
# summary_SCT_1pct_wide_wide <- summary_SCT_1pct_wide_wide[which(summary_SCT_1pct_wide_wide$Time.x != summary_SCT_1pct_wide_wide$Time.y),]
# summary_SCT_1pct_wide_wide <- summary_SCT_1pct_wide_wide[which(summary_SCT_1pct_wide_wide$Time.x == "Baseline"),]

# as.data.frame(summary_SCT_1pct_wide_wide[which(summary_SCT_1pct_wide_wide$Line == "TOB0421" & summary_SCT_1pct_wide_wide$Mean.y < -2.5), ]) ### ENSG00000129824
# as.data.frame(summary_SCT_1pct_wide_wide[which(summary_SCT_1pct_wide_wide$Gene == "ENSG00000129824"),])
# as.data.frame(summary_SCT_1pct_wide_wide[which(summary_SCT_1pct_wide_wide$Line == "TOB0421" & summary_SCT_1pct_wide_wide$Mean.y >4), ]) ### ENSG00000106153 CHCHD2
# as.data.frame(summary_SCT_1pct_wide_wide[which(summary_SCT_1pct_wide_wide$Gene == "ENSG00000106153"),])


# pExpression_Correlation_SCT <- list()
# for (location in unique(summary_SCT_1pct_wide_wide$Location)){
# 	pExpression_Correlation_SCT[[location]] <- ggplot(summary_SCT_1pct_wide_wide[which(summary_SCT_1pct_wide_wide$Location == location),], aes(Mean.x, Mean.y)) +
# 		geom_point() +
# 		theme_classic() +
# 		facet_wrap(vars(Line), nrow = 1, scales = "free")
	
# 	ggsave(pExpression_Correlation_SCT[[location]], filename = paste0(outdir,location,"_Expression_Correlation_norm_sep.png"), width = 15)
# }


# Rsquared_1pct <- unique(summary_SCT_1pct_wide_wide[,c("Line","Location","Time.x","Time.y")])
# Rsquared_1pct$Rsquared <- NA
# Rsquared_1pct$pearson
# Rsquared_1pct$spearman

# for (row in 1:nrow(Rsquared_1pct)){
# 	Rsquared_1pct$Rsquared[row] <- summary(lm(Mean.y ~ Mean.x, summary_SCT_1pct_wide_wide[which(summary_SCT_1pct_wide_wide$Line == Rsquared_1pct$Line[row] & 
# 																						summary_SCT_1pct_wide_wide$Location == Rsquared_1pct$Location[row] & 
# 																						summary_SCT_1pct_wide_wide$Time.y == Rsquared_1pct$Time.y[row] &
# 																						summary_SCT_1pct_wide_wide$Time.x == Rsquared_1pct$Time.x[row]),]))$r.squared
# 	Rsquared_1pct$pearson[row] <- cor(summary_SCT_1pct_wide_wide[which(summary_SCT_1pct_wide_wide$Line == Rsquared_1pct$Line[row] & 
# 																						summary_SCT_1pct_wide_wide$Location == Rsquared_1pct$Location[row] & 
# 																						summary_SCT_1pct_wide_wide$Time.y == Rsquared_1pct$Time.y[row] &
# 																						summary_SCT_1pct_wide_wide$Time.x == Rsquared_1pct$Time.x[row]),]$Mean.y, 
# 							summary_SCT_1pct_wide_wide[which(summary_SCT_1pct_wide_wide$Line == Rsquared_1pct$Line[row] & 
# 																						summary_SCT_1pct_wide_wide$Location == Rsquared_1pct$Location[row] & 
# 																						summary_SCT_1pct_wide_wide$Time.y == Rsquared_1pct$Time.y[row] &
# 																						summary_SCT_1pct_wide_wide$Time.x == Rsquared_1pct$Time.x[row]),]$Mean.x)
# 	Rsquared_1pct$spearman[row] <- cor(summary_SCT_1pct_wide_wide[which(summary_SCT_1pct_wide_wide$Line == Rsquared_1pct$Line[row] & 
# 																						summary_SCT_1pct_wide_wide$Location == Rsquared_1pct$Location[row] & 
# 																						summary_SCT_1pct_wide_wide$Time.y == Rsquared_1pct$Time.y[row] &
# 																						summary_SCT_1pct_wide_wide$Time.x == Rsquared_1pct$Time.x[row]),]$Mean.y, 
# 							summary_SCT_1pct_wide_wide[which(summary_SCT_1pct_wide_wide$Line == Rsquared_1pct$Line[row] & 
# 																						summary_SCT_1pct_wide_wide$Location == Rsquared_1pct$Location[row] & 
# 																						summary_SCT_1pct_wide_wide$Time.y == Rsquared_1pct$Time.y[row] &
# 																						summary_SCT_1pct_wide_wide$Time.x == Rsquared_1pct$Time.x[row]),]$Mean.x, method = "spearman")
# }

# Rsquared_1pct







##### Filter out the genes not expressed in 10% of cells #####
SCT_combined_10pct <- lapply(SCT_combined, function(x){
	lapply(x, function(y){
		subset(y, features = rownames(y)[which((rowSums(y[["SCT"]]@counts > 0)/ncol(y[["SCT"]]@counts)) > 0.1)])
	})
})



summary_SCT_list_10pct <- lapply(names(SCT_combined_10pct), function(village){
	temp <- lapply(names(SCT_combined_10pct[[village]]), function(location){
		summary_SCT[which(summary_SCT$Location == location & summary_SCT$Time == village & summary_SCT$Gene %in% rownames(SCT_combined_10pct[[village]][[location]])), ]
	})
	do.call(rbind, temp)
})


summary_SCT_10pct <- do.call(rbind, summary_SCT_list_10pct)
summary_SCT_10pct$Replicate <- gsub("Brisbane", "Replicate", summary_SCT_10pct$Replicate) %>% gsub("Melbourne", "Replicate", .) %>% gsub("Sydney", "Replicate", .)
summary_SCT_10pct_wide <- pivot_wider(summary_SCT_10pct, names_from = Replicate, values_from = c(Mean, N))
summary_SCT_10pct_wide$SD <- rowSds(as.matrix(summary_SCT_10pct_wide[,c("Mean_Replicate1", "Mean_Replicate2", "Mean_Replicate3")]))
summary_SCT_10pct_wide$Mean <- (summary_SCT_10pct_wide$Mean_Replicate1 * summary_SCT_10pct_wide$N_Replicate1 + summary_SCT_10pct_wide$Mean_Replicate2 * summary_SCT_10pct_wide$N_Replicate2 + summary_SCT_10pct_wide$Mean_Replicate3 * summary_SCT_10pct_wide$N_Replicate3)/(summary_SCT_10pct_wide$N_Replicate1 + summary_SCT_10pct_wide$N_Replicate2 + summary_SCT_10pct_wide$N_Replicate3)


summary_SCT_10pct_wide_wide <- inner_join(summary_SCT_10pct_wide, summary_SCT_10pct_wide, by = c("Gene", "Line", "Location"))
summary_SCT_10pct_wide_wide <- summary_SCT_10pct_wide_wide[which(summary_SCT_10pct_wide_wide$Time.x != summary_SCT_10pct_wide_wide$Time.y),]
summary_SCT_10pct_wide_wide <- summary_SCT_10pct_wide_wide[which(summary_SCT_10pct_wide_wide$Time.x == "Baseline"),]

as.data.frame(summary_SCT_10pct_wide_wide[which(summary_SCT_10pct_wide_wide$Line == "TOB0421" & summary_SCT_10pct_wide_wide$Mean.y < -2.5), ]) ### ENSG00000129824
as.data.frame(summary_SCT_10pct_wide_wide[which(summary_SCT_10pct_wide_wide$Gene == "ENSG00000129824"),])
as.data.frame(summary_SCT_10pct_wide_wide[which(summary_SCT_10pct_wide_wide$Line == "TOB0421" & summary_SCT_10pct_wide_wide$Mean.y >4), ]) ### ENSG00000106153 CHCHD2
as.data.frame(summary_SCT_10pct_wide_wide[which(summary_SCT_10pct_wide_wide$Gene == "ENSG00000106153"),])


pExpression_Correlation_SCT <- list()
for (location in unique(summary_SCT_10pct_wide_wide$Location)){
	pExpression_Correlation_SCT[[location]] <- ggplot(summary_SCT_10pct_wide_wide[which(summary_SCT_10pct_wide_wide$Location == location),], aes(Mean.x, Mean.y)) +
		geom_point(size = 0.5, alpha = 0.5) +
		theme_classic() +
		facet_wrap(vars(Line), nrow = 1, scales = "free")
	
	ggsave(pExpression_Correlation_SCT[[location]], filename = paste0(outdir,location,"_Expression_Correlation_norm_sep_10pct.png"), width = 15)
}


pExpression_Correlation_Combined <- ggplot(summary_SCT_10pct_wide_wide, aes(Mean.x, Mean.y)) +
		geom_point(size = 0.5, alpha = 0.3) +
		theme_classic() +
		facet_grid(Location ~ Line, scales = "free")


pExpression_Correlation_Combined <- ggscatter(data = summary_SCT_10pct_wide_wide, x = "Mean.x", y = "Mean.y",
												color = "Line",
												facet.by = c("Location", "Line"),
												scales = "free",
												size = 0.5,
												alpha = 0.25, 
												ylab = "Average Baseline Normalized Expression",
												xlab = "Average Village Normalized Expression",
												add = "reg.line", 
												conf.int = TRUE,
												scales = "free") +
												stat_cor(aes(color = Line), method = "pearson")
ggsave(pExpression_Correlation_Combined, filename = paste0(outdir,"Expression_Correlation_norm_sep_10pct_all_locations.png"))


summary_SCT_10pct_wide_wide_fresh <- summary_SCT_10pct_wide_wide[which(summary_SCT_10pct_wide_wide$Location != "Sydney_Cryopreserved"),] 
summary_SCT_10pct_wide_wide_fresh$Location <- gsub("_Fresh", "", summary_SCT_10pct_wide_wide_fresh$Location)


for (location in names(site_updates)){
	summary_SCT_10pct_wide_wide_fresh$Location <- gsub(location, site_updates[location], summary_SCT_10pct_wide_wide_fresh$Location)
}

summary_SCT_10pct_wide_wide_fresh$Location_Line <- paste0(summary_SCT_10pct_wide_wide_fresh$Location, "_", summary_SCT_10pct_wide_wide_fresh$Line)

pExpression_Correlation_Combined_fresh <- ggscatter(data = summary_SCT_10pct_wide_wide_fresh, x = "Mean.x", y = "Mean.y",
												color = "Line",
												facet.by = c("Location_Line"),
												palette = line_colors,
												size = 0.5,
												alpha = 0.25, 
												ylab = "Baseline Average Normalized Expression",
												xlab = "Village Average Normalized Expression",
												add = "reg.line", 
												add.params = list(color = "black", fill = "lightgray", size = 0.5),
												conf.int = TRUE,
												scales = "free") +
												stat_cor(aes(label = ..r.label..), method = "pearson")
save_figs(pExpression_Correlation_Combined_fresh, paste0(outdir,"Expression_Correlation_norm_sep_10pct_fresh"), height = 15, width = 15)


### Make accompanying table ###
expression_pearson_df <- unique(summary_SCT_wide_wide[,c("Line","Location","Time.x","Time.y")])
expression_pearson_df$`Pearson R` <- NA
expression_pearson_df$`Pearson P` <- NA

for (row in 1:nrow(expression_pearson_df)){
	expression_pearson_df$`Pearson R`[row] <- cor.test(summary_SCT_wide_wide[which(summary_SCT_wide_wide$Line == expression_pearson_df$Line[row] & 
																						summary_SCT_wide_wide$Location == expression_pearson_df$Location[row] & 
																						summary_SCT_wide_wide$Time.y == expression_pearson_df$Time.y[row] &
																						summary_SCT_wide_wide$Time.x == expression_pearson_df$Time.x[row]),]$Mean.y, 
							summary_SCT_wide_wide[which(summary_SCT_wide_wide$Line == expression_pearson_df$Line[row] & 
																						summary_SCT_wide_wide$Location == expression_pearson_df$Location[row] & 
																						summary_SCT_wide_wide$Time.y == expression_pearson_df$Time.y[row] &
																						summary_SCT_wide_wide$Time.x == expression_pearson_df$Time.x[row]),]$Mean.x, exact = TRUE, method = "pearson")$estimate
																					
	expression_pearson_df$`Pearson P`[row] <- cor.test(summary_SCT_wide_wide[which(summary_SCT_wide_wide$Line == expression_pearson_df$Line[row] & 
																						summary_SCT_wide_wide$Location == expression_pearson_df$Location[row] & 
																						summary_SCT_wide_wide$Time.y == expression_pearson_df$Time.y[row] &
																						summary_SCT_wide_wide$Time.x == expression_pearson_df$Time.x[row]),]$Mean.y, 
							summary_SCT_wide_wide[which(summary_SCT_wide_wide$Line == expression_pearson_df$Line[row] & 
																						summary_SCT_wide_wide$Location == expression_pearson_df$Location[row] & 
																						summary_SCT_wide_wide$Time.y == expression_pearson_df$Time.y[row] &
																						summary_SCT_wide_wide$Time.x == expression_pearson_df$Time.x[row]),]$Mean.x, exact = TRUE, method = "pearson")$p.value
}

for (location in names(site_updates)){
	expression_pearson_df$Location <- gsub(location, site_updates[location], expression_pearson_df$Location)
}
expression_pearson_df <- data.table(expression_pearson_df)

expression_pearson_df_fresh <- expression_pearson_df[Location != "Site 3_Cryopreserved"]
expression_pearson_df_fresh$Location <- gsub("_Fresh", "", expression_pearson_df_fresh$Location)

fwrite(expression_pearson_df_fresh, paste0(outdir, "Fresh_expression_pearson.tsv"), sep = "\t")


expression_pearson_df_cryo <- expression_pearson_df[grepl("Site 3", expression_pearson_df$Location)]
expression_pearson_df_cryo$Location <- gsub("_", " - ", expression_pearson_df_cryo$Location)

fwrite(expression_pearson_df_cryo, paste0(outdir, "Cryopreserved_expression_pearson.tsv"), sep = "\t")




summary_SCT_10pct_wide_wide_cryo <- summary_SCT_10pct_wide_wide[grepl("Sydney", summary_SCT_10pct_wide_wide$Location),] 
summary_SCT_10pct_wide_wide_cryo$Location <- gsub("Sydney_", "", summary_SCT_10pct_wide_wide_cryo$Location)
summary_SCT_10pct_wide_wide_cryo$Location <- factor(summary_SCT_10pct_wide_wide_cryo$Location, levels = c("Fresh", "Cryopreserved"))

for (location in names(site_updates)){
	summary_SCT_10pct_wide_wide_cryo$Location <- gsub(location, site_updates[location], summary_SCT_10pct_wide_wide_cryo$Location)
}

summary_SCT_10pct_wide_wide_cryo$Location_Line <- paste0(summary_SCT_10pct_wide_wide_cryo$Location, "_", summary_SCT_10pct_wide_wide_cryo$Line)


pExpression_Correlation_Combined_cryo <- ggscatter(data = summary_SCT_10pct_wide_wide_cryo, x = "Mean.x", y = "Mean.y",
												color = "Line",
												facet.by = c("Location_Line"),
												palette = line_colors,
												size = 0.5,
												alpha = 0.25, 
												ylab = "Baseline Average Normalized Expression",
												xlab = "Village Average Normalized Expression",
												add = "reg.line", 
												add.params = list(color = "black", fill = "lightgray", size = 0.5),
												conf.int = TRUE,
												scales = "free") +
												stat_cor(aes(label = ..r.label..), method = "pearson")
save_figs(pExpression_Correlation_Combined_cryo, paste0(outdir,"Expression_Correlation_norm_sep_10pct_cryo"), height = 10, width = 12)



Rsquared_10pct <- unique(summary_SCT_10pct_wide_wide[,c("Line","Location","Time.x","Time.y")])
Rsquared_10pct$Rsquared <- NA
Rsquared_10pct$pearson
Rsquared_10pct$spearman

for (row in 1:nrow(Rsquared_10pct)){
	Rsquared_10pct$Rsquared[row] <- summary(lm(Mean.y ~ Mean.x, summary_SCT_10pct_wide_wide[which(summary_SCT_10pct_wide_wide$Line == Rsquared_10pct$Line[row] & 
																						summary_SCT_10pct_wide_wide$Location == Rsquared_10pct$Location[row] & 
																						summary_SCT_10pct_wide_wide$Time.y == Rsquared_10pct$Time.y[row] &
																						summary_SCT_10pct_wide_wide$Time.x == Rsquared_10pct$Time.x[row]),]))$r.squared
	Rsquared_10pct$pearson[row] <- cor(summary_SCT_10pct_wide_wide[which(summary_SCT_10pct_wide_wide$Line == Rsquared_10pct$Line[row] & 
																						summary_SCT_10pct_wide_wide$Location == Rsquared_10pct$Location[row] & 
																						summary_SCT_10pct_wide_wide$Time.y == Rsquared_10pct$Time.y[row] &
																						summary_SCT_10pct_wide_wide$Time.x == Rsquared_10pct$Time.x[row]),]$Mean.y, 
							summary_SCT_10pct_wide_wide[which(summary_SCT_10pct_wide_wide$Line == Rsquared_10pct$Line[row] & 
																						summary_SCT_10pct_wide_wide$Location == Rsquared_10pct$Location[row] & 
																						summary_SCT_10pct_wide_wide$Time.y == Rsquared_10pct$Time.y[row] &
																						summary_SCT_10pct_wide_wide$Time.x == Rsquared_10pct$Time.x[row]),]$Mean.x)
	Rsquared_10pct$spearman[row] <- cor(summary_SCT_10pct_wide_wide[which(summary_SCT_10pct_wide_wide$Line == Rsquared_10pct$Line[row] & 
																						summary_SCT_10pct_wide_wide$Location == Rsquared_10pct$Location[row] & 
																						summary_SCT_10pct_wide_wide$Time.y == Rsquared_10pct$Time.y[row] &
																						summary_SCT_10pct_wide_wide$Time.x == Rsquared_10pct$Time.x[row]),]$Mean.y, 
							summary_SCT_10pct_wide_wide[which(summary_SCT_10pct_wide_wide$Line == Rsquared_10pct$Line[row] & 
																						summary_SCT_10pct_wide_wide$Location == Rsquared_10pct$Location[row] & 
																						summary_SCT_10pct_wide_wide$Time.y == Rsquared_10pct$Time.y[row] &
																						summary_SCT_10pct_wide_wide$Time.x == Rsquared_10pct$Time.x[row]),]$Mean.x, method = "spearman")
}

Rsquared_10pct
Rsquared_1pct
Rsquared





