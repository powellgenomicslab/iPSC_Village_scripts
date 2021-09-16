library(ggplot2)
library(tidyverse)
library(Seurat)
library(lme4)
library(hier.part)
library(scater)
library(data.table)
library(Seurat)
library(gridExtra)
library(ggpmisc)


save_figs <- function(plot, basename, width = 17, height = 17, units = "cm"){
    ggsave(plot, filename = paste0(basename,".png"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".pdf"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".eps"), height = height, width = width, units = units)
}



##### Set up colors #####
variable_colors <- c(Village = "#A2B0D0", Replicate = "#64A66B", Line = "#68319B") 
line_colors <- c(FSA0006 = "#F79E29", MBE1006 = "#9B2C99", TOB0421 = "#35369C")


##### Set Up Directories #####
dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
outdir <- paste0(dir,"output/Variance_contributions_lm_scale_cluster/")
datadir <- paste0(dir,"output/All_data_integrated_remove_bad/")

dir.create(outdir)



seurat <- readRDS(paste0(datadir, "seurat_integrated_all_times_clustered_rb_regressed.rds"))
seurat$Location <- gsub("_.+", "", seurat$Location_Time)
seurat$Location <- ifelse((seurat$Time == "Thawed Village Day 0" | seurat$Time == "Thawed Village Day 7"), "Sydney-Cryopreserved", ifelse(((seurat$Location == "Sydney") & (seurat$Time == "Baseline" | seurat$Time == "Village Day 4")), "Sydney-Fresh", seurat$Location))


seurat_list <- lapply(unique(seurat@meta.data$Location), function(x){
	temp <- lapply(unique(Idents(seurat)), function(y){
		subset(seurat, subset = Location == x, idents = y)
	})
	names(temp) <- unique(Idents(seurat))
	return(temp)
})
names(seurat_list) <- unique(seurat@meta.data$Location)



#### Filter out Genes expressed in low percentage of cells #####
pct_expressing <- lapply(seurat_list, function(x){
	lapply(x, function(y){
		temp <- data.frame(rowSums(as.matrix(y[["SCT"]]@counts) > 0)/ncol(y[["SCT"]]@counts))
		colnames(temp) <- "Cell_proportion"
		return(temp)
	})
})

### 10% to start
seurat_list <- lapply(names(seurat_list), function(x){
	temp <- lapply(names(seurat_list[[x]]), function(y){
		subset(seurat_list[[x]][[y]], features = rownames(pct_expressing[[x]][[y]])[which(pct_expressing[[x]][[y]]$Cell_proportion >= 0.1)])
	})
	names(temp) <- unique(Idents(seurat))
	return(temp)
})
names(seurat_list) <- unique(seurat@meta.data$Location)




if (!file.exists(paste0(outdir, "heirarchical_partitioning_results_clusters.rds"))){
	hier_models <- list()
	i <- 1
	for(x in names(seurat_list)){
		j <- 1
		for (y in names(seurat_list[[x]])){
			k <- 1
			for (gene in unique(rownames(seurat_list[[x]][[y]]))){
				print(paste(i,j,k))
				df <- data.frame("Expression" = data.frame(seurat_list[[x]][[y]][["SCT"]]@scale.data[gene,]), "Village" = ifelse((seurat_list[[x]][[y]]@meta.data$Time == "Baseline" | seurat_list[[x]][[y]]@meta.data$Time == "Thawed Village Day 0"), 0, 1), "Line" = seurat_list[[x]][[y]]@meta.data$Final_Assignment, "Replicate" = gsub("[A-Z][a-z]+", "", seurat_list[[x]][[y]]@meta.data$MULTI_ID))
				colnames(df)[1] <- "Expression"
				hier_models[[x]][[y]][[gene]] <- hier.part(df$Expression, df[,2:ncol(df)])
				k <- k + 1
			}
			j <- j + 1
		}
		i <- 1 + i
	}
	saveRDS(hier_models, paste0(outdir, "heirarchical_partitioning_results_clusters.rds"))
} else {
	hier_models <- readRDS(paste0(outdir, "heirarchical_partitioning_results_clusters.rds"))
}



##### Make dataframes with all the values
contributions_df_list <- lapply(names(hier_models), function(x){
	temp3 <- lapply(names(hier_models[[x]]), function(y){
		temp2 <- lapply(names(hier_models[[x]][[y]]), function(z){
			temp <-hier_models[[x]][[y]][[z]]$IJ
			temp$Variable <- rownames(temp)
			temp$Gene <- z
			temp$Location <- x
			temp$Cluster <- y
			temp$Independendent_Pct_Explained <- hier_models[[x]][[y]][[z]]$I.perc$ind.exp.var
			colnames(temp) <- c("Independent_Contribution", "Joint_Contribution", "Total_Contribution", "Variable", "Gene", "Location", "Cluster", "Independendent_Pct_Explained")
			rownames(temp) <- paste0(temp$Variable, "_", temp$Gene, "_", temp$Location)
			return(temp)
		})
		out_list <- do.call(rbind,temp2)
		return(out_list)
	})
	out_list2 <- do.call(rbind,temp3)
	return(out_list2)
})

contributions_df <- do.call(rbind, contributions_df_list)
contributions_df$Location <- factor(contributions_df$Location, levels = c("Brisbane", "Melbourne", "Sydney-Fresh", "Sydney-Cryopreserved"))



##### Make some figures!!! #####
pTotal_Cont <- ggplot(contributions_df, aes(Total_Contribution, color = Variable)) +
	geom_density() +
	facet_grid(Location ~ Cluster) +
	theme_classic()  +
	# scale_y_continuous(trans = "log10") +
	scale_x_continuous(trans = "log10") +
	scale_color_manual(values = variable_colors) +
	geom_vline(xintercept = 0.01, linetype="dashed")

save_figs(pTotal_Cont, paste0(outdir, "Total_Contribution_Histogram"), width = 40, height = 20)


### Variance explained correlation between sites ###
##### Correlate the results #####
contributions_df_joined <- full_join(contributions_df, contributions_df, by = c("Variable", "Gene", "Cluster"))
contributions_df_joined <- contributions_df_joined[which(contributions_df_joined$Location.x != contributions_df_joined$Location.y),]

Rsquared_list <- lapply(unique(contributions_df_joined$Cluster), function(x){
	temp <- unique(contributions_df_joined[,c("Location.x","Location.y")])
	temp$Line_Rsquared <- NA
	temp$Replicate_Rsquared <- NA
	temp$Village_Rsquared <- NA
	for (row in 1:nrow(temp)){
		temp$Line_Rsquared[row] <- summary(lm(Total_Contribution.y ~ Total_Contribution.x, contributions_df_joined[which(contributions_df_joined$Variable == "Line" & contributions_df_joined$Location.x == temp$Location.x[row] & contributions_df_joined$Location.y == temp$Location.y[row]),]))$r.squared
		temp$Replicate_Rsquared[row] <- summary(lm(Total_Contribution.y ~ Total_Contribution.x, contributions_df_joined[which(contributions_df_joined$Variable == "Replicate" & contributions_df_joined$Location.x == temp$Location.x[row] & contributions_df_joined$Location.y == temp$Location.y[row]),]))$r.squared
		temp$Village_Rsquared[row] <- summary(lm(Total_Contribution.y ~ Total_Contribution.x, contributions_df_joined[which(contributions_df_joined$Variable == "Village" & contributions_df_joined$Location.x == temp$Location.x[row] & contributions_df_joined$Location.y == temp$Location.y[row]),]))$r.squared
	}
	return(temp)
}) 
names(Rsquared_list) <- unique(contributions_df_joined$Cluster)


Rsquared_long_list <- lapply(Rsquared_list, function(x){
	temp <- pivot_longer(x, cols = c( "Line_Rsquared", "Replicate_Rsquared", "Village_Rsquared"), values_to = "Rsquared", names_to = "Variable")
	temp$Variable <- gsub("_Rsquared", "", temp$Variable)
	temp$Rsquared <- as.character(round(temp$Rsquared, 2))
	return(temp)
})



##### Figures #####
pCorr_Location_Point <- list()
for (cluster in names(Rsquared_long_list)){
	pCorr_Location_Point[[cluster]] <- ggplot(contributions_df_joined[which(contributions_df_joined$Cluster == cluster),], aes(Total_Contribution.x, Total_Contribution.y, color = Variable)) +
		geom_point(size = 0.5, alpha = 0.25) +
		facet_grid(rows = vars(Location.y), cols = vars(Location.x)) +
		theme_classic() +
		scale_color_manual(values = variable_colors) +
		xlab("Contribution to Gene Variation") +
		ylab("Contribution to Gene Variation") +
		geom_text(x = 0, y = 0.75, aes(label = paste0("R^2 = ", Rsquared)), data = Rsquared_long_list[[cluster]][which(Rsquared_long_list[[cluster]]$Variable == "Line"),], hjust = 0) +
		geom_text(x = 0, y = 0.7, aes(label = paste0("R^2 = ", Rsquared)), data = Rsquared_long_list[[cluster]][which(Rsquared_long_list[[cluster]]$Variable == "Village"),], hjust = 0) +
		geom_text(x = 0, y = 0.65, aes(label = paste0("R^2 = ", Rsquared)), data = Rsquared_long_list[[cluster]][which(Rsquared_long_list[[cluster]]$Variable == "Replicate"),], hjust = 0)

	save_figs(pCorr_Location_Point[[cluster]], paste0(outdir, "Correlation_Point_Location_combined_cluster",cluster), width = 25, height = 20)
}




##### RERUN WITH 1/N CELLS AS COVARIATE TO SEE EFFECT OF NUMBER OF CELLS #####
seurat_list <- lapply(seurat_list, function(x){
	lapply(x, function(y){
		y$Ncorrection <- NA
		for (rep in unique(y$Site_rep)){
			print(rep)
			y$Ncorrection <- ifelse(y$Site_rep == rep, 1/nrow(y@meta.data[which(y$Site_rep == rep),]), y$Ncorrection)
		}
		return(y)
	})
})



###
if (!file.exists(paste0(outdir, "heirarchical_partitioning_results_clusters_Ncorrected.rds"))){
	hier_models_Ncorrect <- list()
	i <- 1
	for(x in names(seurat_list)){
		j <- 1
		for (y in names(seurat_list[[x]])){
			k <- 1
			for (gene in unique(rownames(seurat_list[[x]][[y]]))){
				print(paste(i,j,k))
				df <- data.frame("Expression" = data.frame(seurat_list[[x]][[y]][["SCT"]]@scale.data[gene,]), "Village" = ifelse((seurat_list[[x]][[y]]@meta.data$Time == "Baseline" | seurat_list[[x]][[y]]@meta.data$Time == "Thawed Village Day 0"), 0, 1), "Line" = seurat_list[[x]][[y]]@meta.data$Final_Assignment, "Replicate" = gsub("[A-Z][a-z]+", "", seurat_list[[x]][[y]]@meta.data$MULTI_ID), "Ncorrection" = seurat_list[[x]][[y]]@meta.data$Ncorrection)
				colnames(df)[1] <- "Expression"
				hier_models_Ncorrect[[x]][[y]][[gene]] <- hier.part(df$Expression, df[,2:ncol(df)])
				k <- k + 1
			}
			j <- j + 1
		}
		i <- 1 + i
	}
	saveRDS(hier_models_Ncorrect, paste0(outdir, "heirarchical_partitioning_results_clusters_Ncorrected.rds"))
} else {
	hier_models_Ncorrect <- readRDS(paste0(outdir, "heirarchical_partitioning_results_clusters_Ncorrected.rds"))
}



##### Make dataframes with all the values
contributions_df_Ncorrect_list <- lapply(names(hier_models_Ncorrect), function(x){
	temp3 <- lapply(names(hier_models_Ncorrect[[x]]), function(y){
		temp2 <- lapply(names(hier_models_Ncorrect[[x]][[y]]), function(z){
			temp <-hier_models_Ncorrect[[x]][[y]][[z]]$IJ
			temp$Variable <- rownames(temp)
			temp$Gene <- z
			temp$Location <- x
			temp$Cluster <- y
			temp$Independendent_Pct_Explained <- hier_models_Ncorrect[[x]][[y]][[z]]$I.perc$ind.exp.var
			colnames(temp) <- c("Independent_Contribution", "Joint_Contribution", "Total_Contribution", "Variable", "Gene", "Location", "Cluster", "Independendent_Pct_Explained")
			rownames(temp) <- paste0(temp$Variable, "_", temp$Gene, "_", temp$Location)
			return(temp)
		})
		out_list <- do.call(rbind,temp2)
		return(out_list)
	})
	out_list2 <- do.call(rbind,temp3)
	return(out_list2)
})

contributions_Ncorrect_df <- do.call(rbind, contributions_df_Ncorrect_list)
contributions_Ncorrect_df$Location <- factor(contributions_Ncorrect_df$Location, levels = c("Brisbane", "Melbourne", "Sydney-Fresh", "Sydney-Cryopreserved"))



##### Make some figures!!! #####
variable_colors_new <- c(variable_colors, "Ncorrection" = "darkorange3")
pTotal_Cont_Ncorrect <- ggplot(contributions_Ncorrect_df, aes(Total_Contribution, color = Variable)) +
	geom_density() +
	facet_grid(Location ~ Cluster) +
	theme_classic()  +
	# scale_y_continuous(trans = "log10") +
	scale_x_continuous(trans = "log10") +
	scale_color_manual(values = variable_colors_new) +
	geom_vline(xintercept = 0.01, linetype="dashed")

save_figs(pTotal_Cont_Ncorrect, paste0(outdir, "Total_Contribution_Histogram_Ncorrect"), width = 40, height = 20)



### Variance explained correlation between sites ###
##### Correlate the results #####
contributions_Ncorrect_df_joined <- full_join(contributions_Ncorrect_df, contributions_Ncorrect_df, by = c("Variable", "Gene", "Cluster"))
contributions_Ncorrect_df_joined <- contributions_Ncorrect_df_joined[which(contributions_Ncorrect_df_joined$Location.x != contributions_Ncorrect_df_joined$Location.y),]

Rsquared_Ncorrect_list <- lapply(unique(contributions_Ncorrect_df_joined$Cluster), function(x){
	temp <- unique(contributions_Ncorrect_df_joined[,c("Location.x","Location.y")])
	temp$Line_Rsquared <- NA
	temp$Replicate_Rsquared <- NA
	temp$Village_Rsquared <- NA
	temp$Ncorrect_Rsquared <- NA
	for (row in 1:nrow(temp)){
		temp$Line_Rsquared[row] <- summary(lm(Total_Contribution.y ~ Total_Contribution.x, contributions_Ncorrect_df_joined[which(contributions_Ncorrect_df_joined$Variable == "Line" & contributions_Ncorrect_df_joined$Location.x == temp$Location.x[row] & contributions_Ncorrect_df_joined$Location.y == temp$Location.y[row]),]))$r.squared
		temp$Replicate_Rsquared[row] <- summary(lm(Total_Contribution.y ~ Total_Contribution.x, contributions_Ncorrect_df_joined[which(contributions_Ncorrect_df_joined$Variable == "Replicate" & contributions_Ncorrect_df_joined$Location.x == temp$Location.x[row] & contributions_Ncorrect_df_joined$Location.y == temp$Location.y[row]),]))$r.squared
		temp$Village_Rsquared[row] <- summary(lm(Total_Contribution.y ~ Total_Contribution.x, contributions_Ncorrect_df_joined[which(contributions_Ncorrect_df_joined$Variable == "Village" & contributions_Ncorrect_df_joined$Location.x == temp$Location.x[row] & contributions_Ncorrect_df_joined$Location.y == temp$Location.y[row]),]))$r.squared
		temp$Ncorrect_Rsquared[row] <- summary(lm(Total_Contribution.y ~ Total_Contribution.x, contributions_Ncorrect_df_joined[which(contributions_Ncorrect_df_joined$Variable == "Ncorrection" & contributions_Ncorrect_df_joined$Location.x == temp$Location.x[row] & contributions_Ncorrect_df_joined$Location.y == temp$Location.y[row]),]))$r.squared
	}
	return(temp)
}) 
names(Rsquared_Ncorrect_list) <- unique(contributions_Ncorrect_df_joined$Cluster)


Rsquared_Ncorrect_long_list <- lapply(Rsquared_Ncorrect_list, function(x){
	temp <- pivot_longer(x, cols = c( "Line_Rsquared", "Replicate_Rsquared", "Village_Rsquared"), values_to = "Rsquared", names_to = "Variable")
	temp$Variable <- gsub("_Rsquared", "", temp$Variable)
	temp$Rsquared <- as.character(round(temp$Rsquared, 2))
	return(temp)
})



##### Figures #####
pCorr_Location_Point <- list()
for (cluster in names(Rsquared_Ncorrect_long_list)){
	pCorr_Location_Point[[cluster]] <- ggplot(contributions_Ncorrect_df_joined[which(contributions_Ncorrect_df_joined$Cluster == cluster),], aes(Total_Contribution.x, Total_Contribution.y, color = Variable)) +
		geom_point(size = 0.5, alpha = 0.25) +
		facet_grid(rows = vars(Location.y), cols = vars(Location.x)) +
		theme_classic() +
		scale_color_manual(values = variable_colors_new) +
		xlab("Contribution to Gene Variation") +
		ylab("Contribution to Gene Variation") +
		geom_text(x = 0, y = 0.825, aes(label = paste0("R^2 = ", Rsquared)), data = Rsquared_Ncorrect_long_list[[cluster]][which(Rsquared_Ncorrect_long_list[[cluster]]$Variable == "Line"),], hjust = 0) +
		geom_text(x = 0, y = 0.775, aes(label = paste0("R^2 = ", Rsquared)), data = Rsquared_Ncorrect_long_list[[cluster]][which(Rsquared_Ncorrect_long_list[[cluster]]$Variable == "Village"),], hjust = 0) +
		geom_text(x = 0, y = 0.625, aes(label = paste0("R^2 = ", Rsquared)), data = Rsquared_Ncorrect_long_list[[cluster]][which(Rsquared_Ncorrect_long_list[[cluster]]$Variable == "Replicate"),], hjust = 0)
		geom_text(x = 0, y = 0.65, aes(label = paste0("R^2 = ", Rsquared)), data = Rsquared_Ncorrect_long_list[[cluster]][which(Rsquared_Ncorrect_long_list[[cluster]]$Variable == "Ncorrection"),], hjust = 0)

	save_figs(pCorr_Location_Point[[cluster]], paste0(outdir, "Correlation_Point_Location_combined_cluster",cluster), width = 25, height = 20)
}
























