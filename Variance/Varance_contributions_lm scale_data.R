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
library(grDevices)
library(plyr)


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
outdir <- paste0(dir,"output/Variance_contributions_lm_scale_data/")
datadir <- paste0(dir,"output/All_data_integrated_remove_bad/")

dir.create(outdir)


##### Read in Data #####
seurat <- readRDS(paste0(datadir,"seurat_integrated_all_times_clustered.rds"))
seurat_sub <- subset(seurat, subset = Time != "Thawed Village Day 0")
seurat_sub <- subset(seurat_sub, subset = Time != "Thawed Village Day 7")
seurat_sub$Location <- gsub("_.+", "", seurat_sub$Location_Time)

seurat_list <- lapply(unique(seurat_sub@meta.data$Location), function(x){
	subset(seurat_sub, subset = Location == x)
})
names(seurat_list) <- unique(seurat_sub@meta.data$Location)


#### Filter out Genes expressed in low percentage of cells #####
pct_expressing <- lapply(seurat_list, function(x){
	temp <- data.frame(rowSums(as.matrix(x[["SCT"]]@counts) > 0)/ncol(x[["SCT"]]@counts))
	colnames(temp) <- "Cell_proportion"
	return(temp)
})

### 10% to start
seurat_list <- lapply(names(seurat_list), function(x){
	subset(seurat_list[[x]], features = rownames(pct_expressing[[x]])[which(pct_expressing[[x]]$Cell_proportion >= 0.1)])
})
names(seurat_list) <- unique(seurat_sub@meta.data$Location)


if (!file.exists(paste0(outdir, "heirarchical_partitioning_results.rds"))){
	hier_models <- list()
	i <- 1
	for(x in names(seurat_list)){
		j <- 1
		for (gene in unique(rownames(seurat_list[[x]]))){
			print(paste(i,j))
			df <- data.frame("Expression" = data.frame(seurat_list[[x]][["SCT"]]@scale.data[gene,]), "Village" = ifelse(seurat_list[[x]]@meta.data$Time == "Baseline", 0, 1), "Line" = seurat_list[[x]]@meta.data$Final_Assignment, "Replicate" = gsub("[A-Z][a-z]+", "", seurat_list[[x]]@meta.data$MULTI_ID))
			colnames(df)[1] <- "Expression"
			hier_models[[x]][[gene]] <- hier.part(df$Expression, df[,2:ncol(df)])
			j <- j + 1
		}
		i <- 1 + i
	}

	saveRDS(hier_models, paste0(outdir, "heirarchical_partitioning_results.rds"))
} else {
	hier_models <- readRDS(paste0(outdir, "heirarchical_partitioning_results.rds"))
}


##### Make dataframes with all the values
joint_df <- data.frame(matrix(nrow = 0, ncol = 3))
colnames("Gene", "Variable", "Variance")
independent_df
total_df
independent_perc_df

temp <- list()


contributions_df_list <- lapply(names(hier_models), function(x){
	temp2 <- lapply(names(hier_models[[x]]), function(y){
		temp <-hier_models[[x]][[y]]$IJ
		temp$Variable <- rownames(temp)
		temp$Gene <- y
		temp$Location <- x
		temp$Independendent_Pct_Explained <- hier_models[[x]][[y]]$I.perc$ind.exp.var
		colnames(temp) <- c("Independent_Contribution", "Joint_Contribution", "Total_Contribution", "Variable", "Gene", "Location", "Independendent_Pct_Explained")
		rownames(temp) <- paste0(temp$Variable, "_", temp$Gene, "_", temp$Location)
		return(temp)
	})
	print(tail(temp2))
	print(length(temp2))
	out_list <- do.call(rbind,temp2)
	return(out_list)
})

contributions_df <- do.call(rbind, contributions_df_list)



##### Make some figures!!! #####
pTotal_Cont <- ggplot(contributions_df, aes(Total_Contribution, color = Variable)) +
	geom_density() +
	facet_wrap(Location ~ ., ncol = 1) +
	theme_classic()  +
	# scale_y_continuous(trans = "log10") +
	scale_x_continuous(trans = "log10")

save_figs(pTotal_Cont, paste0(outdir, "Total_Contribution_Histogram"))



pIndep_Cont <- ggplot(contributions_df, aes(Independent_Contribution, color = Variable)) +
	geom_density() +
	facet_wrap(Location ~ ., ncol = 1) +
	theme_classic()  +
	# scale_y_continuous(trans = "log10") +
	scale_x_continuous(trans = "log10")

save_figs(pIndep_Cont, paste0(outdir, "Independent_Contribution_Histogram"))


pJoint_Cont <- ggplot(contributions_df, aes(Joint_Contribution, color = Variable)) +
	geom_density() +
	facet_wrap(Location ~ ., ncol = 1) +
	theme_classic()  +
	# scale_y_continuous(trans = "log10") +
	scale_x_continuous(trans = "log10")

save_figs(pJoint_Cont, paste0(outdir, "Joint_Contribution_Histogram"))


pIndep_pct_Cont <- ggplot(contributions_df, aes(Independendent_Pct_Explained, color = Variable)) +
	geom_density() +
	facet_wrap(Location ~ ., ncol = 1) +
	theme_classic()  +
	# scale_y_continuous(trans = "log10") +
	scale_x_continuous(trans = "log10")

save_figs(pIndep_pct_Cont, paste0(outdir, "Independendent_Pct_Explained_Histogram"))



### Variance explained correlation between sites ###
contributions_df_joined <- full_join(contributions_df, contributions_df, by = c("Variable", "Gene"))
pCorr_Location_Point <- ggplot(contributions_df_joined, aes(Independent_Contribution.x, Independent_Contribution.y, color = Variable)) +
	geom_point(size = 0.5, alpha = 0.25) +
	facet_grid(rows = vars(Location.y), cols = vars(Location.x)) +
	theme_classic() +
	scale_y_continuous(trans = "log10") +
	scale_x_continuous(trans = "log10")
save_figs(pCorr_Location_Point, paste0(outdir, "Correlation_Point_Location"))


pCorr_Location_Point_Line <- ggplot(contributions_df_joined[which(contributions_df_joined$Variable == "Line"),], aes(Independent_Contribution.x, Independent_Contribution.y)) +
	geom_point(size = 0.5, alpha = 0.25) +
	facet_grid(rows = vars(Location.y), cols = vars(Location.x)) +
	theme_classic() +
	scale_y_continuous(trans = "log10") +
	scale_x_continuous(trans = "log10")
save_figs(pCorr_Location_Point_Line, paste0(outdir, "Correlation_Point_Location_Line"))


pCorr_Location_Point_Village <- ggplot(contributions_df_joined[which(contributions_df_joined$Variable == "Village"),], aes(Independent_Contribution.x, Independent_Contribution.y)) +
	geom_point(size = 0.5, alpha = 0.25) +
	facet_grid(rows = vars(Location.y), cols = vars(Location.x)) +
	theme_classic() +
	scale_y_continuous(trans = "log10") +
	scale_x_continuous(trans = "log10")
save_figs(pCorr_Location_Point_Village, paste0(outdir, "Correlation_Point_Location_Village"))


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


contributions_df_pluri_genes <- left_join(pluri_genes, contributions_df, by = c("ENSG_ID" = "Gene"))
contributions_df_pluri_genes <- na.omit(contributions_df_pluri_genes)


pPluri_Genes_Cont <- ggplot(contributions_df_pluri_genes, aes(Location, Independent_Contribution, fill = Variable)) +
						geom_bar(position = "stack", stat = "identity") +
						theme_classic() +
						facet_wrap(Gene_ID ~ ., nrow = 1)
save_figs(pPluri_Genes_Cont, paste0(outdir, "Pluripotent_Gene_Variable_Contributions"), width = 50, height = 10)






##### Redo with rb % regressed out #####
##### Read in Data #####
seurat_rb <- readRDS(paste0(datadir,"seurat_integrated_all_times_clustered_rb_regressed.rds"))
seurat_rb_sub <- subset(seurat_rb, subset = Time != "Thawed Village Day 0")
seurat_rb_sub <- subset(seurat_rb_sub, subset = Time != "Thawed Village Day 7")
seurat_rb_sub$Location <- gsub("_.+", "", seurat_rb_sub$Location_Time)

seurat_rb_list <- lapply(unique(seurat_rb_sub@meta.data$Location), function(x){
	subset(seurat_rb_sub, subset = Location == x)
})
names(seurat_rb_list) <- unique(seurat_rb_sub@meta.data$Location)


#### Filter out Genes expressed in low percentage of cells #####
pct_expressing_rb <- lapply(seurat_rb_list, function(x){
	temp <- data.frame(rowSums(as.matrix(x[["SCT"]]@counts) > 0)/ncol(x[["SCT"]]@counts))
	colnames(temp) <- "Cell_proportion"
	return(temp)
})

### 10% to start
seurat_rb_list <- lapply(names(seurat_rb_list), function(x){
	subset(seurat_rb_list[[x]], features = rownames(pct_expressing_rb[[x]])[which(pct_expressing_rb[[x]]$Cell_proportion >= 0.1)])
})
names(seurat_rb_list) <- unique(seurat_rb_sub@meta.data$Location)


if (!file.exists(paste0(outdir, "heirarchical_partitioning_results_rb_regressed.rds"))){
	hier_models_rb <- list()
	i <- 1
	for(x in names(seurat_rb_list)){
		j <- 1
		for (gene in unique(rownames(seurat_rb_list[[x]]))){
			print(paste(i,j))
			df <- data.frame("Expression" = data.frame(seurat_rb_list[[x]][["SCT"]]@scale.data[gene,]), "Village" = ifelse(seurat_rb_list[[x]]@meta.data$Time == "Baseline", 0, 1), "Line" = seurat_rb_list[[x]]@meta.data$Final_Assignment, "Replicate" = gsub("[A-Z][a-z]+", "", seurat_rb_list[[x]]@meta.data$MULTI_ID))
			colnames(df)[1] <- "Expression"
			hier_models_rb[[x]][[gene]] <- hier.part(df$Expression, df[,2:ncol(df)])
			j <- j + 1
		}
		i <- 1 + i
	}

	saveRDS(hier_models_rb, paste0(outdir, "heirarchical_partitioning_results_rb_regressed.rds"))
} else {
	hier_models_rb <- readRDS(paste0(outdir, "heirarchical_partitioning_results_rb_regressed.rds"))
}



##### Make dataframes with all the values
contributions_rb_regress_df_list <- lapply(names(hier_models_rb), function(x){
	temp2 <- lapply(names(hier_models_rb[[x]]), function(y){
		temp <-hier_models_rb[[x]][[y]]$IJ
		temp$Variable <- rownames(temp)
		temp$Gene <- y
		temp$Location <- x
		temp$Independendent_Pct_Explained <- hier_models_rb[[x]][[y]]$I.perc$ind.exp.var
		colnames(temp) <- c("Independent_Contribution", "Joint_Contribution", "Total_Contribution", "Variable", "Gene", "Location", "Independendent_Pct_Explained")
		rownames(temp) <- paste0(temp$Variable, "_", temp$Gene, "_", temp$Location)
		return(temp)
	})
	print(tail(temp2))
	print(length(temp2))
	out_list <- do.call(rbind,temp2)
	return(out_list)
})

contributions_rb_regress_df <- do.call(rbind, contributions_rb_regress_df_list)



##### Make some figures!!! #####
pTotal_Cont <- ggplot(contributions_rb_regress_df, aes(Total_Contribution, color = Variable)) +
	geom_density() +
	facet_wrap(Location ~ ., ncol = 1) +
	theme_classic()  +
	# scale_y_continuous(trans = "log10") +
	scale_x_continuous(trans = "log10")

save_figs(pTotal_Cont, paste0(outdir, "Total_Contribution_Histogram_rb_regressed"))




pIndep_Cont <- ggplot(contributions_rb_regress_df, aes(Independent_Contribution, color = Variable)) +
	geom_density() +
	facet_wrap(Location ~ ., ncol = 1) +
	theme_classic()  +
	# scale_y_continuous(trans = "log10") +
	scale_x_continuous(trans = "log10")

save_figs(pIndep_Cont, paste0(outdir, "Independent_Contribution_Histogram_rb_regressed"))


pJoint_Cont <- ggplot(contributions_rb_regress_df, aes(Joint_Contribution, color = Variable)) +
	geom_density() +
	facet_wrap(Location ~ ., ncol = 1) +
	theme_classic()  +
	# scale_y_continuous(trans = "log10") +
	scale_x_continuous(trans = "log10")

save_figs(pJoint_Cont, paste0(outdir, "Joint_Contribution_Histogram_rb_regressed"))


pIndep_pct_Cont <- ggplot(contributions_rb_regress_df, aes(Independendent_Pct_Explained, color = Variable)) +
	geom_density() +
	facet_wrap(Location ~ ., ncol = 1) +
	theme_classic()  +
	# scale_y_continuous(trans = "log10") +
	scale_x_continuous(trans = "log10")

save_figs(pIndep_pct_Cont, paste0(outdir, "Independendent_Pct_Explained_Histogram_rb_regressed"))



### Variance explained correlation between sites ###
contributions_rb_regress_df_joined <- full_join(contributions_rb_regress_df, contributions_rb_regress_df, by = c("Variable", "Gene"))
pCorr_Location_Point <- ggplot(contributions_rb_regress_df_joined, aes(Independent_Contribution.x, Independent_Contribution.y, color = Variable)) +
	geom_point(size = 0.5, alpha = 0.25) +
	facet_grid(rows = vars(Location.y), cols = vars(Location.x)) +
	theme_classic() +
	scale_y_continuous(trans = "log10") +
	scale_x_continuous(trans = "log10")
save_figs(pCorr_Location_Point, paste0(outdir, "Correlation_Point_Location_rb_regressed"))


pCorr_Location_Point_Line <- ggplot(contributions_rb_regress_df_joined[which(contributions_rb_regress_df_joined$Variable == "Line"),], aes(Independent_Contribution.x, Independent_Contribution.y)) +
	geom_point(size = 0.5, alpha = 0.25) +
	facet_grid(rows = vars(Location.y), cols = vars(Location.x)) +
	theme_classic() +
	scale_y_continuous(trans = "log10") +
	scale_x_continuous(trans = "log10")
save_figs(pCorr_Location_Point_Line, paste0(outdir, "Correlation_Point_Location_Line_rb_regressed"))


pCorr_Location_Point_Village <- ggplot(contributions_rb_regress_df_joined[which(contributions_rb_regress_df_joined$Variable == "Village"),], aes(Independent_Contribution.x, Independent_Contribution.y)) +
	geom_point(size = 0.5, alpha = 0.25) +
	facet_grid(rows = vars(Location.y), cols = vars(Location.x)) +
	theme_classic() +
	scale_y_continuous(trans = "log10") +
	scale_x_continuous(trans = "log10")
save_figs(pCorr_Location_Point_Village, paste0(outdir, "Correlation_Point_Location_Village_rb_regressed"))


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


contributions_rb_regress_df_pluri_genes <- left_join(pluri_genes, contributions_rb_regress_df, by = c("ENSG_ID" = "Gene"))
contributions_rb_regress_df_pluri_genes <- na.omit(contributions_rb_regress_df_pluri_genes)


pPluri_Genes_Cont <- ggplot(contributions_rb_regress_df_pluri_genes, aes(Location, Independent_Contribution, fill = Variable)) +
						geom_bar(position = "stack", stat = "identity") +
						theme_classic() +
						facet_wrap(Gene_ID ~ ., nrow = 1)
save_figs(pPluri_Genes_Cont, paste0(outdir, "Pluripotent_Gene_Variable_Contributions_rb_regressed"), width = 50, height = 10)


##### Surprisingly, the village effect seems to be somewhat variable between the locations - want to try normalizing without the freeze-thaw data #####
### Previously did and saw that regressing rb% did help
### Promising = the line effect seems to be more consistent than the village effect 
if (!file.exists(paste0(outdir, "seurat_SCT_wRB_base_village_only.rds"))){
	seurat_rb_sub_norm <- SCTransform(seurat_rb_sub, verbose = TRUE, vars.to.regress = c("scores.G1", "scores.S", "scores.G2M", "percent.mt", "percent.rb"), return.only.var.genes = FALSE)
	saveRDS(seurat_rb_sub_norm, paste0(outdir, "seurat_SCT_wRB_base_village_only.rds"))
} else {
	seurat_rb_sub_norm <- readRDS(paste0(outdir, "seurat_SCT_wRB_base_village_only.rds"))
}



seurat_rb_norm_list <- lapply(unique(seurat_rb_sub_norm@meta.data$Location), function(x){
	subset(seurat_rb_sub_norm, subset = Location == x)
})
names(seurat_rb_norm_list) <- unique(seurat_rb_sub_norm@meta.data$Location)


#### Filter out Genes expressed in low percentage of cells #####
pct_expressing_rb <- lapply(seurat_rb_norm_list, function(x){
	temp <- data.frame(rowSums(as.matrix(x[["SCT"]]@counts) > 0)/ncol(x[["SCT"]]@counts))
	colnames(temp) <- "Cell_proportion"
	return(temp)
})

### 10% to start
seurat_rb_norm_list <- lapply(names(seurat_rb_norm_list), function(x){
	subset(seurat_rb_norm_list[[x]], features = rownames(pct_expressing_rb[[x]])[which(pct_expressing_rb[[x]]$Cell_proportion >= 0.1)])
})
names(seurat_rb_norm_list) <- unique(seurat_rb_sub_norm@meta.data$Location)


if (!file.exists(paste0(outdir, "heirarchical_partitioning_results_rb_regressed_norm.rds"))){
	hier_models_rb <- list()
	i <- 1
	for(x in names(seurat_rb_norm_list)){
		j <- 1
		for (gene in unique(rownames(seurat_rb_norm_list[[x]]))){
			print(paste(i,j))
			df <- data.frame("Expression" = data.frame(seurat_rb_norm_list[[x]][["SCT"]]@scale.data[gene,]), "Village" = ifelse(seurat_rb_norm_list[[x]]@meta.data$Time == "Baseline", 0, 1), "Line" = seurat_rb_norm_list[[x]]@meta.data$Final_Assignment, "Replicate" = gsub("[A-Z][a-z]+", "", seurat_rb_norm_list[[x]]@meta.data$MULTI_ID))
			colnames(df)[1] <- "Expression"
			hier_models_rb[[x]][[gene]] <- hier.part(df$Expression, df[,2:ncol(df)])
			j <- j + 1
		}
		i <- 1 + i
	}
	saveRDS(hier_models_rb, paste0(outdir, "heirarchical_partitioning_results_rb_regressed_norm.rds"))
} else {
	hier_models_rb <- readRDS(paste0(outdir, "heirarchical_partitioning_results_rb_regressed_norm.rds"))
}




##### Make dataframes with all the values
contributions_rb_regress_norm_df_list <- lapply(names(hier_models_rb), function(x){
	temp2 <- lapply(names(hier_models_rb[[x]]), function(y){
		temp <-hier_models_rb[[x]][[y]]$IJ
		temp$Variable <- rownames(temp)
		temp$Gene <- y
		temp$Location <- x
		temp$Independendent_Pct_Explained <- hier_models_rb[[x]][[y]]$I.perc$ind.exp.var
		colnames(temp) <- c("Independent_Contribution", "Joint_Contribution", "Total_Contribution", "Variable", "Gene", "Location", "Independendent_Pct_Explained")
		rownames(temp) <- paste0(temp$Variable, "_", temp$Gene, "_", temp$Location)
		return(temp)
	})
	print(tail(temp2))
	print(length(temp2))
	out_list <- do.call(rbind,temp2)
	return(out_list)
})

contributions_rb_regress_norm_df <- do.call(rbind, contributions_rb_regress_norm_df_list)



##### Make some figures!!! #####
pTotal_Cont <- ggplot(contributions_rb_regress_norm_df, aes(Total_Contribution, color = Variable)) +
	geom_density() +
	facet_wrap(Location ~ ., ncol = 1) +
	theme_classic()  +
	# scale_y_continuous(trans = "log10") +
	scale_x_continuous(trans = "log10") +
	scale_color_manual(values = variable_colors) +
	geom_vline(xintercept = 0.01, linetype="dashed")

save_figs(pTotal_Cont, paste0(outdir, "Total_Contribution_Histogram_rb_regressed_norm"))




pIndep_Cont <- ggplot(contributions_rb_regress_norm_df, aes(Independent_Contribution, color = Variable)) +
	geom_density() +
	facet_wrap(Location ~ ., ncol = 1) +
	theme_classic()  +
	# scale_y_continuous(trans = "log10") +
	scale_x_continuous(trans = "log10")

save_figs(pIndep_Cont, paste0(outdir, "Independent_Contribution_Histogram_rb_regressed_norm"))


pJoint_Cont <- ggplot(contributions_rb_regress_norm_df, aes(Joint_Contribution, color = Variable)) +
	geom_density() +
	facet_wrap(Location ~ ., ncol = 1) +
	theme_classic()  +
	# scale_y_continuous(trans = "log10") +
	scale_x_continuous(trans = "log10")

save_figs(pJoint_Cont, paste0(outdir, "Joint_Contribution_Histogram_rb_regressed_norm"))


pIndep_pct_Cont <- ggplot(contributions_rb_regress_norm_df, aes(Independendent_Pct_Explained, color = Variable)) +
	geom_density() +
	facet_wrap(Location ~ ., ncol = 1) +
	theme_classic()  +
	# scale_y_continuous(trans = "log10") +
	scale_x_continuous(trans = "log10")

save_figs(pIndep_pct_Cont, paste0(outdir, "Independendent_Pct_Explained_Histogram_rb_regressed_norm"))



### Variance explained correlation between sites ###
contributions_rb_regress_norm_df_joined <- full_join(contributions_rb_regress_norm_df, contributions_rb_regress_norm_df, by = c("Variable", "Gene"))
contributions_rb_regress_norm_df_joined <- contributions_rb_regress_norm_df_joined[which(contributions_rb_regress_norm_df_joined$Location.x != contributions_rb_regress_norm_df_joined$Location.y),]

### identify the combinatinons to use
pairs_df <- as.data.frame(t(combn(unique(contributions_rb_regress_norm_df_joined$Location.x), 2, simplify = TRUE)))
colnames(pairs_df) <- c("Location.x", "Location.y")

contributions_rb_regress_norm_df_joined <- left_join(pairs_df, contributions_rb_regress_norm_df_joined)


Rsquared <- unique(contributions_rb_regress_norm_df_joined[,c("Location.x","Location.y")])
Rsquared$Village_Rsquared <- NA
Rsquared$Line_Rsquared <- NA
Rsquared$Replicate_Rsquared <- NA

for (row in 1:nrow(Rsquared)){
	Rsquared$Village_Rsquared[row] <- summary(lm(Total_Contribution.y ~ Total_Contribution.x, contributions_rb_regress_norm_df_joined[which(contributions_rb_regress_norm_df_joined$Variable == "Village" & contributions_rb_regress_norm_df_joined$Location.x == Rsquared$Location.x[row] & contributions_rb_regress_norm_df_joined$Location.y == Rsquared$Location.y[row]),]))$r.squared
	Rsquared$Line_Rsquared[row] <- summary(lm(Total_Contribution.y ~ Total_Contribution.x, contributions_rb_regress_norm_df_joined[which(contributions_rb_regress_norm_df_joined$Variable == "Line" & contributions_rb_regress_norm_df_joined$Location.x == Rsquared$Location.x[row] & contributions_rb_regress_norm_df_joined$Location.y == Rsquared$Location.y[row]),]))$r.squared
	Rsquared$Replicate_Rsquared[row] <- summary(lm(Total_Contribution.y ~ Total_Contribution.x, contributions_rb_regress_norm_df_joined[which(contributions_rb_regress_norm_df_joined$Variable == "Replicate" & contributions_rb_regress_norm_df_joined$Location.x == Rsquared$Location.x[row] & contributions_rb_regress_norm_df_joined$Location.y == Rsquared$Location.y[row]),]))$r.squared
}

Rsquared_long <- pivot_longer(Rsquared, cols = c("Village_Rsquared", "Line_Rsquared", "Replicate_Rsquared"), values_to = "Rsquared", names_to = "Variable")
Rsquared_long$Variable <- gsub("_Rsquared", "", Rsquared_long$Variable)
Rsquared_long$Rsquared <- as.character(round(Rsquared_long$Rsquared, 2))



pCorr_Location_Point <- ggplot(contributions_rb_regress_norm_df_joined, aes(Total_Contribution.x, Total_Contribution.y, color = Variable)) +
	geom_point(size = 0.5, alpha = 0.25) +
	facet_grid(rows = vars(Location.y), cols = vars(Location.x)) +
	theme_classic() +
	scale_y_continuous(trans = "log10") +
	scale_x_continuous(trans = "log10") +
	scale_color_manual(values = variable_colors) +
save_figs(pCorr_Location_Point, paste0(outdir, "Correlation_Point_Location_log_rb_regressed_norm"))


pCorr_Location_Point <- ggplot(contributions_rb_regress_norm_df_joined, aes(Total_Contribution.x, Total_Contribution.y, color = Variable)) +
	geom_point(size = 0.5, alpha = 0.25) +
	facet_grid(rows = vars(Location.y), cols = vars(Location.x), switch="both") +
	theme_classic() +
	scale_color_manual(values = variable_colors) +
	xlab("Contribution to Gene Variation") +
	ylab("Contribution to Gene Variation") +
	geom_text(x = 0, y = 0.6, aes(label = paste0("R^2 = ", Rsquared)), data = Rsquared_long[which(Rsquared_long$Variable == "Village"),], hjust = 0)+
	geom_text(x = 0, y = 0.7, aes(label = paste0("R^2 = ", Rsquared)), data = Rsquared_long[which(Rsquared_long$Variable == "Line"),], hjust = 0) +
	geom_text(x = 0, y = 0.65, aes(label = paste0("R^2 = ", Rsquared)), data = Rsquared_long[which(Rsquared_long$Variable == "Replicate"),], hjust = 0)
save_figs(pCorr_Location_Point, paste0(outdir, "Correlation_Point_Location_rb_regressed_norm"), width = 15, height = 12)



lm_eqn = function(df, x, y){
    m = lm(y ~ x, df);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
         list(a = format(coef(m)[1], digits = 2), 
              b = format(coef(m)[2], digits = 2), 
             r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));                 
}


pCorr_Location_Point_Line <- ggplot(contributions_rb_regress_norm_df_joined[which(contributions_rb_regress_norm_df_joined$Variable == "Line"),], aes(Total_Contribution.x, Total_Contribution.y)) +
	geom_point(size = 0.5, alpha = 0.25, color = variable_colors["Line"]) +
	facet_grid(rows = vars(Location.y), cols = vars(Location.x)) +
	theme_classic() +
	xlab("Contribution to Gene Variation") +
	ylab("Contribution to Gene Variation") +
	geom_text(x = 0, y = 0.75, aes(label = paste0("R^2 = ", Rsquared)), data = Rsquared_long[which(Rsquared_long$Variable == "Line"),], hjust = 0, , color = variable_colors["Line"])
save_figs(pCorr_Location_Point_Line, paste0(outdir, "Correlation_Point_Location_Line_rb_regressed_norm"))


pCorr_Location_Point_Line_log <- ggplot(contributions_rb_regress_norm_df_joined[which(contributions_rb_regress_norm_df_joined$Variable == "Line"),], aes(Total_Contribution.x, Total_Contribution.y)) +
	geom_point(size = 0.5, alpha = 0.25, color = variable_colors["Line"]) +
	facet_grid(rows = vars(Location.y), cols = vars(Location.x)) +
	theme_classic() +
	scale_y_continuous(trans = "log10") +
	scale_x_continuous(trans = "log10")
save_figs(pCorr_Location_Point_Line_log, paste0(outdir, "Correlation_Point_Location_Line_log_rb_regressed_norm"))



pCorr_Location_Point_Village <- ggplot(contributions_rb_regress_norm_df_joined[which(contributions_rb_regress_norm_df_joined$Variable == "Village"),], aes(Total_Contribution.x, Total_Contribution.y)) +
	geom_point(size = 0.5, alpha = 0.25, color = variable_colors["Village"]) +
	facet_grid(rows = vars(Location.y), cols = vars(Location.x)) +
	xlab("Contribution to Gene Variation") +
	ylab("Contribution to Gene Variation") +
	theme_classic() +
	geom_text(x = 0.25, y = 0.45, aes(label = paste0("R^2 = ", Rsquared)), data = Rsquared_long[which(Rsquared_long$Variable == "Village"),], hjust = 0, , color = variable_colors["Village"])
save_figs(pCorr_Location_Point_Village, paste0(outdir, "Correlation_Point_Location_Village_rb_regressed_norm"))

pCorr_Location_Point_Village_log <- ggplot(contributions_rb_regress_norm_df_joined[which(contributions_rb_regress_norm_df_joined$Variable == "Village"),], aes(Total_Contribution.x, Total_Contribution.y)) +
	geom_point(size = 0.5, alpha = 0.25, color = variable_colors["Village"]) +
	facet_grid(rows = vars(Location.y), cols = vars(Location.x)) +
	theme_classic() +
	scale_y_continuous(trans = "log10") +
	scale_x_continuous(trans = "log10")
save_figs(pCorr_Location_Point_Village_log, paste0(outdir, "Correlation_Point_Location_Village_log_rb_regressed_norm"))


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


contributions_rb_regress_norm_df_pluri_genes <- left_join(pluri_genes, contributions_rb_regress_norm_df, by = c("ENSG_ID" = "Gene"))
contributions_rb_regress_norm_df_pluri_genes <- na.omit(contributions_rb_regress_norm_df_pluri_genes)


pPluri_Genes_Cont <- ggplot(contributions_rb_regress_norm_df_pluri_genes, aes(Location, Total_Contribution, fill = Variable)) +
						geom_bar(position = "dodge", stat = "identity") +
						theme_classic() +
						facet_wrap(Gene_ID ~ ., nrow = 1) +
						scale_fill_manual(values = variable_colors) +
						theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1)) +
						ylab("Contribution to Gene Expression Variation")
save_figs(pPluri_Genes_Cont, paste0(outdir, "Pluripotent_Gene_Variable_Contributions_rb_regressed_norm"), width = 35, height = 15)



### Make a figure for the most line-specific pluripotency gene ###
ZFP42_list <- lapply(names(seurat_rb_norm_list), function(x){
	df <- data.frame("Expression" = data.frame(seurat_rb_norm_list[[x]][["SCT"]]@scale.data[pluri_genes$ENSG_ID[which(pluri_genes$Gene_ID == "ZFP42")],]), "Time" = seurat_rb_norm_list[[x]]@meta.data$Time, "Line" = seurat_rb_norm_list[[x]]@meta.data$Final_Assignment, "Replicate" = gsub("[A-Z][a-z]+", "", seurat_rb_norm_list[[x]]@meta.data$MULTI_ID))
	colnames(df)[1] <- "Expression"
	df$Location <- x
	return(df)
})


ZFP42_df <- do.call(rbind, ZFP42_list)


ZFP42_box <- ggplot(ZFP42_df, aes(Time, Expression, color = Line)) +
				geom_boxplot(outlier.size = 0.5) +
				theme_classic() +
				facet_wrap(vars(Location), ncol = 1) +
				scale_color_manual(values = line_colors) +
				ylab("Normalized ZFP42 Expression") +
				theme(axis.title.x=element_blank())
save_figs(ZFP42_box, paste0(outdir, "ZFP42_Expression_facet"), width = 12, height = 12)


### Make a figure for the most line-specific gene: ENSG00000106 CHCHD2 ###
CHCHD2_list <- lapply(names(seurat_rb_norm_list), function(x){
	df <- data.frame("Expression" = data.frame(seurat_rb_norm_list[[x]][["SCT"]]@scale.data["ENSG00000106153",]), "Time" = seurat_rb_norm_list[[x]]@meta.data$Time, "Line" = seurat_rb_norm_list[[x]]@meta.data$Final_Assignment, "Replicate" = gsub("[A-Z][a-z]+", "", seurat_rb_norm_list[[x]]@meta.data$MULTI_ID))
	colnames(df)[1] <- "Expression"
	df$Location <- x
	return(df)
})


CHCHD2_df <- do.call(rbind, CHCHD2_list)


CHCHD2_box <- ggplot(CHCHD2_df, aes(Time, Expression, color = Line)) +
				geom_boxplot(outlier.size = 0.5) +
				theme_classic() +
				facet_wrap(vars(Location), ncol = 1) +
				scale_color_manual(values = line_colors) +
				ylab("Normalized CHCHD2 Expression") +
				theme(axis.title.x=element_blank())
save_figs(CHCHD2_box, paste0(outdir, "CHCHD2_Expression_facet"), width = 12, height = 12)



##### THIS WON'T WORK - NEED AT LEAST TWO INDEPENDENT VARIABLES FOR HEIRARCHICAL PARTITIONING #####
##### Maybe try QR decomposition since provide the same results for joint #####
##### Heirarchical Partitioning for just line effect - all conditions separated #####
##### Read in Data #####
seurat_rb <- readRDS(paste0(datadir,"seurat_integrated_all_times_clustered_rb_regressed.rds"))
seurat_rb_sub <- subset(seurat_rb, subset = Time != "Thawed Village Day 0")
seurat_rb_sub <- subset(seurat_rb_sub, subset = Time != "Thawed Village Day 7")
seurat_rb_sub$Location_Time_Rep <- paste0(seurat_rb_sub$Location_Time, "_", gsub("\\D","",seurat_rb_sub$MULTI_ID))


seurat_rb_separate_list <- lapply(unique(seurat_rb_sub@meta.data$Location_Time_Rep), function(x){
	subset(seurat_rb_sub, subset = Location_Time_Rep == x)
})
names(seurat_rb_separate_list) <- unique(seurat_rb_sub@meta.data$Location_Time_Rep)


#### Filter out Genes expressed in low percentage of cells #####
pct_expressing_rb_separate <- lapply(seurat_rb_separate_list, function(x){
	temp <- data.frame(rowSums(as.matrix(x[["SCT"]]@counts) > 0)/ncol(x[["SCT"]]@counts))
	colnames(temp) <- "Cell_proportion"
	return(temp)
})

### 10% to start
seurat_rb_separate_list <- lapply(names(seurat_rb_separate_list), function(x){
	subset(seurat_rb_separate_list[[x]], features = rownames(pct_expressing_rb_separate[[x]])[which(pct_expressing_rb_separate[[x]]$Cell_proportion >= 0.1)])
})
names(seurat_rb_separate_list) <- unique(seurat_rb_sub@meta.data$Location_Time_Rep)


sce_rb_separate_list <- lapply(seurat_rb_separate_list, function(x){
	SingleCellExperiment(assays = list(counts = x[["SCT"]]@counts, logcounts = x[["SCT"]]@data, normcounts = x[["SCT"]]@scale.data), colData = x@meta.data)
})


if (!file.exists(paste0(outdir, "heirarchical_partitioning_results_rb_regressed_separate.rds"))){
	hier_models_rb <- list()
	i <- 1
	for(x in names(seurat_rb_separate_list)){
		j <- 1
		# for (gene in unique(rownames(seurat_rb_separate_list[[x]]))){
			# print(paste(i,j))
			# df <- data.frame("Expression" = data.frame(seurat_rb_separate_list[[x]][["SCT"]]@scale.data[gene,]), "Line" = factor(seurat_rb_separate_list[[x]]@meta.data$Final_Assignment))
			# colnames(df)[1] <- "Expression"
			# hier_models_rb[[x]][[gene]] <- hier.part(df$Expression, df[,2:ncol(df)])
			df <- data.frame("Line" = factor(seurat_rb_separate_list[[x]]@meta.data$Final_Assignment))
			hier_models_rb[[x]] <- getVarianceExplained(sce_rb_separate_list[[x]], variables = df, exprs_values = "normcounts")
			j <- j + 1
		# }
		i <- 1 + i
	}

	saveRDS(hier_models_rb, paste0(outdir, "heirarchical_partitioning_results_rb_regressed_separate.rds"))
} else {
	hier_models_rb <- readRDS(paste0(outdir, "heirarchical_partitioning_results_rb_regressed_separate.rds"))
}



##### Make dataframes with all the values
contributions_rb_regress_df_list <- lapply(names(hier_models_rb), function(x){
	hier_models_rb[[x]] <- as.data.frame(hier_models_rb[[x]])
	hier_models_rb[[x]]$Gene <- rownames(hier_models_rb[[x]])
	hier_models_rb[[x]]$Location_Time_Rep <- x
	return(hier_models_rb[[x]])
})

contributions_rb_regress_df <- do.call(rbind, contributions_rb_regress_df_list)



### Variance explained correlation between sites ###
contributions_rb_regress_df_joined <- full_join(contributions_rb_regress_df, contributions_rb_regress_df, by = c("Gene"))
contributions_rb_regress_df_joined <- contributions_rb_regress_df_joined[which(contributions_rb_regress_df_joined$Location_Time_Rep.x != contributions_rb_regress_df_joined$Location_Time_Rep.y),]

Rsquared <- unique(contributions_rb_regress_df_joined[,c("Location_Time_Rep.x","Location_Time_Rep.y")])
Rsquared$Line_Rsquared <- NA

for (row in 1:nrow(Rsquared)){
	Rsquared$Line_Rsquared[row] <- summary(lm(Line.y ~ Line.x, contributions_rb_regress_df_joined[which(contributions_rb_regress_df_joined$Location_Time_Rep.x == Rsquared$Location_Time_Rep.x[row] & contributions_rb_regress_df_joined$Location_Time_Rep.y == Rsquared$Location_Time_Rep.y[row]),]))$r.squared
}

# Rsquared_long <- pivot_longer(Rsquared, cols = c("Village_Rsquared", "Line_Rsquared", "Replicate_Rsquared"), values_to = "Rsquared", names_to = "Variable")
# Rsquared_long$Variable <- gsub("_Rsquared", "", Rsquared_long$Variable)
# Rsquared_long$Rsquared <- as.character(round(Rsquared_long$Rsquared, 2))


# pCorr_Location_Point <- ggplot(contributions_rb_regress_df_joined, aes(Total_Contribution.x, Total_Contribution.y, color = Variable)) +
# 	geom_point(size = 0.5, alpha = 0.25) +
# 	facet_grid(rows = vars(Location.y), cols = vars(Location.x)) +
# 	theme_classic() +
# 	scale_y_continuous(trans = "log10") +
# 	scale_x_continuous(trans = "log10") +
# 	scale_color_manual(values = variable_colors) +
# save_figs(pCorr_Location_Point, paste0(outdir, "Correlation_Point_Location_log_rb_regressed_norm"))


pCorr_Location_Point <- ggplot(contributions_rb_regress_df_joined, aes(Line.x, Line.y)) +
	geom_point(size = 0.5, alpha = 0.25) +
	facet_grid(rows = vars(Location_Time_Rep.y), cols = vars(Location_Time_Rep.x)) +
	theme_classic() +
	scale_color_manual(values = variable_colors) +
	xlab("Contribution to Gene Variation") +
	ylab("Contribution to Gene Variation") +
	geom_text(x = 0, y = 0.6, aes(label = paste0("R^2 = ", Line_Rsquared)), data = Rsquared, hjust = 0)
save_figs(pCorr_Location_Point, paste0(outdir, "Correlation_Point_Location_rb_regressed_norm_separated"), width = 20)



lm_eqn = function(df, x, y){
    m = lm(y ~ x, df);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
         list(a = format(coef(m)[1], digits = 2), 
              b = format(coef(m)[2], digits = 2), 
             r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));                 
}


pCorr_Location_Point_Line <- ggplot(contributions_rb_regress_norm_df_joined[which(contributions_rb_regress_norm_df_joined$Variable == "Line"),], aes(Total_Contribution.x, Total_Contribution.y)) +
	geom_point(size = 0.5, alpha = 0.25, color = variable_colors["Line"]) +
	facet_grid(rows = vars(Location.y), cols = vars(Location.x)) +
	theme_classic() +
	xlab("Contribution to Gene Variation") +
	ylab("Contribution to Gene Variation") +
	geom_text(x = 0, y = 0.75, aes(label = paste0("R^2 = ", Rsquared)), data = Rsquared_long[which(Rsquared_long$Variable == "Line"),], hjust = 0, , color = variable_colors["Line"])
save_figs(pCorr_Location_Point_Line, paste0(outdir, "Correlation_Point_Location_Line_rb_regressed_norm"))




##### Try with just most variable genes #####
##### Read in Data #####
seurat_rb_sub_norm <- readRDS(paste0(outdir, "seurat_SCT_wRB_base_village_only.rds"))


seurat_rb_norm_list <- lapply(unique(seurat_rb_sub_norm@meta.data$Location), function(x){
	subset(seurat_rb_sub_norm, subset = Location == x)
})
names(seurat_rb_norm_list) <- unique(seurat_rb_sub_norm@meta.data$Location)


seurat_rb_norm_list <- lapply(seurat_rb_norm_list, function(x){
	SCTransform(x, verbose = TRUE, vars.to.regress = c("scores.G1", "scores.S", "scores.G2M", "percent.mt", "percent.rb"), return.only.var.genes = FALSE)
})



#### Filter out Genes expressed in low percentage of cells #####
pct_expressing_rb <- lapply(seurat_rb_norm_list, function(x){
	temp <- data.frame(rowSums(as.matrix(x[["SCT"]]@counts) > 0)/ncol(x[["SCT"]]@counts))
	colnames(temp) <- "Cell_proportion"
	return(temp)
})

### 10% to start
seurat_rb_norm_list <- lapply(names(seurat_rb_norm_list), function(x){
	subset(seurat_rb_norm_list[[x]], features = rownames(pct_expressing_rb[[x]])[which(pct_expressing_rb[[x]]$Cell_proportion >= 0.1)])
})
names(seurat_rb_norm_list) <- unique(seurat_rb_sub_norm@meta.data$Location)


#### Get just most variable genes #####
seurat_rb_norm_list <- lapply(seurat_rb_norm_list, function(x){
	subset(x, features = seurat_rb_sub_norm@assays$SCT@var.features)
})



if (!file.exists(paste0(outdir, "heirarchical_partitioning_results_rb_regressed_norm_var.rds"))){
	hier_models_rb_var <- list()
	i <- 1
	for(x in names(seurat_rb_norm_list)){
		j <- 1
		for (gene in unique(rownames(seurat_rb_norm_list[[x]]))){
			print(paste(i,j))
			df <- data.frame("Expression" = data.frame(seurat_rb_norm_list[[x]][["SCT"]]@scale.data[gene,]), "Village" = ifelse(seurat_rb_norm_list[[x]]@meta.data$Time == "Baseline", 0, 1), "Line" = seurat_rb_norm_list[[x]]@meta.data$Final_Assignment, "Replicate" = gsub("[A-Z][a-z]+", "", seurat_rb_norm_list[[x]]@meta.data$MULTI_ID))
			colnames(df)[1] <- "Expression"
			hier_models_rb_var[[x]][[gene]] <- hier.part(df$Expression, df[,2:ncol(df)])
			j <- j + 1
		}
		i <- 1 + i
	}
	saveRDS(hier_models_rb_var, paste0(outdir, "heirarchical_partitioning_results_rb_regressed_norm_var.rds"))
} else {
	hier_models_rb_var <- readRDS(paste0(outdir, "heirarchical_partitioning_results_rb_regressed_norm_var.rds"))
}


##### Make dataframes with all the values
contributions_rb_regress_norm_var_df_list <- lapply(names(hier_models_rb_var), function(x){
	temp2 <- lapply(names(hier_models_rb_var[[x]]), function(y){
		temp <-hier_models_rb_var[[x]][[y]]$IJ
		temp$Variable <- rownames(temp)
		temp$Gene <- y
		temp$Location <- x
		temp$Independendent_Pct_Explained <- hier_models_rb_var[[x]][[y]]$I.perc$ind.exp.var
		colnames(temp) <- c("Independent_Contribution", "Joint_Contribution", "Total_Contribution", "Variable", "Gene", "Location", "Independendent_Pct_Explained")
		rownames(temp) <- paste0(temp$Variable, "_", temp$Gene, "_", temp$Location)
		return(temp)
	})
	print(tail(temp2))
	print(length(temp2))
	out_list <- do.call(rbind,temp2)
	return(out_list)
})

contributions_rb_regress_norm_var_df <- do.call(rbind, contributions_rb_regress_norm_var_df_list)



##### Make some figures!!! #####
pTotal_Cont <- ggplot(contributions_rb_regress_norm_var_df, aes(Total_Contribution, color = Variable)) +
	geom_density() +
	facet_wrap(Location ~ ., ncol = 1) +
	theme_classic()  +
	# scale_y_continuous(trans = "log10") +
	scale_x_continuous(trans = "log10") +
	scale_color_manual(values = variable_colors) +
	geom_vline(xintercept = 0.01, linetype="dashed")

save_figs(pTotal_Cont, paste0(outdir, "Total_Contribution_Histogram_rb_regressed_norm_var"))




### Variance explained correlation between sites ###
contributions_rb_regress_norm_var_df_joined <- full_join(contributions_rb_regress_norm_var_df, contributions_rb_regress_norm_var_df, by = c("Variable", "Gene"))
contributions_rb_regress_norm_var_df_joined <- contributions_rb_regress_norm_var_df_joined[which(contributions_rb_regress_norm_var_df_joined$Location.x != contributions_rb_regress_norm_var_df_joined$Location.y),]

### identify the combinatinons to use
pairs_df <- as.data.frame(t(combn(unique(contributions_rb_regress_norm_var_df_joined$Location.x), 2, simplify = TRUE)))
colnames(pairs_df) <- c("Location.x", "Location.y")

contributions_rb_regress_norm_var_df_joined <- left_join(pairs_df, contributions_rb_regress_norm_var_df_joined)


Rsquared_var <- unique(contributions_rb_regress_norm_var_df_joined[,c("Location.x","Location.y")])
Rsquared_var$Village_Rsquared <- NA
Rsquared_var$Line_Rsquared <- NA
Rsquared_var$Replicate_Rsquared <- NA

for (row in 1:nrow(Rsquared_var)){
	Rsquared_var$Village_Rsquared[row] <- summary(lm(Total_Contribution.y ~ Total_Contribution.x, contributions_rb_regress_norm_var_df_joined[which(contributions_rb_regress_norm_var_df_joined$Variable == "Village" & contributions_rb_regress_norm_var_df_joined$Location.x == Rsquared_var$Location.x[row] & contributions_rb_regress_norm_var_df_joined$Location.y == Rsquared_var$Location.y[row]),]))$r.squared
	Rsquared_var$Line_Rsquared[row] <- summary(lm(Total_Contribution.y ~ Total_Contribution.x, contributions_rb_regress_norm_var_df_joined[which(contributions_rb_regress_norm_var_df_joined$Variable == "Line" & contributions_rb_regress_norm_var_df_joined$Location.x == Rsquared_var$Location.x[row] & contributions_rb_regress_norm_var_df_joined$Location.y == Rsquared_var$Location.y[row]),]))$r.squared
	Rsquared_var$Replicate_Rsquared[row] <- summary(lm(Total_Contribution.y ~ Total_Contribution.x, contributions_rb_regress_norm_var_df_joined[which(contributions_rb_regress_norm_var_df_joined$Variable == "Replicate" & contributions_rb_regress_norm_var_df_joined$Location.x == Rsquared_var$Location.x[row] & contributions_rb_regress_norm_var_df_joined$Location.y == Rsquared_var$Location.y[row]),]))$r.squared
}

Rsquared_var_long <- pivot_longer(Rsquared_var, cols = c("Village_Rsquared", "Line_Rsquared", "Replicate_Rsquared"), values_to = "Rsquared", names_to = "Variable")
Rsquared_var_long$Variable <- gsub("_Rsquared", "", Rsquared_var_long$Variable)
Rsquared_var_long$Rsquared <- as.character(round(Rsquared_var_long$Rsquared, 2))


pCorr_Location_Point <- ggplot(contributions_rb_regress_norm_var_df_joined, aes(Total_Contribution.x, Total_Contribution.y, color = Variable)) +
	geom_point(size = 0.5, alpha = 0.25) +
	facet_grid(rows = vars(Location.y), cols = vars(Location.x), switch="both") +
	theme_classic() +
	scale_color_manual(values = variable_colors) +
	xlab("Contribution to Gene Variation") +
	ylab("Contribution to Gene Variation") +
	geom_text(x = 0, y = 0.6, aes(label = paste0("R^2 = ", Rsquared)), data = Rsquared_var_long[which(Rsquared_var_long$Variable == "Village"),], hjust = 0)+
	geom_text(x = 0, y = 0.7, aes(label = paste0("R^2 = ", Rsquared)), data = Rsquared_var_long[which(Rsquared_var_long$Variable == "Line"),], hjust = 0) +
	geom_text(x = 0, y = 0.65, aes(label = paste0("R^2 = ", Rsquared)), data = Rsquared_var_long[which(Rsquared_var_long$Variable == "Replicate"),], hjust = 0)
save_figs(pCorr_Location_Point, paste0(outdir, "Correlation_Point_Location_rb_regressed_norm_var"), width = 15, height = 12)






##### Regress out other effects to compare line effect in each locationXtime #####
seurat <- readRDS(paste0(datadir,"seurat_integrated_all_times_clustered_rb_regressed.rds"))
seurat$Location <- gsub("_.+", "", seurat$Location_Time)
seurat_sub <- subset(seurat, subset = Time != "Thawed Village Day 0")
seurat_sub <- subset(seurat_sub, subset = Time != "Thawed Village Day 7")

seurat_list <- lapply(unique(seurat_sub@meta.data$Time), function(x){
	temp <- lapply(unique(seurat_sub@meta.data$Location), function(y){
		subset(subset(seurat_sub, subset = Time == x), subset = Location == y)
	})
	names(temp) <- unique(seurat_sub@meta.data$Location)
	return(temp)
})
names(seurat_list) <- unique(seurat_sub@meta.data$Time)


### Calculate the number of cells ###
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


if (!file.exists(paste0(outdir, "seurat_SCT_wRB_sydney_regress_all_covs.rds"))){
	seurat_norm_list <- lapply(seurat_list, function(x){
		lapply(x, function(y){
			SCTransform(y, verbose = TRUE, vars.to.regress = c("scores.G1", "scores.S", "scores.G2M", "percent.mt", "percent.rb","MULTI_ID"), return.only.var.genes = FALSE)
		})
	})
	saveRDS(seurat_norm_list, paste0(outdir, "seurat_SCT_wRB_sydney_regress_all_covs.rds"))
} else {
	seurat_norm_list <- readRDS(paste0(outdir, "seurat_SCT_wRB_sydney_regress_all_covs.rds"))
}



#### Filter out Genes expressed in low percentage of cells #####
pct_expressing <- lapply(seurat_norm_list, function(x){
	lapply(x, function(y){
		temp <- data.frame(rowSums(as.matrix(y[["SCT"]]@counts) > 0)/ncol(y[["SCT"]]@counts))
		colnames(temp) <- "Cell_proportion"
		return(temp)
	})
})


### 10% to start
seurat_norm_list <- lapply(names(seurat_norm_list), function(x){
	temp <- lapply(names(seurat_norm_list[[x]]), function(y){
		subset(seurat_norm_list[[x]][[y]], features = rownames(pct_expressing[[x]][[y]])[which(pct_expressing[[x]][[y]]$Cell_proportion >= 0.1)])
	})
	names(temp) <- names(seurat_norm_list[[x]])
	return(temp)
})
names(seurat_norm_list) <- unique(seurat_sub@meta.data$Time)


sce_separate_list <- lapply(seurat_norm_list, function(x){
	lapply(x, function(y){
		SingleCellExperiment(assays = list(counts = y[["SCT"]]@counts, logcounts = y[["SCT"]]@data, normcounts = y[["SCT"]]@scale.data), colData = y@meta.data)
	})
})



if (!file.exists(paste0(outdir, "heirarchical_partitioning_results_sep_control.rds"))){
	scran <- list()
	i <- 1
	for(x in names(seurat_norm_list)){
		for (y in names(seurat_norm_list[[x]])){
			print(i)
			df <- data.frame("Line" = seurat_norm_list[[x]][[y]]@meta.data$Final_Assignment)
			scran[[x]][[y]] <- getVarianceExplained(sce_separate_list[[x]][[y]], variables = df, exprs_values = "normcounts")
			i <- 1 + i
		}
	}
	saveRDS(scran, paste0(outdir, "heirarchical_partitioning_results_sep_control.rds"))
} else {
	scran <- readRDS(paste0(outdir, "heirarchical_partitioning_results_sep_control.rds"))
}



##### Make dataframes with all the values
contributions_cov_df_list <- lapply(names(scran), function(x){
	temp <- lapply(names(scran[[x]]), function(y){
		scran[[x]][[y]] <- as.data.frame(scran[[x]][[y]])
		scran[[x]][[y]]$Gene <- rownames(scran[[x]][[y]])
		scran[[x]][[y]]$Time <- x
		scran[[x]][[y]]$Location <- y
		return(scran[[x]][[y]])
	})
	do.call(rbind, temp)
})
contributions_cov_df <- do.call(rbind, contributions_cov_df_list)
contributions_cov_df$Location_Time <- paste0(contributions_cov_df$Location, ".", contributions_cov_df$Time)




##### Make some figures!!! #####
pTotal_Cont <- ggplot(contributions_cov_df, aes(Line)) +
	geom_density() +
	facet_grid(Time ~ Location) +
	theme_classic()  +
	# scale_y_continuous(trans = "log10") +
	scale_x_continuous(trans = "log10") +
	scale_color_manual(values = variable_colors[["Line"]]) +
	geom_vline(xintercept = 0.01, linetype="dashed")

save_figs(pTotal_Cont, paste0(outdir, "Total_Contribution_Histogram_cov"))



### Variance explained correlation between sites ###
contributions_cov_df_joined <- full_join(contributions_cov_df, contributions_cov_df, by = c("Gene"))
contributions_cov_df_joined <- contributions_cov_df_joined[which(contributions_cov_df_joined$Location_Time.x != contributions_cov_df_joined$Location_Time.y),]

Rsquared <- unique(contributions_cov_df_joined[,c("Location.x","Time.x","Location_Time.x", "Location.y", "Time.y", "Location_Time.y")])
Rsquared$Rsquared <- NA

for (row in 1:nrow(Rsquared)){
	Rsquared$Rsquared[row] <- summary(lm(Line.y ~ Line.x, contributions_cov_df_joined[which(contributions_cov_df_joined$Location_Time.x == Rsquared$Location_Time.x[row] & contributions_cov_df_joined$Location_Time.y == Rsquared$Location_Time.y[row]),]))$r.squared
}



pCorr_Location_Point_Line <- ggplot(contributions_cov_df_joined, aes(Line.x, Line.y)) +
	geom_point(size = 0.5, alpha = 0.25, color = variable_colors["Line"]) +
	facet_grid(rows = vars(Location_Time.y), cols = vars(Location_Time.x)) +
	theme_classic() +
	xlab("Contribution to Gene Variation") +
	ylab("Contribution to Gene Variation") +
	geom_text(x = 0, y = 0.75, aes(label = paste0("R^2 = ", Rsquared)), data = Rsquared, hjust = 0, color = variable_colors["Line"])
save_figs(pCorr_Location_Point_Line, paste0(outdir, "Correlation_Point_Location_Line_cov"))



### Figure with just the most variable genes ###
variable_gene_list <- lapply(seurat_norm_list, function(x){
	temp <- lapply(x, function(y){
		data.frame("Gene" = y@assays$SCT@var.features)
	})
	join_all(temp, type = "inner")
})
variable_gene <- join_all(variable_gene_list, type = "inner")


Rsquared_var <- unique(contributions_cov_df_joined[,c("Location.x","Time.x","Location_Time.x", "Location.y", "Time.y", "Location_Time.y")])
Rsquared_var$Rsquared <- NA

for (row in 1:nrow(Rsquared)){
	Rsquared_var$Rsquared[row] <- summary(lm(Line.y ~ Line.x, contributions_cov_df_joined[which(contributions_cov_df_joined$Location_Time.x == Rsquared$Location_Time.x[row] & contributions_cov_df_joined$Location_Time.y == Rsquared$Location_Time.y[row] & contributions_cov_df_joined$Gene %in% variable_gene$Gene),]))$r.squared
}




pCorr_Location_Point_Line_var <- ggplot(contributions_cov_df_joined[which(contributions_cov_df_joined$Gene %in% variable_gene$Gene),], aes(Line.x, Line.y)) +
	geom_point(size = 0.5, alpha = 0.25, color = variable_colors["Line"]) +
	facet_grid(rows = vars(Location_Time.y), cols = vars(Location_Time.x)) +
	theme_classic() +
	xlab("Contribution to Gene Variation") +
	ylab("Contribution to Gene Variation") +
	geom_text(x = 0, y = 0.75, aes(label = paste0("R^2 = ", Rsquared)), data = Rsquared_var, hjust = 0, color = variable_colors["Line"])
save_figs(pCorr_Location_Point_Line_var, paste0(outdir, "Correlation_Point_Location_Line_cov_variable"))



























# ##### Regress out village effect
# ### Start  by regressing out for each gene separately
# if (!file.exists(paste0(outdir, "heirarchical_partitioning_results_rb_regressed_norm_reg_timeXgene.rds"))){
# 	hier_models_reg_timeXgene <- list()
# 	i <- 1
# 	for(x in names(seurat_rb_norm_list)){
# 		j <- 1
# 		for (gene in unique(rownames(seurat_rb_norm_list[[x]]))){
# 			print(paste(i,j))
# 			df <- data.frame("Expression" = data.frame(seurat_rb_norm_list[[x]][["SCT"]]@scale.data[gene,]), "Village" = ifelse(seurat_rb_norm_list[[x]]@meta.data$Time == "Baseline", 0, 1), "Line" = seurat_rb_norm_list[[x]]@meta.data$Final_Assignment, "Replicate" = gsub("[A-Z][a-z]+", "", seurat_rb_norm_list[[x]]@meta.data$MULTI_ID))
# 			colnames(df)[1] <- "Expression"
# 			df$Residuals <- resid(lm(Expression ~ Village, df))
# 			hier_models_reg_timeXgene[[x]][[gene]] <- hier.part(df$Residuals, df[,c("Line","Replicate")])
# 			j <- j + 1
# 		}
# 		i <- 1 + i
# 	}

# 	saveRDS(hier_models_reg_timeXgene, paste0(outdir, "heirarchical_partitioning_results_rb_regressed_norm_reg_timeXgene.rds"))
# } else {
# 	hier_models_reg_timeXgene <- readRDS(paste0(outdir, "heirarchical_partitioning_results_rb_regressed_norm_reg_timeXgene.rds"))
# }




# ##### Make dataframes with all the values
# contributions_reg_timeXgene_df_list <- lapply(names(hier_models_reg_timeXgene), function(x){
# 	temp2 <- lapply(names(hier_models_reg_timeXgene[[x]]), function(y){
# 		temp <- hier_models_reg_timeXgene[[x]][[y]]$IJ
# 		temp$Variable <- rownames(temp)
# 		temp$Gene <- y
# 		temp$Location <- x
# 		temp$Independendent_Pct_Explained <- hier_models_reg_timeXgene[[x]][[y]]$I.perc$ind.exp.var
# 		colnames(temp) <- c("Independent_Contribution", "Joint_Contribution", "Total_Contribution", "Variable", "Gene", "Location", "Independendent_Pct_Explained")
# 		rownames(temp) <- paste0(temp$Variable, "_", temp$Gene, "_", temp$Location)
# 		return(temp)
# 	})
# 	print(tail(temp2))
# 	print(length(temp2))
# 	out_list <- do.call(rbind,temp2)
# 	return(out_list)
# })

# contributions_reg_timeXgene_df <- do.call(rbind, contributions_reg_timeXgene_df_list)


# ### Join this dataframe with the one that didn't regress time
# joined_contributions_timeXgene <- full_join(contributions_rb_regress_df, contributions_reg_timeXgene_df, by = c("Variable", "Gene", "Location"))


# ### Plot the relationship of the original to the post village-regressed
# pCorr_Point_Independent <- ggplot(joined_contributions_timeXgene[which(joined_contributions_timeXgene$Variable == "Line"),], aes(Independent_Contribution.x, Independent_Contribution.y)) +
# 	geom_point(size = 0.5, alpha = 0.25) +
# 	facet_wrap(Location ~., nrow = 1) +
# 	theme_classic() +
# 	# scale_y_continuous(trans = "log10") +
# 	# scale_x_continuous(trans = "log10") 
# 	# geom_smooth(method="lm", se = TRUE) +
# 	# stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
#     #            label.x.npc = "right", label.y.npc = 0.15,
#     #            formula = formula, parse = TRUE, size = 3)
# 	# annotate(x=1e-05, y =1 ,  geom='text',
#     #           label = lm_eqn(joined_contributions_timeXgene[which(joined_contributions_timeXgene$Variable == "Line"),]), size=5,parse=TRUE)
# save_figs(pCorr_Point_Independent, paste0(outdir, "Correlation_Point_Pre_Post_timeXgene_Regression"), height = 6)


# pCorr_Point_Total <- ggplot(joined_contributions_timeXgene[which(joined_contributions_timeXgene$Variable == "Line"),], aes(Total_Contribution.x, Total_Contribution.y)) +
# 	geom_point(size = 0.5, alpha = 0.25) +
# 	facet_wrap(Location ~., nrow = 1) +
# 	theme_classic() +
# 	scale_y_continuous(trans = "log10") +
# 	scale_x_continuous(trans = "log10") 
# 	# geom_smooth(method="lm", se = TRUE) +
# 	# stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
#     #            label.x.npc = "right", label.y.npc = 0.15,
#     #            formula = formula, parse = TRUE, size = 3)
# 	# annotate(x=1e-05, y =1 ,  geom='text',
#     #           label = lm_eqn(joined_contributions_timeXgene[which(joined_contributions_timeXgene$Variable == "Line"),]), size=5,parse=TRUE)
# save_figs(pCorr_Point_Total, paste0(outdir, "Correlation_Point_Total_Pre_Post_timeXgene_Regression"), height = 6)



# contributions_reg_timeXgene_df_joined <- full_join(contributions_reg_timeXgene_df, contributions_reg_timeXgene_df, by = c("Variable", "Gene"))
# pCorr_Location_Point_Line <- ggplot(contributions_reg_timeXgene_df_joined[which(contributions_reg_timeXgene_df_joined$Variable == "Line"),], aes(Total_Contribution.x, Independent_Contribution.y)) +
# 	geom_point(size = 0.5, alpha = 0.25) +
# 	facet_grid(rows = vars(Location.y), cols = vars(Location.x)) +
# 	theme_classic() 
# 	# scale_y_continuous(trans = "log10") +
# 	# scale_x_continuous(trans = "log10") +
# 	# geom_smooth(method="lm", se = TRUE) +
# 	# stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
#     #            label.x.npc = "right", label.y.npc = 0.15,
#     #            formula = formula, parse = TRUE, size = 3)
# 	# annotate(x=1e-05, y =1 ,  geom='text',
#     #           label = lm_eqn(contributions_reg_timeXgene_df_joined[which(contributions_reg_timeXgene_df_joined$Variable == "Line"),]), size=5,parse=TRUE)
# save_figs(pCorr_Location_Point_Line, paste0(outdir, "Correlation_Point_Location_Line_regression_timeXgene"))




# ##### Regress out effect of the villages and compare the effect of the individual with and without #####
# seurat_rb_time_sub_norm <- SCTransform(seurat_rb_sub, verbose = TRUE, vars.to.regress = c("scores.G1", "scores.S", "scores.G2M", "percent.mt", "percent.rb", "Time"), return.only.var.genes = FALSE)
# saveRDS(seurat_rb_time_sub_norm, paste0(outdir, "seurat_time_regressed.rds"))



# seurat_rb_time_norm_list <- lapply(unique(seurat_rb_time_sub_norm@meta.data$Location), function(x){
# 	subset(seurat_rb_time_sub_norm, subset = Location == x)
# })
# names(seurat_rb_time_norm_list) <- unique(seurat_rb_time_sub_norm@meta.data$Location)


# #### Filter out Genes expressed in low percentage of cells #####
# pct_expressing_rb <- lapply(seurat_rb_time_norm_list, function(x){
# 	temp <- data.frame(rowSums(as.matrix(x[["SCT"]]@counts) > 0)/ncol(x[["SCT"]]@counts))
# 	colnames(temp) <- "Cell_proportion"
# 	return(temp)
# })

# ### 10% to start
# seurat_rb_time_norm_list <- lapply(names(seurat_rb_time_norm_list), function(x){
# 	subset(seurat_rb_time_norm_list[[x]], features = rownames(pct_expressing_rb[[x]])[which(pct_expressing_rb[[x]]$Cell_proportion >= 0.1)])
# })
# names(seurat_rb_time_norm_list) <- unique(seurat_rb_time_sub_norm@meta.data$Location)


# if (!file.exists(paste0(outdir, "heirarchical_partitioning_results_rb_time_regressed_norm.rds"))){
# 	hier_models_rb_time <- list()
# 	i <- 1
# 	for(x in names(seurat_rb_time_norm_list)){
# 		j <- 1
# 		for (gene in unique(rownames(seurat_rb_time_norm_list[[x]]))){
# 			print(paste(i,j))
# 			df <- data.frame("Expression" = data.frame(seurat_rb_time_norm_list[[x]][["SCT"]]@scale.data[gene,]), "Line" = seurat_rb_time_norm_list[[x]]@meta.data$Final_Assignment, "Replicate" = gsub("[A-Z][a-z]+", "", seurat_rb_time_norm_list[[x]]@meta.data$MULTI_ID))
# 			colnames(df)[1] <- "Expression"
# 			hier_models_rb_time[[x]][[gene]] <- hier.part(df$Expression, df[,2:ncol(df)])
# 			j <- j + 1
# 		}
# 		i <- 1 + i
# 	}

# 	saveRDS(hier_models_rb_time, paste0(outdir, "heirarchical_partitioning_results_rb_time_regressed_norm.rds"))
# } else {
# 	hier_models_rb_time <- readRDS(paste0(outdir, "heirarchical_partitioning_results_rb_time_regressed_norm.rds"))
# }




# ##### Make dataframes with all the values
# contributions_rb_time_regress_df_list <- lapply(names(hier_models_rb_time), function(x){
# 	temp2 <- lapply(names(hier_models_rb_time[[x]]), function(y){
# 		temp <-hier_models_rb_time[[x]][[y]]$IJ
# 		temp$Variable <- rownames(temp)
# 		temp$Gene <- y
# 		temp$Location <- x
# 		temp$Independendent_Pct_Explained <- hier_models_rb_time[[x]][[y]]$I.perc$ind.exp.var
# 		colnames(temp) <- c("Independent_Contribution", "Joint_Contribution", "Total_Contribution", "Variable", "Gene", "Location", "Independendent_Pct_Explained")
# 		rownames(temp) <- paste0(temp$Variable, "_", temp$Gene, "_", temp$Location)
# 		return(temp)
# 	})
# 	print(tail(temp2))
# 	print(length(temp2))
# 	out_list <- do.call(rbind,temp2)
# 	return(out_list)
# })

# contributions_rb_time_regress_df <- do.call(rbind, contributions_rb_time_regress_df_list)


# ### Join this dataframe with the one that didn't regress time
# joined_contributions <- full_join(contributions_rb_regress_df, contributions_rb_time_regress_df, by = c("Variable", "Gene", "Location"))


# ### Plot the relationship of the original to the post village-regressed
# pCorr_Point_Independent <- ggplot(joined_contributions[which(joined_contributions$Variable == "Line"),], aes(Independent_Contribution.x, Independent_Contribution.y)) +
# 	geom_point(size = 0.5, alpha = 0.25) +
# 	facet_wrap(Location ~., nrow = 1) +
# 	theme_classic() +
# 	scale_y_continuous(trans = "log10") +
# 	scale_x_continuous(trans = "log10") 
# 	# geom_smooth(method="lm", se = TRUE) +
# 	# stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
#     #            label.x.npc = "right", label.y.npc = 0.15,
#     #            formula = formula, parse = TRUE, size = 3)
# 	# annotate(x=1e-05, y =1 ,  geom='text',
#     #           label = lm_eqn(joined_contributions[which(joined_contributions$Variable == "Line"),]), size=5,parse=TRUE)
# save_figs(pCorr_Point_Independent, paste0(outdir, "Correlation_Point_Pre_Post_Village_Regression"), height = 6)


# pCorr_Point_Total <- ggplot(joined_contributions[which(joined_contributions$Variable == "Line"),], aes(Total_Contribution.x, Total_Contribution.y)) +
# 	geom_point(size = 0.5, alpha = 0.25) +
# 	facet_wrap(Location ~., nrow = 1) +
# 	theme_classic() +
# 	scale_y_continuous(trans = "log10") +
# 	scale_x_continuous(trans = "log10") 
# 	# geom_smooth(method="lm", se = TRUE) +
# 	# stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
#     #            label.x.npc = "right", label.y.npc = 0.15,
#     #            formula = formula, parse = TRUE, size = 3)
# 	# annotate(x=1e-05, y =1 ,  geom='text',
#     #           label = lm_eqn(joined_contributions[which(joined_contributions$Variable == "Line"),]), size=5,parse=TRUE)
# save_figs(pCorr_Point_Total, paste0(outdir, "Correlation_Point_Total_Pre_Post_Village_Regression"), height = 6)
