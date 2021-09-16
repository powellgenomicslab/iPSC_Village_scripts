library(ggplot2)
library(tidyverse)
library(Seurat)
library(lme4)
library(hier.part)
library(scater)
library(data.table)


save_figs <- function(plot, basename, width = 17, height = 17, units = "cm"){
    ggsave(plot, filename = paste0(basename,".png"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".pdf"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".eps"), height = height, width = width, units = units)
}


##### Set Up Directories #####
dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
outdir <- paste0(dir,"output/Variance_contributions_lm/")
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
	# models <- list()
	hier_models <- list()
	# scater_result <- list()
	i <- 1
	for(x in names(seurat_list)){
		j <- 1
		for (gene in unique(rownames(seurat_list[[x]]))){
			print(paste(i,j))
			df <- data.frame("Expression" = t(data.frame(seurat_list[[x]][["SCT"]][gene,])), "Village" = ifelse(seurat_list[[x]]@meta.data$Time == "Baseline", 0, 1), "Line" = seurat_list[[x]]@meta.data$Final_Assignment, "Replicate" = gsub("[A-Z][a-z]+", "", seurat_list[[x]]@meta.data$MULTI_ID))
			colnames(df)[1] <- "Expression"
			# models[[x]][[gene]] <- summary(lmer(Expression ~ (1 | Village) + (1 | Line) + (1 | Replicate), data = df))
			hier_models[[x]][[gene]] <- hier.part(df$Expression, df[,2:ncol(df)])
			j <- j + 1
		}
		# scater_result[[x]] <- getVarianceExplained(seuat_list[[x]][["SCT"]][VariableFeatures(seuat_list[[x]])[1:10],], variables = df[,2:ncol(df)])
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
pCorr_Location_Point <- ggplot(contributions_df_joined, aes(Independendent_Pct_Explained.x, Independendent_Pct_Explained.y, color = Variable)) +
	geom_point(size = 0.5, alpha = 0.25) +
	facet_grid(rows = vars(Location.y), cols = vars(Location.x)) +
	theme_classic() +
	scale_y_continuous(trans = "log10") +
	scale_x_continuous(trans = "log10")
save_figs(pCorr_Location_Point, paste0(outdir, "Correlation_Point_Location"))


pCorr_Location_Point_Line <- ggplot(contributions_df_joined[which(contributions_df_joined$Variable == "Line"),], aes(Independendent_Pct_Explained.x, Independendent_Pct_Explained.y)) +
	geom_point(size = 0.5, alpha = 0.25) +
	facet_grid(rows = vars(Location.y), cols = vars(Location.x)) +
	theme_classic() +
	scale_y_continuous(trans = "log10") +
	scale_x_continuous(trans = "log10")
save_figs(pCorr_Location_Point_Line, paste0(outdir, "Correlation_Point_Location_Line"))


pCorr_Location_Point_Village <- ggplot(contributions_df_joined[which(contributions_df_joined$Variable == "Village"),], aes(Independendent_Pct_Explained.x, Independendent_Pct_Explained.y)) +
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
	# models <- list()
	hier_models_rb <- list()
	# scater_result <- list()
	i <- 1
	for(x in names(seurat_rb_list)){
		j <- 1
		for (gene in unique(rownames(seurat_rb_list[[x]]))){
			print(paste(i,j))
			df <- data.frame("Expression" = t(data.frame(seurat_rb_list[[x]][["SCT"]][gene,])), "Village" = ifelse(seurat_rb_list[[x]]@meta.data$Time == "Baseline", 0, 1), "Line" = seurat_rb_list[[x]]@meta.data$Final_Assignment, "Replicate" = gsub("[A-Z][a-z]+", "", seurat_rb_list[[x]]@meta.data$MULTI_ID))
			colnames(df)[1] <- "Expression"
			# models[[x]][[gene]] <- summary(lmer(Expression ~ (1 | Village) + (1 | Line) + (1 | Replicate), data = df))
			hier_models_rb[[x]][[gene]] <- hier.part(df$Expression, df[,2:ncol(df)])
			j <- j + 1
		}
		# scater_result[[x]] <- getVarianceExplained(seuat_list[[x]][["SCT"]][VariableFeatures(seuat_list[[x]])[1:10],], variables = df[,2:ncol(df)])
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

