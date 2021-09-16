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
outdir <- paste0(dir,"output/Variance_contributions_lm_scale_data_Sydney/")
datadir <- paste0(dir,"output/All_data_integrated_remove_bad/")

dir.create(outdir)




##### Redo with rb % regressed out #####
##### Read in Data #####
seurat_rb <- readRDS(paste0(datadir,"seurat_integrated_all_times_clustered_rb_regressed.rds"))
seurat_rb$Location <- gsub("_.+", "", seurat_rb$Location_Time)
seurat_rb_sub <- subset(seurat_rb, subset = Location != "Brisbane")
seurat_rb_sub <- subset(seurat_rb_sub, subset = Location != "Melbourne")

seurat_rb_list <- lapply(unique(seurat_rb_sub@meta.data$Time), function(x){
	subset(seurat_rb_sub, subset = Time == x)
})
names(seurat_rb_list) <- unique(seurat_rb_sub@meta.data$Time)





if (!file.exists(paste0(outdir, "seurat_SCT_wRB_sydney.rds"))){
	seurat_rb_sub_norm <- SCTransform(seurat_rb_sub, verbose = TRUE, vars.to.regress = c("scores.G1", "scores.S", "scores.G2M", "percent.mt", "percent.rb"), return.only.var.genes = FALSE)
	saveRDS(seurat_rb_sub_norm, paste0(outdir, "seurat_SCT_wRB_sydney.rds"))
} else {
	seurat_rb_sub_norm <- readRDS(paste0(outdir, "seurat_SCT_wRB_sydney.rds"))
}



seurat_rb_norm_list <- lapply(unique(seurat_rb_sub_norm@meta.data$Time), function(x){
	subset(seurat_rb_sub_norm, subset = Time == x)
})
names(seurat_rb_norm_list) <- unique(seurat_rb_sub_norm@meta.data$Time)


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
names(seurat_rb_norm_list) <- unique(seurat_rb_sub_norm@meta.data$Time)


if (!file.exists(paste0(outdir, "heirarchical_partitioning_results_rb_regressed_norm_syd.rds"))){
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
	saveRDS(hier_models_rb, paste0(outdir, "heirarchical_partitioning_results_rb_regressed_norm_syd.rds"))
} else {
	hier_models_rb <- readRDS(paste0(outdir, "heirarchical_partitioning_results_rb_regressed_norm_syd.rds"))
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
contributions_rb_regress_norm_df <- contributions_rb_regress_norm_df[which(contributions_rb_regress_norm_df$Variable != "Village"),]


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

Rsquared <- unique(contributions_rb_regress_norm_df_joined[,c("Location.x","Location.y")])
Rsquared$Line_Rsquared <- NA
Rsquared$Replicate_Rsquared <- NA

for (row in 1:nrow(Rsquared)){
	Rsquared$Line_Rsquared[row] <- summary(lm(Total_Contribution.y ~ Total_Contribution.x, contributions_rb_regress_norm_df_joined[which(contributions_rb_regress_norm_df_joined$Variable == "Line" & contributions_rb_regress_norm_df_joined$Location.x == Rsquared$Location.x[row] & contributions_rb_regress_norm_df_joined$Location.y == Rsquared$Location.y[row]),]))$r.squared
	Rsquared$Replicate_Rsquared[row] <- summary(lm(Total_Contribution.y ~ Total_Contribution.x, contributions_rb_regress_norm_df_joined[which(contributions_rb_regress_norm_df_joined$Variable == "Replicate" & contributions_rb_regress_norm_df_joined$Location.x == Rsquared$Location.x[row] & contributions_rb_regress_norm_df_joined$Location.y == Rsquared$Location.y[row]),]))$r.squared
}

Rsquared_long <- pivot_longer(Rsquared, cols = c( "Line_Rsquared", "Replicate_Rsquared"), values_to = "Rsquared", names_to = "Variable")
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
	facet_grid(rows = vars(Location.y), cols = vars(Location.x)) +
	theme_classic() +
	scale_color_manual(values = variable_colors) +
	xlab("Contribution to Gene Variation") +
	ylab("Contribution to Gene Variation") +
	geom_text(x = 0, y = 0.75, aes(label = paste0("R^2 = ", Rsquared)), data = Rsquared_long[which(Rsquared_long$Variable == "Line"),], hjust = 0) +
	geom_text(x = 0, y = 0.675, aes(label = paste0("R^2 = ", Rsquared)), data = Rsquared_long[which(Rsquared_long$Variable == "Replicate"),], hjust = 0)
save_figs(pCorr_Location_Point, paste0(outdir, "Correlation_Point_Location_rb_regressed_norm"), width = 20)


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



# pCorr_Location_Point_Village <- ggplot(contributions_rb_regress_norm_df_joined[which(contributions_rb_regress_norm_df_joined$Variable == "Village"),], aes(Total_Contribution.x, Total_Contribution.y)) +
# 	geom_point(size = 0.5, alpha = 0.25, color = variable_colors["Village"]) +
# 	facet_grid(rows = vars(Location.y), cols = vars(Location.x)) +
# 	xlab("Contribution to Gene Variation") +
# 	ylab("Contribution to Gene Variation") +
# 	theme_classic() +
# 	geom_text(x = 0.25, y = 0.45, aes(label = paste0("R^2 = ", Rsquared)), data = Rsquared_long[which(Rsquared_long$Variable == "Village"),], hjust = 0, , color = variable_colors["Village"])
# save_figs(pCorr_Location_Point_Village, paste0(outdir, "Correlation_Point_Location_Village_rb_regressed_norm"))

# pCorr_Location_Point_Village_log <- ggplot(contributions_rb_regress_norm_df_joined[which(contributions_rb_regress_norm_df_joined$Variable == "Village"),], aes(Total_Contribution.x, Total_Contribution.y)) +
# 	geom_point(size = 0.5, alpha = 0.25, color = variable_colors["Village"]) +
# 	facet_grid(rows = vars(Location.y), cols = vars(Location.x)) +
# 	theme_classic() +
# 	scale_y_continuous(trans = "log10") +
# 	scale_x_continuous(trans = "log10")
# save_figs(pCorr_Location_Point_Village_log, paste0(outdir, "Correlation_Point_Location_Village_log_rb_regressed_norm"))


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









##### Try normalizing each timepoint separately #####
if (!file.exists(paste0(outdir, "seurat_SCT_wRB_sydney_separate.rds"))){
	seurat_rb_sub_norm_list <- lapply(seurat_rb_list, function(x){
		SCTransform(x, verbose = TRUE, vars.to.regress = c("scores.G1", "scores.S", "scores.G2M", "percent.mt", "percent.rb"), return.only.var.genes = FALSE)
	})
	saveRDS(seurat_rb_sub_norm_list, paste0(outdir, "seurat_SCT_wRB_sydney_separate.rds"))
} else {
	seurat_rb_sub_norm_list <- readRDS(paste0(outdir, "seurat_SCT_wRB_sydney_separate.rds"))
}



# seurat_rb_norm_list <- lapply(unique(seurat_rb_sub_norm@meta.data$Time), function(x){
# 	subset(seurat_rb_sub_norm, subset = Time == x)
# })
# names(seurat_rb_norm_list) <- unique(seurat_rb_sub_norm@meta.data$Time)


#### Filter out Genes expressed in low percentage of cells #####
pct_expressing_rb <- lapply(seurat_rb_sub_norm_list, function(x){
	temp <- data.frame(rowSums(as.matrix(x[["SCT"]]@counts) > 0)/ncol(x[["SCT"]]@counts))
	colnames(temp) <- "Cell_proportion"
	return(temp)
})

### 10% to start
seurat_rb_sub_norm_list <- lapply(names(seurat_rb_sub_norm_list), function(x){
	subset(seurat_rb_sub_norm_list[[x]], features = rownames(pct_expressing_rb[[x]])[which(pct_expressing_rb[[x]]$Cell_proportion >= 0.1)])
})
names(seurat_rb_sub_norm_list) <- unique(seurat_rb_sub_norm@meta.data$Time)


if (!file.exists(paste0(outdir, "heirarchical_partitioning_results_rb_regressed_norm_sep_syd.rds"))){
	hier_models_rb <- list()
	i <- 1
	for(x in names(seurat_rb_sub_norm_list)){
		j <- 1
		for (gene in unique(rownames(seurat_rb_sub_norm_list[[x]]))){
			print(paste(i,j))
			df <- data.frame("Expression" = data.frame(seurat_rb_sub_norm_list[[x]][["SCT"]]@scale.data[gene,]), "Line" = seurat_rb_sub_norm_list[[x]]@meta.data$Final_Assignment, "Replicate" = gsub("[A-Z][a-z]+", "", seurat_rb_sub_norm_list[[x]]@meta.data$MULTI_ID))
			colnames(df)[1] <- "Expression"
			hier_models_rb[[x]][[gene]] <- hier.part(df$Expression, df[,2:ncol(df)])
			j <- j + 1
		}
		i <- 1 + i
	}
	saveRDS(hier_models_rb, paste0(outdir, "heirarchical_partitioning_results_rb_regressed_norm_sep_syd.rds"))
} else {
	hier_models_rb <- readRDS(paste0(outdir, "heirarchical_partitioning_results_rb_regressed_norm_sep_syd.rds"))
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

save_figs(pTotal_Cont, paste0(outdir, "Total_Contribution_Histogram_rb_regressed_norm_sep"))




pIndep_Cont <- ggplot(contributions_rb_regress_norm_df, aes(Independent_Contribution, color = Variable)) +
	geom_density() +
	facet_wrap(Location ~ ., ncol = 1) +
	theme_classic()  +
	# scale_y_continuous(trans = "log10") +
	scale_x_continuous(trans = "log10")

save_figs(pIndep_Cont, paste0(outdir, "Independent_Contribution_Histogram_rb_regressed_norm_sep"))


pJoint_Cont <- ggplot(contributions_rb_regress_norm_df, aes(Joint_Contribution, color = Variable)) +
	geom_density() +
	facet_wrap(Location ~ ., ncol = 1) +
	theme_classic()  +
	# scale_y_continuous(trans = "log10") +
	scale_x_continuous(trans = "log10")

save_figs(pJoint_Cont, paste0(outdir, "Joint_Contribution_Histogram_rb_regressed_norm_sep"))


pIndep_pct_Cont <- ggplot(contributions_rb_regress_norm_df, aes(Independendent_Pct_Explained, color = Variable)) +
	geom_density() +
	facet_wrap(Location ~ ., ncol = 1) +
	theme_classic()  +
	# scale_y_continuous(trans = "log10") +
	scale_x_continuous(trans = "log10")

save_figs(pIndep_pct_Cont, paste0(outdir, "Independendent_Pct_Explained_Histogram_rb_regressed_norm_sep"))



### Variance explained correlation between sites ###
contributions_rb_regress_norm_df_joined <- full_join(contributions_rb_regress_norm_df, contributions_rb_regress_norm_df, by = c("Variable", "Gene"))
contributions_rb_regress_norm_df_joined <- contributions_rb_regress_norm_df_joined[which(contributions_rb_regress_norm_df_joined$Location.x != contributions_rb_regress_norm_df_joined$Location.y),]

Rsquared <- unique(contributions_rb_regress_norm_df_joined[,c("Location.x","Location.y")])
Rsquared$Line_Rsquared <- NA
Rsquared$Replicate_Rsquared <- NA

for (row in 1:nrow(Rsquared)){
	Rsquared$Line_Rsquared[row] <- summary(lm(Total_Contribution.y ~ Total_Contribution.x, contributions_rb_regress_norm_df_joined[which(contributions_rb_regress_norm_df_joined$Variable == "Line" & contributions_rb_regress_norm_df_joined$Location.x == Rsquared$Location.x[row] & contributions_rb_regress_norm_df_joined$Location.y == Rsquared$Location.y[row]),]))$r.squared
	Rsquared$Replicate_Rsquared[row] <- summary(lm(Total_Contribution.y ~ Total_Contribution.x, contributions_rb_regress_norm_df_joined[which(contributions_rb_regress_norm_df_joined$Variable == "Replicate" & contributions_rb_regress_norm_df_joined$Location.x == Rsquared$Location.x[row] & contributions_rb_regress_norm_df_joined$Location.y == Rsquared$Location.y[row]),]))$r.squared
}

Rsquared_long <- pivot_longer(Rsquared, cols = c( "Line_Rsquared", "Replicate_Rsquared"), values_to = "Rsquared", names_to = "Variable")
Rsquared_long$Variable <- gsub("_Rsquared", "", Rsquared_long$Variable)
Rsquared_long$Rsquared <- as.character(round(Rsquared_long$Rsquared, 2))


pCorr_Location_Point <- ggplot(contributions_rb_regress_norm_df_joined, aes(Total_Contribution.x, Total_Contribution.y, color = Variable)) +
	geom_point(size = 0.5, alpha = 0.25) +
	facet_grid(rows = vars(Location.y), cols = vars(Location.x)) +
	theme_classic() +
	scale_y_continuous(trans = "log10") +
	scale_x_continuous(trans = "log10") +
	scale_color_manual(values = variable_colors) +
save_figs(pCorr_Location_Point, paste0(outdir, "Correlation_Point_Location_log_rb_regressed_norm_sep"))


pCorr_Location_Point <- ggplot(contributions_rb_regress_norm_df_joined, aes(Total_Contribution.x, Total_Contribution.y, color = Variable)) +
	geom_point(size = 0.5, alpha = 0.25) +
	facet_grid(rows = vars(Location.y), cols = vars(Location.x)) +
	theme_classic() +
	scale_color_manual(values = variable_colors) +
	xlab("Contribution to Gene Variation") +
	ylab("Contribution to Gene Variation") +
	geom_text(x = 0, y = 0.75, aes(label = paste0("R^2 = ", Rsquared)), data = Rsquared_long[which(Rsquared_long$Variable == "Line"),], hjust = 0) +
	geom_text(x = 0, y = 0.675, aes(label = paste0("R^2 = ", Rsquared)), data = Rsquared_long[which(Rsquared_long$Variable == "Replicate"),], hjust = 0)
save_figs(pCorr_Location_Point, paste0(outdir, "Correlation_Point_Location_rb_regressed_norm_sep"), width = 20)


pCorr_Location_Point_Line <- ggplot(contributions_rb_regress_norm_df_joined[which(contributions_rb_regress_norm_df_joined$Variable == "Line"),], aes(Total_Contribution.x, Total_Contribution.y)) +
	geom_point(size = 0.5, alpha = 0.25, color = variable_colors["Line"]) +
	facet_grid(rows = vars(Location.y), cols = vars(Location.x)) +
	theme_classic() +
	xlab("Contribution to Gene Variation") +
	ylab("Contribution to Gene Variation") +
	geom_text(x = 0, y = 0.75, aes(label = paste0("R^2 = ", Rsquared)), data = Rsquared_long[which(Rsquared_long$Variable == "Line"),], hjust = 0, , color = variable_colors["Line"])
save_figs(pCorr_Location_Point_Line, paste0(outdir, "Correlation_Point_Location_Line_rb_regressed_norm_sep"))


pCorr_Location_Point_Line_log <- ggplot(contributions_rb_regress_norm_df_joined[which(contributions_rb_regress_norm_df_joined$Variable == "Line"),], aes(Total_Contribution.x, Total_Contribution.y)) +
	geom_point(size = 0.5, alpha = 0.25, color = variable_colors["Line"]) +
	facet_grid(rows = vars(Location.y), cols = vars(Location.x)) +
	theme_classic() +
	scale_y_continuous(trans = "log10") +
	scale_x_continuous(trans = "log10")
save_figs(pCorr_Location_Point_Line_log, paste0(outdir, "Correlation_Point_Location_Line_log_rb_regressed_norm_sep"))



##### Try including village effect as well #####
seurat_rb_sub$Cryopreservation <- ifelse((seurat_rb_sub$Time == "Thawed Village Day 0" | seurat_rb_sub$Time == "Thawed Village Day 7"), "Post-cryopreservation", "Fresh")

seurat_rb_list <- lapply(unique(seurat_rb_sub@meta.data$Cryopreservation), function(x){
	subset(seurat_rb_sub, subset = Cryopreservation == x)
})
names(seurat_rb_list) <- unique(seurat_rb_sub@meta.data$Cryopreservation)


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
names(seurat_rb_list) <- unique(seurat_rb_sub$Cryopreservation)


if (!file.exists(paste0(outdir, "heirarchical_partitioning_results_rb_regressed_syd_combined.rds"))){
	hier_models <- list()
	i <- 1
	for(x in names(seurat_rb_list)){
		j <- 1
		for (gene in unique(rownames(seurat_rb_list[[x]]))){
			print(paste(i,j))
			df <- data.frame("Expression" = data.frame(seurat_rb_list[[x]][["SCT"]]@scale.data[gene,]), "Line" = seurat_rb_list[[x]]@meta.data$Final_Assignment, "Replicate" = gsub("[A-Z][a-z]+", "", seurat_rb_list[[x]]@meta.data$MULTI_ID), "Village" = ifelse(seurat_rb_list[[x]]@meta.data$Time == "Baseline", 0, ifelse(seurat_rb_list[[x]]@meta.data$Time == "Thawed Village Day 0", 0, 1)))
			colnames(df)[1] <- "Expression"
			hier_models[[x]][[gene]] <- hier.part(df$Expression, df[,2:ncol(df)])
			j <- j + 1
		}
		i <- 1 + i
	}
	saveRDS(hier_models, paste0(outdir, "heirarchical_partitioning_results_rb_regressed_syd_combined.rds"))
} else {
	hier_models <- readRDS(paste0(outdir, "heirarchical_partitioning_results_rb_regressed_syd_combined.rds"))
}



##### Make dataframes with all the values
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
	theme_classic()  +
	facet_wrap(Location ~ ., ncol = 1) +
	# scale_y_continuous(trans = "log10") +
	scale_x_continuous(trans = "log10") +
	scale_color_manual(values = variable_colors) +
	geom_vline(xintercept = 0.01, linetype="dashed")

save_figs(pTotal_Cont, paste0(outdir, "Total_Contribution_Histogram_sydney_combined"))




### Variance explained correlation between sites ###
contributions_df_joined <- full_join(contributions_df, contributions_df, by = c("Variable", "Gene"))
contributions_df_joined <- contributions_df_joined[which(contributions_df_joined$Location.x != contributions_df_joined$Location.y),]

Rsquared <- unique(contributions_df_joined[,c("Location.x","Location.y")])
Rsquared$Village_Rsquared <- NA
Rsquared$Line_Rsquared <- NA
Rsquared$Replicate_Rsquared <- NA

for (row in 1:nrow(Rsquared)){
	Rsquared$Village_Rsquared[row] <- summary(lm(Total_Contribution.y ~ Total_Contribution.x, contributions_df_joined[which(contributions_df_joined$Variable == "Village" & contributions_df_joined$Location.x == Rsquared$Location.x[row] & contributions_df_joined$Location.y == Rsquared$Location.y[row]),]))$r.squared
	Rsquared$Line_Rsquared[row] <- summary(lm(Total_Contribution.y ~ Total_Contribution.x, contributions_df_joined[which(contributions_df_joined$Variable == "Line" & contributions_df_joined$Location.x == Rsquared$Location.x[row] & contributions_df_joined$Location.y == Rsquared$Location.y[row]),]))$r.squared
	Rsquared$Replicate_Rsquared[row] <- summary(lm(Total_Contribution.y ~ Total_Contribution.x, contributions_df_joined[which(contributions_df_joined$Variable == "Replicate" & contributions_df_joined$Location.x == Rsquared$Location.x[row] & contributions_df_joined$Location.y == Rsquared$Location.y[row]),]))$r.squared
}

Rsquared_long <- pivot_longer(Rsquared, cols = c("Village_Rsquared", "Line_Rsquared", "Replicate_Rsquared"), values_to = "Rsquared", names_to = "Variable")
Rsquared_long$Variable <- gsub("_Rsquared", "", Rsquared_long$Variable)
Rsquared_long$Rsquared <- as.character(round(Rsquared_long$Rsquared, 2))

Rsquared_long <- Rsquared_long[which(Rsquared_long$Location.x == "Fresh"),]
contributions_df_joined <- contributions_df_joined[which(contributions_df_joined$Location.x == "Fresh"),]

##### Figures #####
pCorr_Location_Point <- ggplot(contributions_df_joined, aes(Total_Contribution.x, Total_Contribution.y, color = Variable)) +
	geom_point(size = 0.5, alpha = 0.25) +
	# facet_grid(rows = vars(Location.y), cols = vars(Location.x)) +
	theme_classic() +
	scale_color_manual(values = variable_colors) +
	xlab("Contribution to Gene Variation - Fresh Samples") +
	ylab("Contribution to Gene Variation - Cryopreserved Samples") +
	geom_text(x = 0, y = 0.75, aes(label = paste0("R^2 = ", Rsquared)), data = Rsquared_long[which(Rsquared_long$Variable == "Line"),], hjust = 0) +
	geom_text(x = 0, y = 0.7, aes(label = paste0("R^2 = ", Rsquared)), data = Rsquared_long[which(Rsquared_long$Variable == "Village"),], hjust = 0) +
	geom_text(x = 0, y = 0.65, aes(label = paste0("R^2 = ", Rsquared)), data = Rsquared_long[which(Rsquared_long$Variable == "Replicate"),], hjust = 0)
save_figs(pCorr_Location_Point, paste0(outdir, "Correlation_Point_Location_combined"), width = 14, height = 11)



### Make a figure for the most line-specific gene: ENSG00000106 CHCHD2 ###
CHCHD2_list <- lapply(names(seurat_rb_list), function(x){
	df <- data.frame("Expression" = data.frame(seurat_rb_list[[x]][["SCT"]]@scale.data["ENSG00000106153",]), "Time" = seurat_rb_list[[x]]@meta.data$Time, "Line" = seurat_rb_list[[x]]@meta.data$Final_Assignment, "Replicate" = gsub("[A-Z][a-z]+", "", seurat_rb_list[[x]]@meta.data$MULTI_ID))
	colnames(df)[1] <- "Expression"
	df$Location <- x
	df$Village = ifelse(seurat_rb_list[[x]]@meta.data$Time == "Baseline", "Baseline", ifelse(seurat_rb_list[[x]]@meta.data$Time == "Thawed Village Day 0", "Baseline", "Village"))
	return(df)
})


CHCHD2_df <- do.call(rbind, CHCHD2_list)


CHCHD2_box <- ggplot(CHCHD2_df, aes(Village, Expression, color = Line)) +
				geom_boxplot(outlier.size = 0.5) +
				theme_classic() +
				facet_wrap(vars(Location), ncol = 1) +
				scale_color_manual(values = line_colors) +
				ylab("Normalized CHCHD2 Expression") +
				theme(axis.title.x=element_blank())
save_figs(CHCHD2_box, paste0(outdir, "CHCHD2_Expression_facet"), width = 12, height = 12)


##### Regress out other effects to compare line effect in each locationXtime #####
seurat <- readRDS(paste0(datadir,"seurat_integrated_all_times_clustered_rb_regressed.rds"))
seurat$Location <- gsub("_.+", "", seurat$Location_Time)
seurat_sub <- subset(seurat, subset = Location != "Brisbane")
seurat_sub <- subset(seurat_sub, subset = Location != "Melbourne")

seurat_list <- lapply(unique(seurat_sub@meta.data$Time), function(x){
	subset(seurat_sub, subset = Time == x)
})
names(seurat_list) <- unique(seurat_sub@meta.data$Time)


### Calculate the number of cells ###
seurat_list <- lapply(seurat_list, function(x){
	x$Ncorrection <- NA
	for (rep in unique(x$Site_rep)){
		print(rep)
		x$Ncorrection <- ifelse(x$Site_rep == rep, 1/nrow(x@meta.data[which(x$Site_rep == rep),]), x$Ncorrection)
	}
	return(x)
})


if (!file.exists(paste0(outdir, "seurat_SCT_wRB_sydney_regress_all_covs.rds"))){
	seurat_norm_list <- lapply(seurat_list, function(x){
		SCTransform(x, verbose = TRUE, vars.to.regress = c("scores.G1", "scores.S", "scores.G2M", "percent.mt", "percent.rb","MULTI_ID"), return.only.var.genes = FALSE)
	})
	saveRDS(seurat_norm_list, paste0(outdir, "seurat_SCT_wRB_sydney_regress_all_covs.rds"))
} else {
	seurat_norm_list <- readRDS(paste0(outdir, "seurat_SCT_wRB_sydney_regress_all_covs.rds"))
}



#### Filter out Genes expressed in low percentage of cells #####
pct_expressing <- lapply(seurat_norm_list, function(x){
	temp <- data.frame(rowSums(as.matrix(x[["SCT"]]@counts) > 0)/ncol(x[["SCT"]]@counts))
	colnames(temp) <- "Cell_proportion"
	return(temp)
})

### 10% to start
seurat_norm_list <- lapply(names(seurat_norm_list), function(x){
	subset(seurat_norm_list[[x]], features = rownames(pct_expressing[[x]])[which(pct_expressing[[x]]$Cell_proportion >= 0.1)])
})
names(seurat_norm_list) <- unique(seurat_sub@meta.data$Time)


sce_rb_separate_list <- lapply(seurat_norm_list, function(x){
	SingleCellExperiment(assays = list(counts = x[["SCT"]]@counts, logcounts = x[["SCT"]]@data, normcounts = x[["SCT"]]@scale.data), colData = x@meta.data)
})



if (!file.exists(paste0(outdir, "heirarchical_partitioning_results_sep_control_syd.rds"))){
	scran <- list()
	i <- 1
	for(x in names(seurat_norm_list)){
		print(i)
		df <- data.frame("Line" = seurat_norm_list[[x]]@meta.data$Final_Assignment)
		scran[[x]] <- getVarianceExplained(sce_rb_separate_list[[x]], variables = df, exprs_values = "normcounts")
		i <- 1 + i
	}
	saveRDS(scran, paste0(outdir, "heirarchical_partitioning_results_sep_control_syd.rds"))
} else {
	scran <- readRDS(paste0(outdir, "heirarchical_partitioning_results_sep_control_syd.rds"))
}



##### Make dataframes with all the values
contributions_cov_df_list <- lapply(names(scran), function(x){
	scran[[x]] <- as.data.frame(scran[[x]])
	scran[[x]]$Gene <- rownames(scran[[x]])
	scran[[x]]$Location_Time_Rep <- x
	return(scran[[x]])
})

contributions_cov_df <- do.call(rbind, contributions_cov_df_list)



##### Make some figures!!! #####
pTotal_Cont <- ggplot(contributions_cov_df, aes(Line)) +
	geom_density() +
	facet_wrap(Location_Time_Rep ~ ., ncol = 1) +
	theme_classic()  +
	# scale_y_continuous(trans = "log10") +
	scale_x_continuous(trans = "log10") +
	scale_color_manual(values = variable_colors[["Line"]]) +
	geom_vline(xintercept = 0.01, linetype="dashed")

save_figs(pTotal_Cont, paste0(outdir, "Total_Contribution_Histogram_cov"))



### Variance explained correlation between sites ###
contributions_cov_df_joined <- full_join(contributions_cov_df, contributions_cov_df, by = c("Gene"))
contributions_cov_df_joined <- contributions_cov_df_joined[which(contributions_cov_df_joined$Location_Time_Rep.x != contributions_cov_df_joined$Location_Time_Rep.y),]

Rsquared <- unique(contributions_cov_df_joined[,c("Location_Time_Rep.x","Location_Time_Rep.y")])
Rsquared$Line_Rsquared <- NA

for (row in 1:nrow(Rsquared)){
	Rsquared$Line_Rsquared[row] <- summary(lm(Line.y ~ Line.x, contributions_cov_df_joined[which(contributions_cov_df_joined$Location_Time_Rep.x == Rsquared$Location_Time_Rep.x[row] & contributions_cov_df_joined$Location_Time_Rep.y == Rsquared$Location_Time_Rep.y[row]),]))$r.squared
}



pCorr_Location_Point_Line <- ggplot(contributions_cov_df_joined, aes(Line.x, Line.y)) +
	geom_point(size = 0.5, alpha = 0.25, color = variable_colors["Line"]) +
	facet_grid(rows = vars(Location_Time_Rep.y), cols = vars(Location_Time_Rep.x)) +
	theme_classic() +
	xlab("Contribution to Gene Variation") +
	ylab("Contribution to Gene Variation") +
	geom_text(x = 0, y = 0.75, aes(label = paste0("R^2 = ", Line_Rsquared)), data = Rsquared, hjust = 0, color = variable_colors["Line"])
save_figs(pCorr_Location_Point_Line, paste0(outdir, "Correlation_Point_Location_Line_cov"))









