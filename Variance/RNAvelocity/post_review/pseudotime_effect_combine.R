library(tidyverse)
library(ggplot2)
library(ggrepel)
library(data.table)
library(Seurat)
library(viridis)
library(colorspace)
library(RColorBrewer)
library(ggsignif)
library(ggridges)



dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
datadir <- paste0(dir,"output/Variance/RNAvelocity/post_review/results/gene_separated/")
icc_dir <- paste0(datadir,"icc/")
icc_interaction_dir <- paste0(datadir,"icc_interaction/")
effect_interaction_dir <- paste0(datadir,"effect_betas/")
outdir <- paste0(dir,"output/Variance/RNAvelocity/post_review/results/combined/")
dir.create(outdir)




save_figs <- function(plot, basename, width = 17, height = 17, units = "cm"){
    ggsave(plot, filename = paste0(basename,".png"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".pdf"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".eps"), height = height, width = width, units = units)
}



##### Set up colors #####
variable_colors <- c(Village = "#a186aa", Replicate = "#f2c1ce", Line = "#4734a9") 
line_colors <- c(FSA0006 = "#F79E29", MBE1006 = "#9B2C99", TOB0421 = "#35369C")


site_updates <- c("Brisbane" = "Site 1", "Melbourne" = "Site 2" ,"Sydney" = "Site 3")


vars <- c("Line", "Village", "Cryopreserved", "Site", "Pseudotime", "Replicate","Line:Village", "Line:Cryopreserved", "Line:Site", "Line:Pseudotime", "Village:Cryopreserved","Village:Site",  "Replicate:Village", "Replicate:Line","Replicate:Cryopreserved",  "Replicate:Site", "Residual")
selected_vars <- c("Line", "Village", "Site",  "Replicate", "Line:Village", "Line:Site", "Village:Site", "Replicate:Village", "Replicate:Line", "Replicate:Site", "Residual")
# var_colors <- c(colorRampPalette(brewer.pal(11, "Spectral"))(14), "grey90")
var_colors <- c("#4734a9", rev(c("#115cc7", "#a4c3c8", "#499090", "#405940", "#685a54", "#f7d312","#85929E",   "#f8bf33", "#e4910e", "#f65d19", "#931519", "#935116", "#f2c1ce", "#e17aab", "#a186aa")), "gray80")


names(var_colors) <- vars


##### Add gene IDs for easy identification downstream #####
GeneConversion1 <- read_delim("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/DRENEA_1/outs/filtered_feature_bc_matrix/features.tsv.gz", col_names = F, delim = "\t")
GeneConversion2 <- read_delim("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/Village_B_1_week/outs/filtered_feature_bc_matrix/features.tsv.gz", col_names = F, delim = "\t")

GeneConversion <- unique(rbind(GeneConversion1, GeneConversion2))
GeneConversion <- GeneConversion[!duplicated(GeneConversion$X1),]
GeneConversion$X3 <- NULL
colnames(GeneConversion) <- c("ENSG_ID", "Gene_ID")




if (!file.exists(paste0(outdir,"variance_partitionin_df.tsv"))){
	##### Get list of files #####
	icc_files <- list.files(icc_dir)
	icc_interaction_files <- list.files(icc_interaction_dir)
	effect_interaction_files <- list.files(effect_interaction_dir)


	##### Read in the files #####
	icc_list <- lapply(icc_files, function(icc){
		readRDS(paste0(icc_dir,icc))
	})
	names(icc_list) <- gsub("_icc.rds", "", icc_files)


	icc_interaction_list <- lapply(icc_interaction_files, function(icc){
		readRDS(paste0(icc_interaction_dir,icc))
	})
	names(icc_interaction_list) <- gsub("_icc.rds", "", icc_interaction_files)

	## Takes too much memory to load all of them - will choose some to pull based on significance
	effect_interaction_list <- lapply(effect_interaction_files, function(icc){
		readRDS(paste0(effect_interaction_dir,icc))
	})
	names(effect_interaction_list) <- gsub("_effects.rds", "", effect_interaction_files)

	effect_interaction_list <- lapply(names(effect_interaction_list), function(x){
		effect_interaction_list[[x]]$Gene <- x
		return(effect_interaction_list[[x]])
	})


	##### Make dataframe of effect sizes from models #####
	effect_interaction_dt <- do.call(rbind, effect_interaction_list)

	effect_interaction_dt[grep(":Pseudotime", grp)][order(Effect)]


	##### Combine fits Results #####
	icc_df <- do.call(rbind, icc_interaction_list)

	write_delim(icc_df, paste0(outdir,"variance_partitionin_df.tsv"), delim = "\t")

} else {
	icc_df <- fread(paste0(outdir,"variance_partitionin_df.tsv"), sep = "\t")
}


##### Data.table for manuscript #####
var_explained_manuscript <- icc_df[,c("grp", "percent", "P", "gene")]

colnames(var_explained_manuscript) <- c("Covariate", "Percent_Variance_Explained", "P", "ENSG")



var_explained_manuscript <- data.table(GeneConversion)[var_explained_manuscript, on = c("ENSG_ID" = "ENSG")]

fwrite(var_explained_manuscript, paste0(outdir, "pseudotime_variance_explained_manuscript.tsv"), sep = "\t")



##### Check if any interactions are pluri genes #####
pluri_genes <- fread(paste0(dir,"data/pluripotency_genes.tsv"), sep = "\t")

icc_df[grp == "Line:Pseudotime"][pluri_genes, on = c("gene" = "ENSG")]



### Read in the seurat data objects
seurat_list <- list()
seurat_files <- list.files("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/RNAvelocity/nb_partitioning/data/")
for (file in seurat_files){
	seurat_list[[gsub("_seurat_1pct.rds", "", file)]] <- readRDS(paste0("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/RNAvelocity/nb_partitioning/data/",file))
}

seurat <- readRDS("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/RNAvelocity/post_review/data/seurat_integrated_all_times_clustered_1pct_expressing_pseudotime.rds")


##### Make long dataframe for plotting #####
# fits_df_long <- pivot_longer(fits_df, cols = c("Line", "Village", "Replicate"), names_to = "Covariate", values_to = "Variance_Explained")
# fits_df_long$Variance_Explained <- round(fits_df_long$Variance_Explained,6)

icc_df <- left_join(icc_df, GeneConversion, by = c("gene" = "ENSG_ID"))

icc_df$grp <- factor(icc_df$grp, levels = rev(c("Line", "Village", "Site", "Cryopreserved", "Replicate", "Pseudotime", "Line:Village", "Line:Site", "Line:Cryopreserved", "Line:Pseudotime", "Village:Site", "Village:Cryopreserved", "Replicate:Village", "Replicate:Cryopreserved", "Replicate:Line", "Replicate:Site", "Residual")))



##### Make some figures!!! #####
pRaincloud <- ggplot(icc_df, aes(x = percent, y = factor(grp, levels = rev(levels(grp))), fill = factor(grp, levels = rev(levels(grp))))) + 
                geom_density_ridges(size = 0.1,stat = "binline", bins = 100, scale = 0.7, draw_baseline = FALSE, aes(height =..ndensity..)) +
                geom_point(size =1, position = position_nudge(y=-0.11), shape = "|", aes(color = factor(grp, levels = rev(levels(grp))))) +
                coord_cartesian(xlim = c(1.2, NA), clip = "off") +
                theme_classic() +
                theme(axis.title.y=element_blank()) +
                xlab("Percent Variance Explained") +
                scale_y_discrete(expand = c(0.03, 0)) +
                scale_fill_manual(values = var_colors, name = "Variable") +
                scale_color_manual(values = var_colors, name = "Variable") +
                geom_vline(xintercept = 1, lty="11", color = "grey50", size = 0.5)


save_figs(pRaincloud, paste0(outdir, "Total_Contribution_Histogram_cov"))




quintile_list <- lapply(seurat_list, function(x){
	unique(data.table(Barcode = rownames(x@meta.data), Quintile = as.numeric(as.character(gsub("Q","", x@meta.data$Quintile))), Location = gsub("Brisbane", "Site 1", x@meta.data$Location) %>% gsub("Melbourne", "Site 2", .) %>% gsub("Sydney", "Site 3", .) %>% gsub("_", " ", .), Expression = x[["SCT"]]@counts["ENSG00000106153",], Normalized_Expression = x[["SCT"]]@data["ENSG00000106153",], Scaled_Expression = x[["SCT"]]@scale.data["ENSG00000106153",], Line = x@meta.data$Final_Assignment, Time = x@meta.data$Time, Pseudotime = x@meta.data$latent_time))
})

quintile_dt <- do.call(rbind, quintile_list)

pseudo_dt <- data.table(Barcode = rownames(seurat@meta.data), Location = gsub("Brisbane", "Site 1", seurat@meta.data$Location) %>% gsub("Melbourne", "Site 2", .) %>% gsub("Sydney", "Site 3", .) %>% gsub("_", " ", .), Expression = seurat[["SCT"]]@counts["ENSG00000106153",], Normalized_Expression = seurat[["SCT"]]@data["ENSG00000106153",], Scaled_Expression = seurat[["SCT"]]@scale.data["ENSG00000106153",], Line = seurat@meta.data$Final_Assignment, Time = seurat@meta.data$Time, Pseudotime = seurat@meta.data$latent_time, Cryopreserved = seurat@meta.data$Cryopreserved)

pseudo_dt <- pseudo_dt[quintile_dt[,c("Barcode", "Quintile")], on = "Barcode"]

pseudo_dt$Location <- ifelse(pseudo_dt$Cryopreserved == "Cryopreserved", paste0(pseudo_dt$Location, " Cryopreserved"), pseudo_dt$Location)

pseudo_dt$box_location <- as.numeric(NA)


for (quin in unique(pseudo_dt$Quintile)){
	pseudo_dt[Quintile == quin]$box_location <- mean(pseudo_dt[Quintile == quin]$Pseudotime)
}



pExpression <- ggplot(pseudo_dt, aes(x = factor(Quintile, levels = seq(1:5)), y = Expression, color = Line, fill = Line)) +
	geom_boxplot(outlier.size = 0.25, width = 0.5, position=position_dodge(0.7)) +
	theme_classic() +
	facet_wrap(vars(Location), ncol = 1, scales = "free_y") +
	scale_color_manual(values = line_colors) +
	scale_fill_manual(values = alpha(line_colors, 0.5)) +
	ylab(paste0("CHCHD2 Expression")) +
	xlab("Quintile")


ggsave(pExpression, filename = paste0(outdir, "CHCHD2_quintile_plot.png"), width = 5, height = 7)
ggsave(pExpression, filename = paste0(outdir, "CHCHD2_quintile_plot.pdf"), width = 5, height = 7)


pExpression <- ggplot(pseudo_dt, aes(x = factor(Quintile, levels = seq(1:5)), y = Scaled_Expression, color = Line, fill = Line)) +
	geom_boxplot(outlier.size = 0.25, width = 0.5, position=position_dodge(0.7)) +
	theme_classic() +
	facet_wrap(vars(Location), ncol = 1, scales = "free_y") +
	scale_color_manual(values = line_colors) +
	scale_fill_manual(values = alpha(line_colors, 0.5)) +
	ylab(paste0("CHCHD2 Expression")) +
	xlab("Quintile")


ggsave(pExpression, filename = paste0(outdir, "CHCHD2_quintile_plot_normalized.png"), width = 5, height = 7)
ggsave(pExpression, filename = paste0(outdir, "CHCHD2_quintile_plot_normalized.pdf"), width = 5, height = 7)






plot_facet <- ggplot(pseudo_dt, aes(Pseudotime, Expression, color = Line)) +
			# geom_point(alpha = 0.05)  +
			facet_wrap(vars(Location), ncol = 1, scales = "free_y") +
			geom_smooth(se = TRUE) +
			# geom_smooth(method = "lm", se = TRUE) +
			theme_classic() +
			scale_color_manual(values = line_colors)

ggsave(plot_facet, filename = paste0(outdir, "CHCHD2_expression_pseudotime_line.png"), width = 4, height = 7)


plot_facet_normalized <- ggplot(pseudo_dt) +
			# geom_boxplot(outlier.size = 0.25, width = 0.5, position=position_dodge(0.7), aes(x = factor(box_location), y = Scaled_Expression, color = Line, fill = Line)) +
			geom_point(aes(x = Pseudotime, y = Scaled_Expression, color = Line, fill = Line),alpha = 0.15, size = 0.1)  +
			facet_wrap(vars(Location), ncol = 1, scales = "free_y") +
			geom_smooth(aes(Pseudotime, Scaled_Expression, group = Line), color = "black", size = 0.95, se = FALSE) +
			geom_smooth(aes(Pseudotime, Scaled_Expression, color = Line), fill = "black", size = 0.5, se = TRUE) +
			# geom_smooth(method = "lm", se = TRUE) +
			theme_classic() +
			scale_color_manual(values = line_colors) +
			ylab("CHCHD2 Normalized Expression")

ggsave(plot_facet_normalized, filename = paste0(outdir, "CHCHD2_expression_pseudotime_line_normalized.png"), width = 3.5, height = 5)
ggsave(plot_facet_normalized, filename = paste0(outdir, "CHCHD2_expression_pseudotime_line_normalized.pdf"), width = 3.5, height = 5)

