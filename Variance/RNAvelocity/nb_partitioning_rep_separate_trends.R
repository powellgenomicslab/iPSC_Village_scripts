library(tidyverse)
library(data.table)
library(ComplexHeatmap)
library(Seurat)
library(ggpubr)
library(grid)


##### Set up Directories #####
datadir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/RNAvelocity/nb_partitioning_rep_separate/combined/"
seuratdir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/RNAvelocity/nb_partitioning_rep_separate/data/"
outdir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/RNAvelocity/nb_partitioning_rep_separate/trends/"
dir.create(outdir, recursive = TRUE)

site_updates <- c("Brisbane" = "Site 1", "Sydneycryopreserved" = "Site 3 Cryopreserved", "Sydney" = "Site 3" ,"Melbourne" = "Site 2")
line_colors <- c(FSA0006 = "#F79E29", MBE1006 = "#9B2C99", TOB0421 = "#35369C")
site_colors <- c("Site 1" = "#C9D8EA", "Site 2" = "#5D59AB", "Site 3 Cryopreserved" = "#A7C9A9", "Site 3" = "#179085")


##### Read in Data #####
fits_df_long <- fread(sep = "\t", paste0(datadir, "var_explained_long.tsv"))



##### Reformat Line dataframe for Trend dettecting #####
fits_df_long_line <- fits_df_long[Covariate == "Line"]
fits_df_long_line$Quintile <- as.numeric(as.character(gsub("Q", "", fits_df_long_line$Quintile)))
fits_df_long_line <- unique(fits_df_long_line)
fits_df_long_line$Location <- gsub("Site 3", "Site 4", fits_df_long_line$Location) %>% gsub("Site 2", "Site 3", .) %>% gsub("Site 4", "Site 2", .)


##### Calculate for linear trends in data #####
fits <- list()
data_frames <- list()
fits_combined <- list()
data_frames_combined <- list()
df_combined2 <- list()
i <- 1

for (gene in unique(fits_df_long_line$Gene)){
	print(i)
	df_combined <- fits_df_long_line[Gene == gene]
	# for (location in unique(fits_df_long_line$Location)){
	# 	df <- fits_df_long_line[Location == location & Gene == gene]
	# 	if (nrow(df) >= 3 & length(unique(df$Quintile)) >= 2){
	# 		fits[[gene]][[location]] <- tryCatch(summary(lm(Variance_Explained ~ Quintile, df)), error=function(e) NULL, warning=function(w) NULL)
	# 		data_frames[[gene]][[location]] <- data.table(Gene = gene, Gene_ID = unique(fits_df_long_line[Gene == gene]$Gene_ID), Location = location, P = fits[[gene]][[location]]$coefficients[2,4], Beta = fits[[gene]][[location]]$coefficients[2,1])
	# 	}
	# }
	# df_combined2[[gene]] <- data.table(Residuals = resid(tryCatch(lm(Variance_Explained ~ Location, df_combined), error=function(e) NULL, warning=function(w) NULL)), Quintile = df_combined$Quintile, Location = df_combined$Location)
	if (nrow(df_combined2[[gene]]) >= 3 & length(unique(df_combined2[[gene]]$Quintile)) >= 2){
		# fits_combined[[gene]] <- tryCatch(summary(lm(Residuals ~ Quintile, df_combined2[[gene]])), error=function(e) NULL, warning=function(w) NULL)
		data_frames_combined[[gene]] <- data.table(Gene = gene, Gene_ID = unique(fits_df_long_line[Gene == gene]$Gene_ID), P = fits_combined[[gene]]$coefficients[2,4], Beta = fits_combined[[gene]]$coefficients[2,1], SE = fits_combined[[gene]]$coefficients[2,2])
	}
	i <- i + 1
}


##### Make dataframe of linear trends #####
fit_df_list <- lapply(data_frames, function(x){
	do.call(rbind, c(x, fill=TRUE))
})

fit_df <- data.table(do.call(rbind, c(fit_df_list, fill=TRUE)))

fit_df_combined_list <- data.table(do.call(rbind, c(data_frames_combined, fill=TRUE)))

## Remove the rows that have Na for Beta and/or P value
fit_df <- fit_df[!(is.na(P) & is.na(Beta))]
fit_df <- fit_df[!(is.na(P))]

fit_df_combined_list <- fit_df_combined_list[!(is.na(P) & is.na(Beta))]
fit_df_combined_list <- fit_df_combined_list[!(is.na(P))]

fit_df <- fit_df[order(fit_df$P),]
fit_df_combined_list <- fit_df_combined_list[order(fit_df_combined_list$P),]
fit_df_combined_list$fdr <- p.adjust(fit_df_combined_list$P, method = "fdr")


### Any significant gene
fit_df_combined_sig <- fit_df_combined_list[fdr < 0.05] ### Decide to use this - where the locations are combined in the same model but the location effect is regressed out

fwrite(fit_df_combined_sig, paste0(outdir, "significant_pseudotime_dynamic_genes.tsv"), sep = "\t")
fit_df_combined_sig <- fread(paste0(outdir, "significant_pseudotime_dynamic_genes.tsv"), sep = "\t")


##### Plot the Significant Results #####
### Across all sites ###
mat_sig <- dcast(fit_df_combined_sig, . ~ Gene_ID, value.var = "Beta")
mat_sig$. <- NULL

mat_sig <- as.matrix(mat_sig)


sig_all_heatmap <- Heatmap(mat_sig, na_col = "grey", show_column_dend = FALSE, show_row_dend = FALSE, row_dend_reorder = FALSE)


pdf(paste0(outdir,"Heatmap_all_sig.pdf"), width = 20, height = 2)
draw(sig_all_heatmap)
dev.off()


png(paste0(outdir,"Heatmap_all_sig.png"), width = 12000, height = 1500, res = 600)
draw(sig_all_heatmap)
dev.off()


### Significant with absolute effect > 0.05 ###
fit_df_combined_sig_05 <- fit_df_combined_sig[abs(Beta) > 0.05]

mat_sig_05 <- dcast(fit_df_combined_sig_05, . ~ Gene_ID, value.var = "Beta")
mat_sig_05$. <- NULL

mat_sig_05 <- as.matrix(mat_sig_05)


sig_heatmap_05 <- Heatmap(mat_sig_05, na_col = "grey", show_column_dend = FALSE, show_row_dend = FALSE)


pdf(paste0(outdir,"Heatmap.pdf"), width = 10, height = 1.75)
draw(sig_heatmap_05)
dev.off()


png(paste0(outdir,"Heatmap.png"), width = 5000, height = 250, res = 600)
draw(sig_heatmap_05)
dev.off()



##### Try a Dot Plot with  SE to visualize #####
fit_df_combined_sig_05_ordered <- fit_df_combined_sig_05[order(fit_df_combined_sig_05$Beta)]
fit_df_combined_sig_05_ordered$Gene_ID <- factor(fit_df_combined_sig_05_ordered$Gene_ID, levels = fit_df_combined_sig_05_ordered$Gene_ID)

pDot <- ggplot(fit_df_combined_sig_05_ordered, aes(Beta, Gene_ID)) +
		geom_errorbar(aes(xmin=Beta-SE, xmax=Beta+SE), colour="black", width=0) +
		geom_point(size = 2) +
		theme_classic() +
		# scale_colour_gradient2(low = "blue",mid = "white",high = "red") +
		geom_vline(xintercept = 0, linetype = "dashed") +
		ylab("Gene")

ggsave(pDot, filename = paste0(outdir, "Beta_dot_plot.png"), width = 3, height = 4.25)
ggsave(pDot, filename = paste0(outdir, "Beta_dot_plot.pdf"), width = 3, height = 4.25)



# ##### Make plots of variance
# pD21S2088E <- ggplot(fits_df_long_line[Gene_ID == "D21S2088E"], aes(Quintile, Variance_Explained)) +
# 		geom_point(size = 0.7, alpha = 0.7) +
# 		theme_classic() +
# 		facet_wrap(vars(Location), scales = "free", nrow = 1)
# ggsave(pD21S2088E, filename = paste0(outdir, "D21S2088E_variance_explained.png"), width = 8, height = 2)


##### Read in seurat objects for expression comparison
seurat_files <- list.files(seuratdir, pattern = ".rds")

seurat_list <- lapply(seurat_files, function(x){
	readRDS(paste0(seuratdir,x))
})

seurat <- merge(seurat_list[[1]], seurat_list[2:length(seurat_list)])

head(seurat@meta.data)


##### Update Site names #####
for (location in names(site_updates)){
	print(location)
	seurat@meta.data$Location <- gsub(location, site_updates[location], seurat@meta.data$Location)
}

seurat@meta.data$Time <- gsub("Village Day 4", "Village", seurat@meta.data$Time)

gene_expression_dir <- paste0(outdir,"var_gene_quintile_expression/")
dir.create(gene_expression_dir)

for (gene_id in colnames(mat_sig_05)){
	gene <- unique(fit_df_combined_sig[Gene_ID == gene_id]$Gene)
	df <- data.table(Quintile = as.numeric(as.character(gsub("Q","", seurat@meta.data$Quintile))), Location = seurat@meta.data$Location, Expression = seurat[["SCT"]]@counts[gene,], Line = seurat@meta.data$Final_Assignment, Time = seurat@meta.data$Time)

	pVar <- ggplot(df_combined2[[gene]], aes(Quintile, Residuals, color = Location)) +
			geom_point() +
			theme_classic() +
			scale_color_manual(values = site_colors) +
			# facet_wrap(vars(Location), scales = "free", nrow = 1) +
			geom_smooth(method = "lm", aes(group = 1), color = "black", size=0.5) +
			ylab("Normalized\nVariance Explained")

	pExpression <- ggplot(df, aes(x = factor(Quintile, levels = seq(1:5)), y = Expression, fill = Line)) +
		geom_boxplot(outlier.size = 0.5) +
		theme_classic() +
		facet_wrap(vars(Location), nrow = 1, scales = "free") +
		# facet_grid(Time ~ Location, scales = "free") +
		scale_fill_manual(values = line_colors) +
		ylab(gene_id) +
		xlab("Quintile") 
		# theme(legend.position="bottom")

	pCombined <- ggarrange(pVar, pExpression, ncol = 1, nrow = 2, heights = c(2,1))
	ggsave(pCombined, filename = paste0(gene_expression_dir,gene_id, "_variance_expression.png"), height = 7, width = 8)
}



##### Final figure fo main fig #####
pVar <- ggplot(df_combined2[["ENSG00000106538"]], aes(x = factor(Quintile), y = Residuals, color = Location)) +
		geom_point(size = 1) +
		theme_classic() +
		scale_color_manual(values = site_colors) +
		# facet_wrap(vars(Location), scales = "free", nrow = 1) +
		geom_smooth(method = "lm", aes(group = 1), color = "black", size=0.5) +
		ylab("Normalized\nVariance Explained") +
		xlab("Quintile")


df <- data.table(Quintile = as.numeric(as.character(gsub("Q","", seurat@meta.data$Quintile))), Location = seurat@meta.data$Location, Expression = seurat[["SCT"]]@counts["ENSG00000106538",], Line = seurat@meta.data$Final_Assignment, Time = seurat@meta.data$Time)

pExpression <- ggplot(df, aes(x = factor(Quintile, levels = seq(1:5)), y = Expression, color = Line, fill = Line)) +
	geom_boxplot(outlier.size = 0.25, width = 0.5, position=position_dodge(0.7)) +
	theme_classic() +
	facet_wrap(vars(Location), ncol = 1, scales = "free_y") +
	scale_color_manual(values = line_colors) +
	scale_fill_manual(values = alpha(line_colors, 0.5)) +
	ylab(paste0(gene_id, " Expression")) +
	xlab("Quintile")

pCombined <- ggarrange(pVar, pExpression, heights = c(1.5, 4),
          ncol = 1, nrow = 2, align = "v")

ggsave(pCombined, filename = paste0(outdir, "RARRES2_combined_plot.png"), width = 5, height = 7)
ggsave(pCombined, filename = paste0(outdir, "RARRES2_combined_plot.pdf"), width = 5, height = 7)


AHR_df <- data.table(Quintile = as.numeric(as.character(gsub("Q","", seurat@meta.data$Quintile))), Location = seurat@meta.data$Location, Expression = seurat[["SCT"]]@counts["ENSG00000106546",], Line = seurat@meta.data$Final_Assignment, Time = seurat@meta.data$Time)

pExpression_AHR <- ggplot(AHR_df, aes(x = factor(Quintile, levels = seq(1:5)), y = Expression)) +
	geom_boxplot(outlier.size = 0.25, width = 0.5, position=position_dodge(0.7)) +
	theme_classic() +
	facet_wrap(vars(Location), ncol = 1, scales = "free_y") +
	ylab(paste0("AHR Expression")) +
	xlab("Quintile")

ggsave(pExpression_AHR, filename = paste0(outdir, "AHR_exression.png"), width = 5, height = 7)
ggsave(pExpression_AHR, filename = paste0(outdir, "AHR_exression.pdf"), width = 5, height = 7)