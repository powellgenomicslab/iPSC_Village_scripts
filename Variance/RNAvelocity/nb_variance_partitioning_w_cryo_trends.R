library(tidyverse)
library(data.table)
library(ComplexHeatmap)



##### Set up Directories #####
datadir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/RNAvelocity/nb_partitioning/combined/"
outdir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/RNAvelocity/nb_partitioning/trends/"
dir.create(outdir, recursive = TRUE)



##### Read in Data #####
fits_df_long <- fread(sep = "\t", paste0(datadir, "var_explained_long.tsv"))



##### Reformat Line dataframe for Trend dettecting #####
fits_df_long_line <- fits_df_long[Covariate == "Line"]
fits_df_long_line$Quartile <- as.numeric(as.character(gsub("Q", "", fits_df_long_line$Quartile)))
fits_df_long_line <- unique(fits_df_long_line)


##### Calculate for linear trends in data #####
fits <- list()
data_frames <- list()
i <- 1

for (gene in unique(fits_df_long_line$Gene)){
	print(i)
	for (location in unique(fits_df_long_line$Location)){
		df <- fits_df_long_line[Location == location & Gene == gene]
		if (nrow(df) >= 3){
			fits[[gene]][[location]] <- summary(lm(Variance_Explained ~ Quartile, df))
			data_frames[[gene]][[location]] <- data.table(Gene = gene, Gene_ID = unique(fits_df_long_line[Gene == gene]$Gene_ID), Location = location, P = fits[[gene]][[location]]$coefficients[2,4], Beta = fits[[gene]][[location]]$coefficients[2,1])
		}
	}
	i <- i + 1
}


##### Make dataframe of linear trends #####
fit_df_list <- lapply(data_frames, function(x){
	do.call(rbind, x)
})


fit_df <- data.table(do.call(rbind, fit_df_list))


fit_df <- fit_df[order(fit_df$P),]


### Any significant gene
fit_df_sig <- fit_df[P < 0.05]
fit_df_1sig <- fit_df_sig[Gene %in% names(which(table(fit_df_sig$Gene) > 2))]

mat_sig <- dcast(fit_df_1sig, Location ~ Gene_ID, value.var = "Beta")
locations <- mat_sig$Location
rownames(mat_sig) <- locations
mat_sig$Location <- NULL

mat_sig <- as.matrix(mat_sig)
rownames(mat_sig) <- locations
rownames(mat_sig) <- c("Site 1", "Site 2 - Fresh", "Site 2 - Cryopreserved", "Site 3")

sig_any_heatmap <- Heatmap(mat_sig, na_col = "grey", show_column_dend = FALSE, show_row_dend = FALSE, row_dend_reorder = FALSE, column_dend_reorder = FALSE)


pdf(paste0(outdir,"Heatmap_3_sig.pdf"), width = 10, height = 5)
draw(sig_any_heatmap)
dev.off()


png(paste0(outdir,"Heatmap_3_sig.png"), width = 4000, height = 1500, res = 600)
draw(sig_any_heatmap)
dev.off()








##### Plot the Significant Results #####
fit_df_sig <- fit_df[P < 0.05]


### Across all sites ###
fit_df_all_sig <- fit_df_sig[Gene %in% names(which(table(fit_df_sig$Gene) > 3))]


mat_sig <- dcast(fit_df_all_sig, Location ~ Gene_ID, value.var = "Beta")
locations <- mat_sig$Location
rownames(mat_sig) <- locations
mat_sig$Location <- NULL

mat_sig <- as.matrix(mat_sig)
rownames(mat_sig) <- locations






## Get just the genes whose diretion is the same for all sites
mat_sig_all <- mat_sig_all[,c(names(which(colSums(mat_sig_all > 0, na.rm = TRUE) == colSums(!is.na(mat_sig_all)))), names(which(colSums(mat_sig_all < 0, na.rm = TRUE) == colSums(!is.na(mat_sig_all)))))]
mat_sig_all <- mat_sig_all[c("Site 1", "Site 2", "Site 2Cryopreserved", "Site 3"), ]
rownames(mat_sig_all) <- c("Site 1", "Site 2 - Fresh", "Site 2 - Cryopreserved", "Site 3")


sig_all_heatmap <- Heatmap(mat_sig_all, na_col = "grey", show_column_dend = FALSE, show_row_dend = FALSE, row_dend_reorder = FALSE)


pdf(paste0(outdir,"Heatmap_all_sig.pdf"), width = 10, height = 5)
draw(sig_all_heatmap)
dev.off()


png(paste0(outdir,"Heatmap_all_sig.png"), width = 4000, height = 1500, res = 600)
draw(sig_all_heatmap)
dev.off()


### Across three sites ###
fit_df_3sig <- fit_df_sig[Gene %in% names(which(table(fit_df_sig$Gene) > 2))]


mat <- dcast(fit_df_3sig, Location ~ Gene_ID, value.var = "Beta")
locations <- mat$Location
rownames(mat) <- locations
mat$Location <- NULL

mat <- as.matrix(mat)
rownames(mat) <- locations

## Get just the genes whose diretion is the same for all sites
mat <- mat[,c(names(which(colSums(mat > 0, na.rm = TRUE) == colSums(!is.na(mat)))), names(which(colSums(mat < 0, na.rm = TRUE) == colSums(!is.na(mat)))))]

sig_heatmap <- Heatmap(mat, na_col = "grey", show_column_dend = FALSE, show_row_dend = FALSE)


pdf(paste0(outdir,"Heatmap.pdf"), width = 15, height = 5)
draw(sig_heatmap)
dev.off()


png(paste0(outdir,"Heatmap.png"), width = 10000, height = 1500, res = 600)
draw(sig_heatmap)
dev.off()