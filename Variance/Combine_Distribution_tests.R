library(tidyverse)
library(hier.part.negbinom)
library(Seurat)


dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
datadir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Distribution_tests_parallel/gene_separated/"
outdir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Distribution_tests_parallel/combined/"
dir.create(outdir)


##### Get list of files #####
fits_files <- list.files(datadir, pattern = "_fits.rds")
names(fits_files) <- gsub("_fits.rds", "", fits_files)

# ##### Read in files #####
# fits_list <- lapply(fits_files, function(x){
# 	readRDS(paste0(datadir,x))
# })


if (!file.exists(paste0(outdir,"aic_df"))){
	aic_files <- list.files(datadir, pattern = "_aic.tsv")
	names(aic_files) <- gsub("_aic.tsv", "", aic_files)


	aic_list <- lapply(aic_files, function(x){
		read_delim(paste0(datadir,x), delim = "\t")
	})


	##### Combine AIC Results #####
	aic_df <- do.call(rbind, aic_list)
	aic_df$Gene <- names(aic_list)

	write_delim(aic_df, paste0(outdir,"aic_df"))

} else {
	aic_df <- read_delim(paste0(outdir,"aic_df"))
}


##### Take a look at the results #####
table(aic_df$better_model)

max(aic_df$difference[which(aic_df$better_model == "normal")])

max(aic_df$difference[which(aic_df$better_model == "negative_binomial")])



##### Read in seurat with genes #####
seurat_sub <- readRDS(paste0(dir,"output/Distribution_tests/seurat_integrated_all_times_clustered_1pct_expressing.rds"))


gene <- 58
df <- data.frame("Expression" = data.frame(seurat_sub[["SCT"]]@counts[gene,]), "Village" = ifelse(seurat_sub@meta.data$Time == "Baseline", 0, 1), "Line" = seurat_sub@meta.data$Final_Assignment, "Replicate" = gsub("[A-Z][a-z]+", "", seurat_sub@meta.data$MULTI_ID))
colnames(df)[1] <- "Expression"

hier_models <- hier.part.negbinom::hier.part(y = df$Expression, xcan = df[,2:ncol(df)], family = "negbinom",gof = "logLik")
hier_models_rmspe <- hier.part.negbinom::hier.part(y = df$Expression, xcan = df[,2:ncol(df)], family = "negbinom",gof = "RMSPE")
# hier_models <- hier.part.negbinom::hier.part(y = df$Expression, xcan = df[,2:ncol(df)], gof = "logLik")

df2 <- data.frame("Expression" = data.frame(seurat_sub[["SCT"]]@counts["ENSG00000106153",]), "Village" = ifelse(seurat_sub@meta.data$Time == "Baseline", 0, 1), "Line" = seurat_sub@meta.data$Final_Assignment, "Replicate" = gsub("[A-Z][a-z]+", "", seurat_sub@meta.data$MULTI_ID))
colnames(df2)[1] <- "Expression"
CHCHD2_log <- hier.part.negbinom::hier.part(y = df2$Expression, xcan = df2[,2:ncol(df2)], family = "negbinom",gof = "logLik")
CHCHD2_rmspe <- hier.part.negbinom::hier.part(y = df2$Expression, xcan = df2[,2:ncol(df2)], family = "negbinom",gof = "RMSPE")


df2_factor <- df2
df2_factor$Village <- as.factor(df2$Village)
CHCHD2_log_fact <- hier.part.negbinom::hier.part(y = df2_factor$Expression, xcan = df2_factor[,2:ncol(df2_factor)], family = "negbinom",gof = "logLik", link = "log")