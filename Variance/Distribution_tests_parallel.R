library(ggplot2)
library(tidyverse)
library(Seurat)
library(gamlss)
library(fitdistrplus)
library(VGAM)
library(emdbook)


##### Bring in variables #####
### Bring in arguments
args <- commandArgs(trailingOnly = TRUE)
outdir <- paste0(args[1])
number <- as.numeric(args[2])



##### Set Up Directories #####
dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
datadir <- paste0(dir,"output/Distribution_tests/")

dir.create(outdir)



##### Redo with all the cells across all times #####
##### Obtain and manage data #####
seurat_sub <- readRDS(paste0(datadir,"seurat_integrated_all_times_clustered_1pct_expressing.rds"))
gene <- rownames(seurat_sub[["SCT"]]@counts)[number]


aic_df <- as.data.frame(matrix(nrow = 0, ncol = 3))
colnames(aic_df) <- c("gene", "normal", "negative_binomial")

fit <- list()


fit[["norm"]] <- fitdist(seurat_sub[["SCT"]]@scale.data[number,], "norm")
fit[["nbinom"]] <- fitdist(seurat_sub[["SCT"]]@counts[number,], "nbinom")
aic_df <- rbind(aic_df, data.frame("gene" = rownames(seurat_sub[["SCT"]])[number], "normal" = fit[["norm"]]$aic, "negative_binomial" = fit[["nbinom"]]$aic))


aic_df$difference <- aic_df$normal - aic_df$negative_binomial
aic_df$better_model <- ifelse(aic_df$normal < aic_df$negative_binomial, "normal","negative_binomial")

saveRDS(fit, paste0(outdir,gene,"_fits.rds"))
write_delim(aic_df, paste0(outdir,gene,"_aic.tsv"), delim = "\t")

