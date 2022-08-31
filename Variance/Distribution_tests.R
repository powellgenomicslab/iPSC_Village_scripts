library(ggplot2)
library(tidyverse)
library(Seurat)
library(gamlss)
library(fitdistrplus)
library(VGAM)
library(emdbook)
library(data.table)


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
outdir <- paste0(dir,"output/Distribution_tests/")
datadir <- paste0(dir,"output/All_data_integrated_remove_bad/")


dir.create(outdir)


##### Originally used the separated seurat objects by condition but need to use the same model for each condition so will just detect the model type across all conditions #####
# location_SCT <- readRDS(paste0(dir,"output/Variance_contributions_lm_scale_data/seurat_SCT_wRB_sydney_regress_all_covs.rds"))
# cryo_SCT <- readRDS(paste0(dir,"output/Variance_contributions_lm_scale_data_Sydney/seurat_SCT_wRB_sydney_regress_all_covs.rds"))


# ### Update names to combine objects ###
# names(location_SCT) <- c("Baseline", "Village")
# location_SCT <- lapply(location_SCT, function(x){
# 	names(x) <- c("Brisbane", "Melbourne", "Sydney_Fresh")
# 	return(x)
# })


# cryo_SCT_updated <- list()
# cryo_SCT_updated[["Baseline"]][["Sydney_Fresh"]] <- cryo_SCT[["Baseline"]]
# cryo_SCT_updated[["Baseline"]][["Sydney_Cryopreserved"]] <- cryo_SCT[["Thawed Village Day 0"]]
# cryo_SCT_updated[["Village"]][["Sydney_Fresh"]] <- cryo_SCT[["Village Day 4"]]
# cryo_SCT_updated[["Village"]][["Sydney_Cryopreserved"]] <- cryo_SCT[["Thawed Village Day 7"]]


# SCT_combined <- location_SCT
# SCT_combined[["Baseline"]][["Sydney_Cryopreserved"]] <- cryo_SCT_updated[["Baseline"]][["Sydney_Cryopreserved"]]
# SCT_combined[["Village"]][["Sydney_Cryopreserved"]] <- cryo_SCT_updated[["Village"]][["Sydney_Cryopreserved"]]



#### Filter out Genes expressed in low percentage of cells #####
pct_zero <- lapply(SCT_combined, function(x){
	lapply(x, function(y){
		temp <- data.frame(rowSums(as.matrix(y[["SCT"]]@counts) == 0)/ncol(y[["SCT"]]@counts))
		colnames(temp) <- "Cell_proportion"
		return(temp)
	})
})


lapply(pct_zero, function(x) lapply(x, head))


pct_zero_df_list <- lapply(names(pct_zero), function(x){
	temp <- lapply(names(pct_zero[[x]]), function(y){
		pct_zero[[x]][[y]]$Condition <- y
		return(pct_zero[[x]][[y]])
	})
	temp2 <- do.call(rbind, temp)
	print(head(temp2))
	temp2$Time <- x
	print(head(temp2))
	return(temp2)
})

pct_zero_df <- do.call(rbind, pct_zero_df_list)


figures <- ggplot(pct_zero_df, aes(Cell_proportion)) +
		geom_histogram() +
		theme_classic() +
		facet_grid(Time ~ Condition)


save_figs(figures, paste0(outdir, "percent_zeros_histogram.png"))


# ##### Use fitDist from gamlss package to compare different distributions for the data and select the best - test first #####
# realAll_fit_list <- lapply(SCT_combined, function(x){
# 	lapply(x, function(y){
# 		fitDist(as.matrix(y[["SCT"]]@scale.data[1089,]), k = 2, type = "realAll")
# 	})
# })



# counts_fit_list <- lapply(SCT_combined, function(x){
# 	lapply(x, function(y){
# 		fitDist(as.matrix(y[["SCT"]]@counts[1089,]), k = 2, type = "counts")
# 	})
# })

# ### conclusion: the best fit for each condition is different for the same gene



# realAll_fit[["fits"]]
# counts_fit[["fits"]]

# ref_dist <- rGT(ncol(SCT_combined[[1]][[1]][["SCT"]]), mu = realAll_fit[["mu"]], sigma = realAll_fit[["sigma"]], nu = realAll_fit[["nu"]], tau = realAll_fit[["tau"]])
# plot_scaled <- ggplot(data.frame("Counts" = SCT_combined[[1]][[1]][["SCT"]]@scale.data[1089,]), aes(Counts)) +
# 					geom_histogram(binwidth = 1) +
# 					theme_classic() +
# 					geom_density(data = data.frame("Counts" = ref_dist), aes(y=1 * ..count..))
# ggsave(plot_scaled, filename = paste0(outdir,"test_dist_histo_scale_data.png"))

# qq_ref_dist <- ggplot(data.frame("Empirical" = SCT_combined[[1]][[1]][["SCT"]]@scale.data[1089,][order(SCT_combined[[1]][[1]][["SCT"]]@scale.data[1089,])], "Simulated" = ref_dist[order(ref_dist)]), aes(Empirical, Simulated)) +
# 					geom_point() +
# 					geom_abline(intercept = 0, slope = 1) +
# 					theme_classic()
# ggsave(qq_ref_dist, filename = paste0(outdir,"empirical_simulated_qq_scale_data.png"))


# ref_count_dist <- rZIP(n=ncol(SCT_combined[[1]][[1]][["SCT"]]), mu = counts_fit[["mu"]], sigma = counts_fit[["sigma"]])
# plot_count <- ggplot(data.frame("Counts" = SCT_combined[[1]][[1]][["SCT"]]@counts[1089,]), aes(Counts)) +
# 					geom_histogram(binwidth = 1) +
# 					theme_classic() +
# 					geom_density(data = data.frame("Counts" = ref_count_dist), aes(y=1 * ..count..))
# ggsave(plot_count, filename = paste0(outdir,"test_dist_histo_count.png"))


# qq_ref_count_dist <- ggplot(data.frame("Empirical" = SCT_combined[[1]][[1]][["SCT"]]@counts[1089,][order(SCT_combined[[1]][[1]][["SCT"]]@counts[1089,])], "Simulated" = ref_count_dist[order(ref_count_dist)]), aes(Empirical, Simulated)) +
# 					geom_point() +
# 					geom_abline(intercept = 0, slope = 1) +
# 					theme_classic()
# ggsave(qq_ref_count_dist, filename = paste0(outdir,"empirical_simulated_qq_count.png"))



##### Redo with all the cells across all times #####
##### Obtain and manage data #####
seurat <- readRDS(paste0(datadir,"seurat_integrated_all_times_clustered.rds"))


### Test fitdist function for looking at just normal on scale.data and negative binomial and zero-inflated negative binomial distributions ###
number <- 16437
# norm_temp <- fitDist(as.matrix(seurat[["SCT"]]@scale.data[number,]), k = 2, type = "realAll")
fit.norm <- fitdist(seurat[["SCT"]]@scale.data[number,], "norm")
fit.nbinom <- fitdist(seurat[["SCT"]]@counts[number,], "nbinom")
fit.zinbinom <- fitdist(seurat[["SCT"]]@counts[number,], "ZINBI", start = list(mu = 1, sigma = 1))

fit.norm$aic
fit.nbinom$aic
fit.zinbinom$aic



ref_dist <- rnorm(length(seurat[["SCT"]]@scale.data[number,]), mean = as.numeric(fit.norm$estimate["mean"]), sd = as.numeric(fit.norm$estimate["sd"]))
plot_scaled <- ggplot(data.frame("Counts" = seurat[["SCT"]]@scale.data[number,]), aes(Counts)) +
					geom_histogram(binwidth = 0.25) +
					theme_classic() +
					geom_density(data = data.frame("Counts" = ref_dist), aes(y=0.25 * ..count..))
ggsave(plot_scaled, filename = paste0(outdir,"all_test_dist_histo_scale_data.png"))

qq_ref_dist <- ggplot(data.frame("Empirical" = seurat[["SCT"]]@scale.data[number,][order(seurat[["SCT"]]@scale.data[number,])], "Simulated" = ref_dist[order(ref_dist)]), aes(Empirical, Simulated)) +
					geom_point() +
					geom_abline(intercept = 0, slope = 1) +
					theme_classic()
ggsave(qq_ref_dist, filename = paste0(outdir,"all_empirical_simulated_qq_scale_data.png"))



ref_nbinom_dist <- rnbinom(length(seurat[["SCT"]]@counts[number,]), size = fit.nbinom$estimate["size"], mu = fit.nbinom$estimate["mu"])
plot_scaled <- ggplot(data.frame("Counts" = seurat[["SCT"]]@counts[number,]), aes(Counts)) +
					geom_histogram(binwidth = 1) +
					theme_classic() +
					geom_density(data = data.frame("Counts" = ref_nbinom_dist), stat = "count")
ggsave(plot_scaled, filename = paste0(outdir,"all_test_dist_histo_nbinom.png"))

qq_ref_nbinom_dist <- ggplot(data.frame("Empirical" = seurat[["SCT"]]@counts[number,][order(seurat[["SCT"]]@counts[number,])], "Simulated" = ref_nbinom_dist[order(ref_nbinom_dist)]), aes(Empirical, Simulated)) +
					geom_point(alpha = 0.1) +
					geom_abline(intercept = 0, slope = 1) +
					theme_classic()
ggsave(qq_ref_nbinom_dist, filename = paste0(outdir,"all_empirical_simulated_qq_nbinom.png"))



ref_zinbinom_dist <- rZINBI(length(seurat[["SCT"]]@counts[number,]), sigma = as.numeric(fit.zinbinom$estimate["sigma"]), mu = as.numeric(fit.zinbinom$estimate["mu"]))
plot_scaled <- ggplot(data.frame("Counts" = seurat[["SCT"]]@counts[number,]), aes(Counts)) +
					geom_histogram(binwidth = 1) +
					theme_classic() +
					geom_density(data = data.frame("Counts" = ref_zinbinom_dist), stat = "count")
ggsave(plot_scaled, filename = paste0(outdir,"all_test_dist_histo_zinbinom.png"))

qq_ref_zinbinom_dist <- ggplot(data.frame("Empirical" = seurat[["SCT"]]@counts[number,][order(seurat[["SCT"]]@counts[number,])], "Simulated" = ref_zinbinom_dist[order(ref_zinbinom_dist)]), aes(Empirical, Simulated)) +
					geom_point(alpha = 0.1) +
					geom_abline(intercept = 0, slope = 1) +
					theme_classic()
ggsave(qq_ref_zinbinom_dist, filename = paste0(outdir,"all_empirical_simulated_qq_zinbinom.png"))




### Test fitdist function for looking at just normal on scale.data and negative binomial and zero-inflated negative binomial distributions ###
seurat <- readRDS(paste0(datadir,"seurat_integrated_all_times_clustered.rds"))


## Keep only genes expressed in at least 1% of cells 
seurat_sub <- subset(seurat, features = rownames(seurat)[which(rowSums(seurat[["SCT"]]@counts > 0)/ncol(seurat[["SCT"]]@counts) >= 0.01)])
saveRDS(seurat_sub, paste0(outdir,"seurat_integrated_all_times_clustered_1pct_expressing.rds"))
seurat_sub <- readRDS(paste0(outdir,"seurat_integrated_all_times_clustered_1pct_expressing.rds"))


### Make a list of the genes to be used by snakemake in variance_partition_post_reviedw.smk ###
genes_dt <- data.table(Gene = rownames(seurat_sub))
fwrite(genes_dt, paste0(outdir,"seurat_integrated_all_times_clustered_1pct_expressing_genelist.tsv"), sep = "\t")


aic_df <- as.data.frame(matrix(nrow = 0, ncol = 3))
colnames(aic_df) <- c("gene", "normal", "negative_binomial")

fit <- list()

for (number in c(1:9,9733)){
	fit[["norm"]] <- fitdist(seurat_sub[["SCT"]]@scale.data[number,], "norm")
	fit[["nbinom"]] <- fitdist(seurat_sub[["SCT"]]@counts[number,], "nbinom")
	# fit[["zinbinom"]] <- tryCatch({
	# 		fit.zinbinom <- fitdist(seurat_sub[["SCT"]]@counts[number,], "ZINBI", start = list(mu = 1, sigma = 1))
	# 	}, error = function(e) {
	# 		print(e)
	# 		fit.zinbinom <- NA
	# 		fit.zinbinom$aic <- NA
	# 		return(fit.zinbinom)
	# 	})
	# fit.zinbinom <- fitdist(seurat_sub[["SCT"]]@counts[number,], "ZINBI", start = list(mu = ifelse((round(mean(seurat_sub[["SCT"]]@counts[number,][which(seurat_sub[["SCT"]]@counts[number,] > 0)])) > 1), round(mean(seurat_sub[["SCT"]]@counts[number,][which(seurat_sub[["SCT"]]@counts[number,] > 0)])), 2), sigma = 1))
	# fit.zinbinom <- fitdist(seurat_sub[["SCT"]]@counts[number,], "ZINBI", start = list(mu = 1, sigma = 1))
	aic_df <- rbind(aic_df, data.frame("gene" = rownames(seurat_sub[["SCT"]])[number], "normal" = fit[["norm"]]$aic, "negative_binomial" = fit[["nbinom"]]$aic))
}


aic_df$difference <- aic_df$normal - aic_df$negative_binomial
aic_df$better_model <- ifelse(aic_df$normal < aic_df$negative_binomial, "normal","negative_binomial")



