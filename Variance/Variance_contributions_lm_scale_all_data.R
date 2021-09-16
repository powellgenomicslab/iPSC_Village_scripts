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
outdir <- paste0(dir,"output/Variance_contributions_lm_scale_all_data/")


dir.create(outdir)


location_SCT <- readRDS(paste0(dir,"output/Variance_contributions_lm_scale_data/seurat_SCT_wRB_sydney_regress_all_covs.rds"))
cryo_SCT <- readRDS(paste0(dir,"output/Variance_contributions_lm_scale_data_Sydney/seurat_SCT_wRB_sydney_regress_all_covs.rds"))


### Update names to combine objects ###
names(location_SCT) <- c("Baseline", "Village")
location_SCT <- lapply(location_SCT, function(x){
	names(x) <- c("Brisbane", "Melbourne", "Sydney_Fresh")
	return(x)
})


cryo_SCT_updated <- list()
cryo_SCT_updated[["Baseline"]][["Sydney_Fresh"]] <- cryo_SCT[["Baseline"]]
cryo_SCT_updated[["Baseline"]][["Sydney_Cryopreserved"]] <- cryo_SCT[["Thawed Village Day 0"]]
cryo_SCT_updated[["Village"]][["Sydney_Fresh"]] <- cryo_SCT[["Village Day 4"]]
cryo_SCT_updated[["Village"]][["Sydney_Cryopreserved"]] <- cryo_SCT[["Thawed Village Day 7"]]


SCT_combined <- location_SCT
SCT_combined[["Baseline"]][["Sydney_Cryopreserved"]] <- cryo_SCT_updated[["Baseline"]][["Sydney_Cryopreserved"]]
SCT_combined[["Village"]][["Sydney_Cryopreserved"]] <- cryo_SCT_updated[["Village"]][["Sydney_Cryopreserved"]]


#### Filter out Genes expressed in low percentage of cells #####
pct_expressing <- lapply(SCT_combined, function(x){
	lapply(x, function(y){
		temp <- data.frame(rowSums(as.matrix(y[["SCT"]]@counts) > 0)/ncol(y[["SCT"]]@counts))
		colnames(temp) <- "Cell_proportion"
		return(temp)
	})
})


### 10% to start
SCT_combined <- lapply(names(SCT_combined), function(x){
	temp <- lapply(names(SCT_combined[[x]]), function(y){
		subset(SCT_combined[[x]][[y]], features = rownames(pct_expressing[[x]][[y]])[which(pct_expressing[[x]][[y]]$Cell_proportion >= 0.1)])
	})
	names(temp) <- names(SCT_combined[[x]])
	return(temp)
})
names(SCT_combined) <- names(cryo_SCT_updated)


sce_separate_list <- lapply(SCT_combined, function(x){
	lapply(x, function(y){
		SingleCellExperiment(assays = list(counts = y[["SCT"]]@counts, logcounts = y[["SCT"]]@data, normcounts = y[["SCT"]]@scale.data), colData = y@meta.data)
	})
})



if (!file.exists(paste0(outdir, "heirarchical_partitioning_results_sep_all.rds"))){
	scran <- list()
	i <- 1
	for(x in names(SCT_combined)){
		for (y in names(SCT_combined[[x]])){
			print(i)
			df <- data.frame("Line" = SCT_combined[[x]][[y]]@meta.data$Final_Assignment)
			scran[[x]][[y]] <- getVarianceExplained(sce_separate_list[[x]][[y]], variables = df, exprs_values = "normcounts")
			i <- 1 + i
		}
	}
	saveRDS(scran, paste0(outdir, "heirarchical_partitioning_results_sep_all.rds"))
} else {
	scran <- readRDS(paste0(outdir, "heirarchical_partitioning_results_sep_all.rds"))
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
	geom_text(x = 0, y = 75, aes(label = paste0("R^2 = ", Rsquared)), data = Rsquared, hjust = 0, color = variable_colors["Line"])
save_figs(pCorr_Location_Point_Line, paste0(outdir, "Correlation_Point_Location_Line_cov"))


