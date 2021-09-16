library(Seurat)
library(tidyverse)
library(lme4)

save_figs <- function(plot, basename, width = 17, height = 17, units = "cm"){
    ggsave(plot, filename = paste0(basename,".png"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".pdf"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".eps"), height = height, width = width, units = units)
}

##### Set Up Directories #####
dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
outdir <- paste0(dir,"output/Variance_contributions/")
datadir <- paste0(dir,"output/All_data_integrated_remove_bad/")

dir.create(outdir)


##### Read in Data #####
seurat <- readRDS(paste0(datadir,"seurat_integrated_all_times_clustered.rds"))


##### Remove low expressing genes and look at expression distribution #####
## Separate by cluster and cell line
seurat_list <- list()
for (cluster in unique(seurat@meta.data$integrated_snn_res.0.28)){
	seurat_list[[cluster]] <- list()
	for (line in unique(seurat@meta.data$Final_Assignment)){
		seurat_list[[cluster]][[line]] <- list()
		for (experiment in unique(seurat@meta.data$Location_Time)){
			seurat_list[[cluster]][[line]][[experiment]] <- subset(seurat, subset = (Final_Assignment == line & integrated_snn_res.0.28 == cluster & Location_Time == experiment))
		}
	}
}


## Count UMIs per cell line per cluster per experimental condition ##
row_sums_list <- lapply(names(seurat_list), function(x){
	temp2 <- lapply(names(seurat_list[[x]]), function(y){
		temp <- data.frame(rowSums(seurat_list[[x]][[y]][["SCT"]]@counts))
		colnames(temp) <- "Corrected_UMIs"
		temp$CellLine <- y
		temp$Cluster <- x
		temp$Gene <- rownames(temp)
		temp$Proportion <- apply(seurat_list[[x]][[y]][["SCT"]]@counts, 1, function(z){
			length(which(z > 0))/length(z)
		})
		temp$N <- ncol(seurat_list[[x]][[y]][["SCT"]]@counts)
		return(temp)
	})
	names(temp2) <- names(seurat_list[[x]])
	return(temp2)
})

names(row_sums_list) <- names(seurat_list)

row_sums_df_list <- lapply(row_sums_list, function(x){
	do.call(rbind, x)
})

row_sums_df <- do.call(rbind, row_sums_df_list)

row_sums_df_wide <- pivot_wider(row_sums_df, values_from = c("Corrected_UMIs", "Proportion", "N"), names_from = "CellLine")
row_sums_df_wide$Total_Counts <- rowSums(row_sums_df_wide[,c("FSA0006", "TOB0421", "MBE1006")])


## Remove genes with 0 counts in any line ##
row_sums_df_wide <- row_sums_df_wide[which(row_sums_df_wide$Corrected_UMIs_FSA0006 > 0 & row_sums_df_wide$Corrected_UMIs_TOB0421 > 0 & row_sums_df_wide$Corrected_UMIs_MBE1006 > 0),]
row_sums_df_wide[which(row_sums_df_wide$Corrected_UMIs_FSA0006 > 0 & row_sums_df_wide$Corrected_UMIs_TOB0421 > 0 & row_sums_df_wide$Corrected_UMIs_MBE1006 > 0),]


pHisto_Facet_Cluster <- ggplot(row_sums_df_wide, aes(Total_Counts)) +
	geom_histogram() +
	facet_wrap(Cluster ~ ., scales = "free") +
	theme_classic()  +
	# scale_y_continuous(trans = "log10") +
	scale_x_continuous(trans = "log10")

save_figs(pHisto_Facet_Cluster, basename = paste0(outdir, "UMI_histogram_cluster_faceted"))


##### Make a dataframe for expression  per cluster #####
seurat_list_subset <- lapply(names(seurat_list), function(x){
	temp <- lapply(names(seurat_list[[x]]), function(y){
		subset(seurat_list[[x]][[y]], features = row_sums_df_wide[which(row_sums_df_wide$Cluster == x), ]$Gene)
	})
	names(temp) <- names(seurat_list[[x]])
	return(temp)
})
names(seurat_list_subset) <- names(seurat_list)

##### Add new data from Maciej and Alice into meta.data object #####
seurat_list_subset <- lapply(seurat_list_subset, function(x){
	lapply(x, function(y){
		y$Age <- ifelse(y@meta.data$Final_Assignment == "MBE1006", 58, ifelse(y@meta.data$Final_Assignment == "TOB0421", 83.3, ifelse(y@meta.data$Final_Assignment == "FSA0006", 51.7, NA)))
		y$Sex <- ifelse(y@meta.data$Final_Assignment == "MBE1006", "Male", ifelse(y@meta.data$Final_Assignment == "TOB0421", "Female", ifelse(y@meta.data$Final_Assignment == "FSA0006", "Male", NA)))
		return(y)
	})
})


seurat_list_cluster <- lapply(seurat_list_subset, function(x){
	merge(x[[1]], c(x[[2]], x[[3]]))
})

temp <- list()
seurat_list_cluster_sample <- lapply(names(seurat_list_cluster), function(x){
	temp[[x]] <- list()
	temp[[x]] <- lapply(unique(seurat_list_cluster[[x]]@meta.data$Location_Time), function(y){
		subset(seurat_list_cluster[[x]], subset = Location_Time == y)
	})
	names(temp[[x]]) <- unique(seurat_list_cluster[[x]]@meta.data$Location_Time)
	return(temp[[x]])
})
names(seurat_list_cluster_sample) <- names(seurat_list_cluster)


##### Run linear model for each gene in each cluster with covariates #####
data_frame_list <- lapply(seurat_list_cluster_sample, function(x){
	lapply(x, function(y){
		data.frame("Barcode" = colnames(y), "Individual" = y@meta.data$Final_Assignment, "Age" = y@meta.data$Age, "Sex" = y@meta.data$Sex, "Cycle_Phase" = y@meta.data$phases, "MT_Percent" = y@meta.data$percent.mt)
	})
})


lapply(data_frame_list, function(x) lapply(x, function(y) table(y$Individual)))



lapply(names(seurat_list_cluster_sample), function(x){
	lapply(names(seurat_list_cluster_sample[[x]]), function(y){
		for (gene in rownames(y)){
			data_frame_list[[x]][[y]]$Expression <- x[["SCT"]]@counts[gene,]
			lmer(Expression ~ (1 | individual, data = data_frame_list[[x]][[y]]
		}
	})
})



