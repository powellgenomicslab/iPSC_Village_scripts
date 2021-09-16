library(haven)
library(ggplot2)
library(lme4)
library(glmmTMB)
library(Seurat)



# ##### Bring in variables #####
# ### Bring in arguments
# args <- commandArgs(trailingOnly = TRUE)
# outdir <- paste0(args[1])
# number <- as.numeric(args[2])




dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
outdir <- paste0(dir,"output/nb_partitioning_village_separate/data/")
dir.create(outdir)


##### Read in Data ######
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


SCT_combined <- lapply(SCT_combined, function(x){
	lapply(x, function(y){
		subset(y, features = rownames(y)[which(rowSums(y[["SCT"]]@counts > 0)/ncol(y[["SCT"]]@counts) >= 0.01)])
	})
})


lapply(names(SCT_combined), function(x){
	lapply(names(SCT_combined[[x]]), function(y){
		saveRDS(SCT_combined[[x]][[y]], paste0(outdir, x, "_", y,"_seurat.rds"))
	})
})


SCT_combined <- lapply(SCT_combined, function(x){
	lapply(x, function(y){
		subset(y, features = rownames(y)[which(rowSums(y[["SCT"]]@counts > 0)/ncol(y[["SCT"]]@counts) >= 0.1)])
	})
})
