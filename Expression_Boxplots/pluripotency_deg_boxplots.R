library(data.table)
library(tidyverse)
library(Seurat)




outdir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Expression_Boxplots/pluri_degs/"
dir.create(outdir, recursive = TRUE)


pluri_genes <- fread("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/pluripotency_genes.tsv")

seurat <- readRDS("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/All_data_integrated_remove_bad/seurat_integrated_all_times_clustered.rds")
seurat@meta.data$Location <- gsub("_.+", "", seurat@meta.data$Location)
seurat@meta.data$Cryopreserved <- ifelse(grepl("Thawed", seurat@meta.data$Location_Time), "Cryopreserved", "Fresh")
seurat@meta.data$Location_Cryopreserved_Line <- paste0(seurat@meta.data$Location, "_", seurat@meta.data$Cryopreserved, "_", seurat@meta.data$Final_Assignment)
seurat@meta.data$Village <- gsub("Baseline", "Uni-Culture", seurat@meta.data$Time) %>% 
        gsub("Thawed Village Day 0", "Uni-Culture", .) %>%
        gsub("Thawed Village Day 7", "Village", .) %>%
        gsub("Village Day 4", "Village", .)


##### Combine replicates from same groups together #####
seurat_list <- lapply(unique(seurat@meta.data$Location_Cryopreserved_Line), function(x){
	subset(seurat, subset = Location_Cryopreserved_Line == x)
})
names(seurat_list) <- unique(seurat@meta.data$Location_Cryopreserved_Line)



seurat_list <- lapply(seurat_list, function(x){
    Idents(x) <- "Village"
    return(x)
})



seurat_list <- lapply(seurat_list, function(x){
    SCTransform(x, vars.to.regress = c("scores.G1", "scores.S", "scores.G2M", "percent.mt", "percent.rb"), return.only.var.genes = FALSE)
})

seurat_list <- lapply(seurat_list, function(x){
    PrepSCTFindMarkers(x)
})

saveRDS(seurat_list, paste0(outdir, "sct_location_cryo_line.rds"))
seurat_list <- readRDS(paste0(outdir, "sct_location_cryo_line.rds"))



##### Combine Village and Baselines together #####
### Try with logistic MAST ###
MAST_DEGs <- lapply(seurat_list, function(x){
    FindMarkers(x, ident.1 = "Uni-Culture", ident.2 = "Village", latent.vars = "MULTI_classification", test.use = "MAST", logfc.threshold = 0, assay = "RNA")
})

saveRDS(MAST_DEGs, paste0(outdir, "MAST_DEGs.rds"))
MAST_DEGs <- readRDS(paste0(outdir, "MAST_DEGs.rds"))


### Update multiple testing for all sites ###
MAST_DEGs_nrow_list <- lapply(MAST_DEGs, function(x) data.table(nrow(x)))
MAST_DEGs_nrow <- sum(do.call(rbind, MAST_DEGs_nrow_list)$V1)

MAST_DEGs <- lapply(MAST_DEGs, function(x){
    x$p_val_adj_updated <- p.adjust(x$p_val, method = "bonferroni", n = MAST_DEGs_nrow)
    return(x)
})


MAST_DEGs_pluri <- lapply(MAST_DEGs, function(x){
    x$ENSG <- rownames(x)
    x <- data.table(x)
    x <- x[pluri_genes, on = c("ENSG")]
    x <- x[!is.na(p_val) & p_val_adj_updated < 0.05]
    return(x)
})

MAST_DEGs_pluri <- lapply(MAST_DEGs, function(x){
    x$ENSG <- rownames(x)
    x <- data.table(x)
    x <- x[pluri_genes, on = c("ENSG")]
    x <- x[!is.na(p_val) & p_val_adj < 0.05]
    return(x)
})

MAST_DEGs_pluri_4 <- lapply(MAST_DEGs_pluri, function(x){
    x[GeneID %in% c("MYC", "NANOG", "POU5F1", "SOX2")]
})


### Try with logistic regression ###
LR_DEGs <- lapply(seurat_list, function(x){
    FindMarkers(x, ident.1 = "Uni-Culture", ident.2 = "Village", latent.vars = "MULTI_classification", test.use = "LR", logfc.threshold = 0)
})

saveRDS(LR_DEGs, paste0(outdir, "LR_DEGs.rds"))
##### Will use LR instead of MAST - MAST finds more DEG and lower p-values which I'm not sure is legitimate #####


### Update multiple testing for all sites ###
LR_DEGs_nrow_list <- lapply(LR_DEGs, function(x) data.table(nrow(x)))
LR_DEGs_nrow <- sum(do.call(rbind, LR_DEGs_nrow_list)$V1)

LR_DEGs <- lapply(LR_DEGs, function(x){
    x$p_val_adj_updated <- p.adjust(x$p_val, method = "bonferroni", n = LR_DEGs_nrow)
    return(x)
})


LR_DEGs_pluri <- lapply(LR_DEGs, function(x){
    x$ENSG <- rownames(x)
    x <- data.table(x)
    x <- x[pluri_genes, on = c("ENSG")]
    x <- x[!is.na(p_val) & p_val_adj_updated < 0.05]
    return(x)
})


LR_DEGs_pluri_4 <- lapply(LR_DEGs_pluri, function(x){
    x[GeneID %in% c("MYC", "NANOG", "POU5F1", "SOX2")]
})



saveRDS(LR_DEGs_pluri_4, paste0(outdir, "LR_DEGs_4pluri_genes.rds"))

