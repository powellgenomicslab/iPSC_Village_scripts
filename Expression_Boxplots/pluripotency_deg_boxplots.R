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
    n = min(table(Idents(x)))
    print(n)
    # FindMarkers(x, ident.1 = "Uni-Culture", ident.2 = "Village", latent.vars = "MULTI_classification", test.use = "MAST", logfc.threshold = 0, assay = "RNA", max.cells.per.ident = n)
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

# MAST_DEGs_pluri <- lapply(MAST_DEGs, function(x){
#     x$ENSG <- rownames(x)
#     x <- data.table(x)
#     x <- x[pluri_genes, on = c("ENSG")]
#     x <- x[!is.na(p_val) & p_val_adj < 0.05]
#     return(x)
# })

MAST_DEGs_pluri_4 <- lapply(MAST_DEGs_pluri, function(x){
    x[GeneID %in% c("MYC", "NANOG", "POU5F1", "SOX2")]
})


### Try with logistic regression ###
LR_DEGs <- lapply(seurat_list, function(x){
    n = min(table(Idents(x)))
    print(n)
    # FindMarkers(x, ident.1 = "Uni-Culture", ident.2 = "Village", latent.vars = "MULTI_classification", test.use = "LR", logfc.threshold = 0, max.cells.per.ident = n, assay = "RNA")
})

saveRDS(LR_DEGs, paste0(outdir, "LR_DEGs.rds"))
LR_DEGs <- readRDS(paste0(outdir, "LR_DEGs.rds"))
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



### Try with logistic regression  downsampled to lowest number across all pools###
LR_DEGs_sub <- lapply(seurat_list, function(x){
    FindMarkers(x, ident.1 = "Uni-Culture", ident.2 = "Village", latent.vars = "MULTI_classification", test.use = "LR", logfc.threshold = 0, max.cells.per.ident = 572, assay = "RNA")
})

saveRDS(LR_DEGs_sub, paste0(outdir, "LR_DEGs_subset.rds"))
LR_DEGs_sub <- readRDS(paste0(outdir, "LR_DEGs_subset.rds"))


### Update multiple testing for all sites ###
LR_DEGs_sub_nrow_list <- lapply(LR_DEGs_sub, function(x) data.table(nrow(x)))
LR_DEGs_sub_nrow <- sum(do.call(rbind, LR_DEGs_sub_nrow_list)$V1)

LR_DEGs_sub <- lapply(LR_DEGs_sub, function(x){
    x$p_val_adj_updated <- p.adjust(x$p_val, method = "bonferroni", n = LR_DEGs_sub_nrow)
    return(x)
})


LR_DEGs_sub_pluri <- lapply(LR_DEGs_sub, function(x){
    x$ENSG <- rownames(x)
    x <- data.table(x)
    x <- x[pluri_genes, on = c("ENSG")]
    return(x)
})


LR_DEGs_sub_pluri_sig <- lapply(LR_DEGs_sub_pluri, function(x){
    x <- x[!is.na(p_val) & p_val_adj_updated < 0.05]
    return(x)
})



LR_DEGs_sub_pluri_sig_4 <- lapply(LR_DEGs_sub_pluri_sig, function(x){
    x[GeneID %in% c("MYC", "NANOG", "POU5F1", "SOX2")]
})


saveRDS(LR_DEGs_sub_pluri_sig_4, paste0(outdir, "LR_DEGs_subset_4pluri_genes.rds"))



##### Make volcano of significants #####
LR_DEGs_sub_pluri <- lapply(names(LR_DEGs_sub_pluri), function(x){
    LR_DEGs_sub_pluri[[x]]$Group <- x
    return(LR_DEGs_sub_pluri[[x]])
})


LR_DEGs_sub_pluri_dt <- do.call(rbind, LR_DEGs_sub_pluri)
LR_DEGs_sub_pluri_dt$p_val_adj_updated <- ifelse(is.na(LR_DEGs_sub_pluri_dt$p_val_adj_updated), 1, LR_DEGs_sub_pluri_dt$p_val_adj_updated)
LR_DEGs_sub_pluri_dt$significant <- ifelse(LR_DEGs_sub_pluri_dt$p_val_adj_updated > 0.05, "not significant", ifelse(abs(LR_DEGs_sub_pluri_dt$avg_log2FC) > 1,"significant large effect size", "significant small effect size"))
LR_DEGs_sub_pluri_dt$significant  <- factor(LR_DEGs_sub_pluri_dt$significant, levels = c("significant large effect size", "significant small effect size", "not significant"))



pVolcano <- ggplot(LR_DEGs_sub_pluri_dt[!grepl("Cryopreserved", Group)], aes(avg_log2FC, -log2(p_val_adj_updated), color = significant)) +
    geom_point() +
    theme_classic() +
    geom_vline(xintercept = -1, linetype = "dashed") +
    geom_vline(xintercept = 1, linetype = "dashed") +
    geom_hline(yintercept = 0.05, linetype = "dashed") +
    scale_color_manual(values = c("black", "grey50", "grey80")) +
    ylab("-log2(P Value)") +
    xlab("log2(Fold Change)") +
    ggtitle("Pluripotent Gene Differential Expression\nBetween Uni-Culture and Village") +
    theme(plot.title = element_text(hjust = 0.5))


ggsave(pVolcano, filename = paste0(outdir, "volcano.png"), width = 4.5, height = 3)



pVolcano_cryo <- ggplot(LR_DEGs_sub_pluri_dt[grepl("Sydney", Group)], aes(avg_log2FC, -log2(p_val_adj_updated), color = significant)) +
    geom_point() +
    theme_classic() +
    geom_vline(xintercept = -1, linetype = "dashed") +
    geom_vline(xintercept = 1, linetype = "dashed") +
    geom_hline(yintercept = 0.05, linetype = "dashed") +
    scale_color_manual(values = c("black", "grey50", "grey80")) +
    ylab("-log2(P Value)") +
    xlab("log2(Fold Change)") +
    ggtitle("Pluripotent Gene Differential Expression\nBetween Uni-Culture and Village") +
    theme(plot.title = element_text(hjust = 0.5))


ggsave(pVolcano_cryo, filename = paste0(outdir, "volcano_cryo.png"), width = 4.5, height = 3)


LR_DEGs_sub_pluri_dt$Group <- gsub("Brisbane", "Site 1", LR_DEGs_sub_pluri_dt$Group) %>%
                                gsub("Melbourne", "Site 2", .) %>%
                                gsub("Sydney", "Site 3", .)

fwrite(LR_DEGs_sub_pluri_dt, paste0(outdir, "DEGs4volcano.tsv"), sep = "\t")