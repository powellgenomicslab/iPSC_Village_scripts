##### Read in Arguments #####
print("Reading and assigning input arguments")
args <- commandArgs(TRUE)
arguments <- read.table(args, header = F)
stage <- arguments[1,]
outdir <- arguments[2,]
dir <- arguments[3,]
sge <- arguments[4,]
print(sge)
sge <- as.numeric(as.character(sge))
print(sge)
print(outdir)

if ( stage[1]=="CellCycle" ){
    library(scran, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/generalR/lib/R/library")
    library(dplyr, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/generalR/lib/R/library")
    library(tidyr, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/generalR/lib/R/library")
    library(tidyverse, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/generalR/lib/R/library")
    library(ggplot2, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/generalR/lib/R/library")
    library(Seurat, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/harmony/lib/R/library")
    library(jcolors, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/harmony/lib/R/library")
    library(cowplot, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/harmony/lib/R/library")
    library(RColorBrewer, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/harmony/lib/R/library")
    library(readr, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/harmony/lib/R/library")
    library(purrr, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/harmony/lib/R/library")
} else {
    library(dplyr, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/harmony/lib/R/library")
    library(tidyr, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/harmony/lib/R/library")
    library(tidyverse, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/harmony/lib/R/library")
    library(ggplot2, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/harmony/lib/R/library")
    library(clustree, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/harmony/lib/R/library")
    library(reticulate, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/harmony/lib/R/library")
    library(Seurat, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/harmony/lib/R/library")
    library(schex, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/harmony/lib/R/library")
    library(BiocParallel, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/harmony/lib/R/library")
    library(coop, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/harmony/lib/R/library")
    library(harmony, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/harmony/lib/R/library")
    library(cowplot, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/harmony/lib/R/library")
    library(RColorBrewer, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/harmony/lib/R/library")
    library(readr, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/harmony/lib/R/library")
    library(purrr, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/harmony/lib/R/library")
    library(jcolors, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/harmony/lib/R/library")
}
print(sessionInfo())

cell_line_colors <- c(jcolors("default"),"sienna")
names(cell_line_colors) <- c("FSA006","MBE1006","doublet","TOB421","unassigned","combination_fail_remove")

mad_function <- function(seurat, column, number_mad){
    mad <- mad(seurat@meta.data[,column])
    low <- median(seurat@meta.data[,column]) - number_mad*mad
    high <- median(seurat@meta.data[,column]) + number_mad*mad
    print("The mad is:")
    print(mad)
    print("The lower bound is:")
    print(low)
    print("The upper bound is:")
    print(high)
    seurat@meta.data[,paste0(column,"_mad")] <- ifelse((seurat@meta.data[,column] > low & seurat@meta.data[,column] < high),"NotOutlier", "Outlier")
    return(seurat)
}

PoolColors <- brewer.pal(n = 6, name = "Dark2")




##### First tried identifying the cells and doublets with all cells in one Seurat object
    ### THIS WAS NOT SUCCESSFUL - always identify sample and doublets from hashing before merging the objects together (as done below)
    # pools <- dir("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/", pattern = "DRENEA")
    # dirs_hash <- paste0("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Barcodes_200213_A00152_0206_AHK3HYDRXX/Hashing/", pools, "/umi_count/")
    # dirs10x <- paste0("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/", pools, "/outs/filtered_feature_bc_matrix/")


    # ## Read in expression data
    # counts_list <- lapply(dirs10x, function(x){
    #     Read10X(x)
    # })
    # names(counts_list) <- pools

    # ## Add poolnames to cell names so can easily match back if there are any cells with same barcoe
    # counts_list <- lapply(names(counts_list), function(x){
    #     colnames(counts_list[[x]]) <- paste0(x, "_", colnames(counts_list[[x]]))
    #     return(counts_list[[x]])
    # })
    # names(counts_list) <- pools


    # ### Read in the hashing data ###
    # hash_list <- lapply(dirs_hash, function(x){
    #     Read10X(x, gene.column =1)
    # })
    # names(hash_list) <- pools
    # print(lapply(hash_list, rownames))

    # hash_list <- lapply(hash_list, function(x){
    #     x <- as.data.frame(t(as.data.frame(x)))
    #     x$Barcode <- rownames(x)
    #     return(x)
    # })
    # names(hash_list) <- pools
    # print(lapply(hash_list, colnames))


    # hash_list <- lapply(names(hash_list), function(x){
    #     tmp <- as.data.frame(t(hash_list[[x]]))
    #     colnames(tmp) <- hash_list[[x]]$Barcode
    #     tmp <- tmp[1:(nrow(tmp)-1),]
    # })
    # names(hash_list) <- pools
    # print(lapply(hash_list, rownames))

    # ### Add pool names to barcodes so all are unique ###
    # hash_list <- lapply(names(hash_list), function(x){
    #     colnames(hash_list[[x]]) <- paste0(x, "_", colnames(hash_list[[x]]))
    #     return(hash_list[[x]])
    # })
    # names(hash_list) <- pools
    # print(lapply(hash_list, rownames))


    # hash_list <- lapply(hash_list, function(x){
    #     x$Hashtag <- rownames(x)
    #     return(x)
    # })
    # names(hash_list) <- pools
    # print(lapply(hash_list, rownames))

    # ## Combine the matrices together
    # counts <- do.call(cbind, counts_list)

    # ### Get the metadata ###
    # meta <- read_delim(paste0(dir,"/output/CompareDemultiplexing/doublet_metadata.txt"), delim = "\t")

    # ## Join with the barcode list ##
    # meta <- as.data.frame(meta)
    # rownames(meta) <- meta$Barcode
    # print(dim(meta))

    # ### Combine hashing tables and correct the direction
    # hash <- hash_list[["DRENEA_1"]]
    # for (name in names(hash_list)[2:length(names(hash_list))]){
    #     hash <- full_join(hash, hash_list[[name]], by = "Hashtag")
    # }
    # print(dim(hash))
    # hash_barcode <- colnames(hash)
    # hash_t_colnames <- hash$Hashtag
    # hash_t <- as.data.frame(t(hash))
    # colnames(hash_t) <- hash_t_colnames
    # hash_t <- hash_t[-c(which(rownames(hash_t) == "Hashtag")),]
    # hash_t <- as.data.frame(data.matrix(hash_t))
    # hash_t[is.na(hash_t)] <- 0
    # hash_t$Barcode <- colnames(hash)[-c(which(colnames(hash) == "Hashtag"))]
    # print(dim(hash_t))
    # rownames(hash_t) <- hash_t$Barcode


    # barcodes <- as.data.frame(hash_t$Barcode)
    # colnames(barcodes) <- "Barcode"
    # print(dim(counts))
    # counts <- counts[,barcodes$Barcode]
    # print(dim(counts))
    # hash_t$Barcode <- NULL

    # ## Make seurat object
    # seurat <- CreateSeuratObject(counts = counts, meta.data = meta)

    # ## Add QC metrics
    # RbGeneList <- read.delim(file = "/directflow/SCCGGroupShare/projects/DrewNeavin/References/RibosomalGeneList.txt",header = F)
    # print("Calculating Mt %")
    # seurat <- PercentageFeatureSet(seurat, pattern = "^MT-", col.name = "percent.mt")
    # print("Calculating Rb %")
    # seurat <- PercentageFeatureSet(seurat, features = RbGeneList$V1, col.name = "percent.rb")

    # hash_t_t <- as.data.frame(t(hash_t))
    # write.table(hash_t_t, paste0(outdir,"hashtag_dataframe.txt"), sep = "\t")
    # # hash_t_t <- read.table(paste0(outdir,"hashtag_dataframe.txt"), sep = "\t")

    # ### Normalize expression data ###
    # seurat <- SCTransform(seurat, verbose = TRUE)

    # # # Normalize RNA data with log normalization
    # # seurat <- NormalizeData(seurat)
    # # # Find and scale variable features
    # # seurat <- FindVariableFeatures(seurat, selection.method = "mean.var.plot")
    # # seurat <- ScaleData(seurat, features = VariableFeatures(seurat))
    # saveRDS(seurat, paste0(outdir,"seurat_object.rds"))

    # ### Add hashtag data ###
    # seurat[["HTO"]] <- CreateAssayObject(counts = hash_t_t)
    # seurat <- NormalizeData(seurat, assay = "HTO", normalization.method = "CLR")

    # seurat <- MULTIseqDemux(seurat, assay = "HTO", autoThresh = TRUE)

    # ### Make a table of the doublet results ###
    # print(table(seurat$MULTI_ID))

    # ### Visualize results ###
    # Idents(seurat) <- "HTO_maxID"
    # p_hto_ridge_all <- RidgePlot(seurat, assay = "HTO", features = rownames(seurat[["HTO"]]), ncol = 2)
    # p_hto_ridge_sng <- RidgePlot(subset(seurat, subset = MULTI_ID == "Singlet" ), assay = "HTO", features = rownames(seurat[["HTO"]]), ncol = 2)
    # p_hto_ridge_neg <- RidgePlot(subset(seurat, subset = MULTI_ID == "Negative" ), assay = "HTO", features = rownames(seurat[["HTO"]]), ncol = 2)
    # p_hto_ridge_dbl <- RidgePlot(subset(seurat, subset = MULTI_ID == "Doublet" ), assay = "HTO", features = rownames(seurat[["HTO"]]), ncol = 2)

    # ggsave(p_hto_ridge_all, filename = paste0(outdir,"RidgePlot_AllCells.png"))
    # ggsave(p_hto_ridge_sng, filename = paste0(outdir,"RidgePlot_htoSinglet.png"))
    # ggsave(p_hto_ridge_neg, filename = paste0(outdir,"RidgePlot_htoNegative.png"))
    # ggsave(p_hto_ridge_dbl, filename = paste0(outdir,"RidgePlot_htoDoublet.png"))



    # saveRDS(seurat, paste0(outdir,"AllCellsSeurat_meta.rds"))
    # # seurat <- readRDS(paste0(outdir,"AllCellsSeurat_meta.rds"))

    # ## Scrublet vs demultiplexing tables ##
    # print(prop.table(table(seurat@meta.data$FinalAssignment, gsub(TRUE,"scrublet_doublet",seurat@meta.data$ScrubletDoublet) %>% gsub(FALSE,"scrublet_singlet",.)), margin = 1))
    # print(prop.table(table(seurat@meta.data$FinalAssignment, gsub(TRUE,"scrublet_doublet",seurat@meta.data$ScrubletDoublet) %>% gsub(FALSE,"scrublet_singlet",.)), margin = 2))
    # print(prop.table(table(seurat@meta.data$FinalAssignment, seurat@meta.data$MULTI_ID), margin = 1))
    # print(prop.table(table(seurat@meta.data$FinalAssignment, seurat@meta.data$MULTI_ID), margin = 2))
    # print(prop.table(table(seurat@meta.data$FinalAssignment, seurat@meta.data$MULTI_classification), margin = 1))
    # print(prop.table(table(seurat@meta.data$FinalAssignment, seurat@meta.data$MULTI_classification), margin = 2))

    # print("The seurat object before removing doublets:")
    # print(seurat)
    # Idents(seurat) <- "ScrubletDoublet"
    # seurat <- subset(seurat, idents = FALSE)
    # Idents(seurat) <- "FinalAssignment"
    # seurat <- subset(seurat,  idents = c("FSA006", "TOB421","MBE1006"))
    # Idents(seurat) <- "HTO_classification.global"
    # seurat <- subset(seurat,  idents != "Doublets")
    # print("The seurat object after removing doublets")
    # print(seurat)

    # saveRDS(seurat, paste0(outdir,"Seurat_noDoublets_meta.rds"))


if ( stage == "Hashtags"){
    pools <- dir("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/", pattern = "DRENEA")
    dirs_hash <- paste0("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Barcodes_200213_A00152_0206_AHK3HYDRXX/Hashing/", pools, "/umi_count/")
    dirs10x <- paste0("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/", pools, "/outs/filtered_feature_bc_matrix/")


    ## Read in expression data
    counts_list <- lapply(dirs10x, function(x){
        Read10X(x)
    })
    names(counts_list) <- pools

    ## Add poolnames to cell names so can easily match back if there are any cells with same barcoe
    counts_list <- lapply(names(counts_list), function(x){
        colnames(counts_list[[x]]) <- paste0(x, "_", colnames(counts_list[[x]]))
        return(counts_list[[x]])
    })
    names(counts_list) <- pools


    ### Read in the hashing data ###
    hash_list <- lapply(dirs_hash, function(x){
        Read10X(x, gene.column =1)
    })
    names(hash_list) <- pools

    hash_list <- lapply(hash_list, function(x){
        x <- as.data.frame(t(as.data.frame(x)))
        x$Barcode <- rownames(x)
        return(x)
    })
    names(hash_list) <- pools

    hash_list <- lapply(names(hash_list), function(x){
        tmp <- as.data.frame(t(hash_list[[x]]))
        colnames(tmp) <- hash_list[[x]]$Barcode
        tmp <- tmp[1:(nrow(tmp)-1),]
    })
    names(hash_list) <- pools

    ### Add pool names to barcodes so all are unique ###
    hash_list <- lapply(names(hash_list), function(x){
        colnames(hash_list[[x]]) <- paste0(x, "_", colnames(hash_list[[x]]))
        return(hash_list[[x]])
    })
    names(hash_list) <- pools

    ### Remove the unmapped row (creates noise when identifying doublets/singlets)
    hash_list <- lapply(hash_list, function(x){
        x[-c(which(rownames(x) == "unmapped")),]
    })

    ### Change the hashtag names for easy interpretation
    hash_list <- lapply(hash_list, function(x){
        rownames(x) <- gsub("A0252-TGATGGCCTATTGGG","Brisbane2", rownames(x)) %>%
            gsub("A0253-TTCCGCCTCTCTTTG", "Brisbane3", .) %>%
            gsub("A0254-AGTAAGTTCAGCGTA", "Sydney1", .) %>%
            gsub("A0255-AAGTATCGTTTCGCA", "Sydney2", .) %>%
            gsub("A0256-GGTTGCCAGATGTCA", "Sydney3", .) %>%
            gsub("A0257-TGTCTTTCCTGCCAG", "Melbourne1", .) %>%
            gsub("A0258-CTCCTCTGCAATTAC", "Melbourne2", .)
        return(x)
    })
    rownames(hash_list[["DRENEA_1"]]) <- gsub("A0251-GTCAACTCTTTAGCG","Brisbane1", rownames(hash_list[["DRENEA_1"]]))
    rownames(hash_list[["DRENEA_4"]]) <- gsub("A0251-GTCAACTCTTTAGCG","Brisbane1", rownames(hash_list[["DRENEA_4"]]))
    rownames(hash_list[["DRENEA_3"]]) <- gsub("A0251-GTCAACTCTTTAGCG","Melbourne3", rownames(hash_list[["DRENEA_3"]]))
    rownames(hash_list[["DRENEA_6"]]) <- gsub("A0251-GTCAACTCTTTAGCG","Melbourne3", rownames(hash_list[["DRENEA_6"]]))


    ## Get just the cells that are also in the hash list ##
    counts_list <- lapply(names(counts_list), function(x){
        counts_list[[x]][,colnames(hash_list[[x]])]
    })

    ##### Tried analysis with SCT (not recommended but not not recommended on Seurat website) #####
        ### saw a lot more negative cells when using either HTODemux or MULTIseqDemux, stick with the recommended normalization for hashtag demultiplexing
    # ## Make seurat object
    # seurat_list <- lapply(counts_list, function(x){
    #     CreateSeuratObject(counts = x)
    # })
    # names(seurat_list) <- pools

    # seurat_list <- lapply(seurat_list, function(x){
    #     SCTransform(x, verbose = TRUE)
    # })
    
    # ### Add hashtag data ###
    # seurat_list <- lapply(names(seurat_list), function(x){
    #     seurat_list[[x]][["HTO"]] <- CreateAssayObject(counts = hash_list[[x]])
    #     return(seurat_list[[x]])
    # })
    # names(seurat_list) <- pools

    # seurat_list <- lapply(seurat_list, function(x){
    #     NormalizeData(x, assay = "HTO", normalization.method = "CLR")
    # })

    # seurat_list <- lapply(seurat_list, function(x){
    #     HTODemux(x, assay = "HTO")
    # })

    # lapply(seurat_list, function(x){
    #     prop.table(table(x@meta.data$MULTI_ID))
    # })

    # lapply(seurat_list, function(x){
    #     prop.table(table(x@meta.data$MULTI_ID))
    # })
    
    #     lapply(seurat_list, function(x){
    #     prop.table(table(x@meta.data$HTO_classification.global))
    # })

    # lapply(seurat_list, function(x){
    #     prop.table(table(x@meta.data$HTO_classification.global))
    # })
    

    # ### Visualize results ###
    # lapply(names(seurat_list), function(x){
    #     Idents(seurat_list[[x]]) <- "MULTI_ID"
    #     p_hto_ridge_all <- RidgePlot(seurat_list[[x]], assay = "HTO", features = rownames(seurat_list[[x]][["HTO"]]), ncol = 1)
    #     p_hto_ridge_sng <- RidgePlot(subset(seurat_list[[x]], subset = MULTI_ID == c("Brisbane1", "Brisbane2", "Brisbane3", "Sydney1","Sydney2","Sydney3", "Melbourne1","Melbourne2","Melbourne3")), assay = "HTO", features = rownames(seurat_list[[x]][["HTO"]]), ncol = 1)
    #     p_hto_ridge_neg <- RidgePlot(subset(seurat_list[[x]], subset = MULTI_ID == "Negative" ), assay = "HTO", features = rownames(seurat_list[[x]][["HTO"]]), ncol = 1)
    #     p_hto_ridge_dbl <- RidgePlot(subset(seurat_list[[x]], subset = MULTI_ID == "Doublet" ), assay = "HTO", features = rownames(seurat_list[[x]][["HTO"]]), ncol = 1)

    #     ggsave(p_hto_ridge_all, filename = paste0(outdir,x,"_RidgePlot_AllCells_SCT.png"))
    #     ggsave(p_hto_ridge_sng, filename = paste0(outdir,x,"_RidgePlot_htoSinglet_SCT.png"))
    #     ggsave(p_hto_ridge_neg, filename = paste0(outdir,x,"_RidgePlot_htoNegative_SCT.png"))
    #     ggsave(p_hto_ridge_dbl, filename = paste0(outdir,x,"_RidgePlot_htoDoublet_SCT.png"))

    # })


    ##### Tried analysis with normal normalization as recommended on website #####
    ### Note: very low identification of doublets or negative cells ~3% per pool but this is not really a surprise since we single cell sorted directly before putting on 10x
    seurat_list_norm <- lapply(counts_list, function(x){
        CreateSeuratObject(counts = x)
    })
    names(seurat_list_norm) <- pools

    seurat_list_norm <- lapply(seurat_list_norm, function(x){
        tmp <- NormalizeData(x, verbose = TRUE)
        tmp <- FindVariableFeatures(tmp, selection.method = "mean.var.plot")
        tmp <- ScaleData(tmp, features = VariableFeatures(tmp))
        return(tmp)
    })
    names(seurat_list_norm) <- pools
    
    ### Add hashtag data ###
    seurat_list_norm <- lapply(names(seurat_list_norm), function(x){
        seurat_list_norm[[x]][["HTO"]] <- CreateAssayObject(counts = hash_list[[x]])
        return(seurat_list_norm[[x]])
    })
    names(seurat_list_norm) <- pools

    seurat_list_norm <- lapply(seurat_list_norm, function(x){
        NormalizeData(x, assay = "HTO", normalization.method = "RC")
    })


    seurat_list_norm <- lapply(seurat_list_norm, function(x){
        MULTIseqDemux(x, assay = "HTO", autoThresh = TRUE)
    })

    lapply(seurat_list_norm, function(x){
        prop.table(table(x@meta.data$MULTI_ID))
    })

    lapply(seurat_list_norm, function(x){
        prop.table(table(x@meta.data$MULTI_ID))
    })
    

    ### Visualize results ###
    lapply(names(seurat_list_norm), function(x){
        Idents(seurat_list_norm[[x]]) <- "MULTI_ID"
        p_hto_ridge_all <- RidgePlot(seurat_list_norm[[x]], assay = "HTO", features = rownames(seurat_list_norm[[x]][["HTO"]]), ncol = 1)
        p_hto_ridge_sng <- RidgePlot(subset(seurat_list_norm[[x]], subset = MULTI_ID == c("Brisbane1", "Brisbane2", "Brisbane3", "Sydney1","Sydney2","Sydney3", "Melbourne1","Melbourne2","Melbourne3")), assay = "HTO", features = rownames(seurat_list_norm[[x]][["HTO"]]), ncol = 1)
        p_hto_ridge_neg <- RidgePlot(subset(seurat_list_norm[[x]], subset = MULTI_ID == "Negative" ), assay = "HTO", features = rownames(seurat_list_norm[[x]][["HTO"]]), ncol = 1)
        p_hto_ridge_dbl <- RidgePlot(subset(seurat_list_norm[[x]], subset = MULTI_ID == "Doublet" ), assay = "HTO", features = rownames(seurat_list_norm[[x]][["HTO"]]), ncol = 1)

        ggsave(p_hto_ridge_all, filename = paste0(outdir,x,"_RidgePlot_AllCells.png"))
        ggsave(p_hto_ridge_sng, filename = paste0(outdir,x,"_RidgePlot_htoSinglet.png"))
        ggsave(p_hto_ridge_neg, filename = paste0(outdir,x,"_RidgePlot_htoNegative.png"))
        ggsave(p_hto_ridge_dbl, filename = paste0(outdir,x,"_RidgePlot_htoDoublet.png"))

    })

    seurat_norm <- merge(seurat_list_norm[["DRENEA_1"]], y = seurat_list_norm[2:length(seurat_list_norm)],project = "Village_Phase1")

    ### Get the metadata ###
    meta <- read_delim(paste0(dir,"/output/CompareDemultiplexing/doublet_metadata.txt"), delim = "\t")
    meta <- as.data.frame(meta)
    rownames(meta) <- meta$Barcode

    seurat_norm <- AddMetaData(seurat_norm, meta)
    seurat_norm@meta.data$Time <- ifelse((seurat_norm@meta.data$Pool == "DRENEA_1" | 
                                            seurat_norm@meta.data$Pool == "DRENEA_2" |
                                            seurat_norm@meta.data$Pool == "DRENEA_3"), "Day_0", "Day_4")

    ## Add QC metrics
    RbGeneList <- read.delim(file = "/directflow/SCCGGroupShare/projects/DrewNeavin/References/RibosomalGeneList.txt",header = F)
    print("Calculating Mt %")
    seurat_norm <- PercentageFeatureSet(seurat_norm, pattern = "^MT-", col.name = "percent.mt")
    print("Calculating Rb %")
    seurat_norm <- PercentageFeatureSet(seurat_norm, features = RbGeneList$V1, col.name = "percent.rb")

    saveRDS(seurat_norm, paste0(outdir,"AllCellsSeurat_norm_meta.rds"))
    seurat_norm <- readRDS(paste0(outdir,"AllCellsSeurat_norm_meta.rds"))

    print(prop.table(table(seurat_norm@meta.data$FinalAssignment, gsub(TRUE,"scrublet_doublet",seurat_norm@meta.data$ScrubletDoublet) %>% gsub(FALSE,"scrublet_singlet",.)), margin = 1))
    print(prop.table(table(seurat_norm@meta.data$FinalAssignment, gsub(TRUE,"scrublet_doublet",seurat_norm@meta.data$ScrubletDoublet) %>% gsub(FALSE,"scrublet_singlet",.)), margin = 2))
    print(prop.table(table(seurat_norm@meta.data$FinalAssignment, seurat_norm@meta.data$MULTI_ID), margin = 1))
    print(prop.table(table(seurat_norm@meta.data$FinalAssignment, seurat_norm@meta.data$MULTI_ID), margin = 2))
    print(prop.table(table(seurat_norm@meta.data$MULTI_ID, gsub(TRUE,"scrublet_doublet",seurat_norm@meta.data$ScrubletDoublet) %>% gsub(FALSE,"scrublet_singlet",.)), margin = 1))
    print(prop.table(table(seurat_norm@meta.data$MULTI_ID, gsub(TRUE,"scrublet_doublet",seurat_norm@meta.data$ScrubletDoublet) %>% gsub(FALSE,"scrublet_singlet",.)), margin = 2))
    print(prop.table(table(seurat_norm@meta.data$FinalAssignment, gsub(TRUE,"scrublet_doublet",seurat_norm@meta.data$ScrubletDoublet) %>% gsub(FALSE,"scrublet_singlet",.)), margin = 2))
    print(prop.table(table(seurat_norm@meta.data$FinalAssignment, paste(seurat_norm@meta.data$MULTI_ID, seurat_norm@meta.data$Time)), margin = 1))
    print(prop.table(table(seurat_norm@meta.data$FinalAssignment, paste(seurat_norm@meta.data$MULTI_ID, seurat_norm@meta.data$Time)), margin = 2))

    p_Site_by_CellNumber <- ggplot(seurat_norm@meta.data, aes(factor(paste(seurat_norm@meta.data$MULTI_ID, seurat_norm@meta.data$Time), 
            levels = c("Brisbane1 Day_0","Brisbane2 Day_0","Brisbane3 Day_0", "Brisbane1 Day_4","Brisbane2 Day_4","Brisbane3 Day_4",
                        "Sydney1 Day_0","Sydney2 Day_0","Sydney3 Day_0", "Sydney1 Day_4","Sydney2 Day_4","Sydney3 Day_4", 
                        "Melbourne1 Day_0","Melbourne2 Day_0","Melbourne3 Day_0",
                        "Melbourne1 Day_4","Melbourne2 Day_4","Melbourne3 Day_4", 
                        "Doublet Day_0","Negative Day_0", "Doublet Day_4","Negative Day_4")), fill = factor(seurat_norm@meta.data$FinalAssignment, levels = c("FSA006","MBE1006","TOB421","doublet","unassigned","combination_fail_remove")))) +
        geom_bar(stat = "count", position = "dodge") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(y="Total Cells", x = "Site Sample", fill = "Cell Line") +
        scale_fill_manual(values = cell_line_colors) 

    ggsave(p_Site_by_CellNumber, filename = paste0(outdir,"Site_by_CellNumber.png"))

    p_Site_by_CellNumber <- ggplot(seurat_norm@meta.data, aes(factor(paste(seurat_norm@meta.data$MULTI_ID, seurat_norm@meta.data$Time), 
            levels = c("Brisbane1 Day_0","Brisbane2 Day_0","Brisbane3 Day_0", "Brisbane1 Day_4","Brisbane2 Day_4","Brisbane3 Day_4",
                        "Sydney1 Day_0","Sydney2 Day_0","Sydney3 Day_0", "Sydney1 Day_4","Sydney2 Day_4","Sydney3 Day_4", 
                        "Melbourne1 Day_0","Melbourne2 Day_0","Melbourne3 Day_0",
                        "Melbourne1 Day_4","Melbourne2 Day_4","Melbourne3 Day_4", 
                        "Doublet Day_0","Negative Day_0", "Doublet Day_4","Negative Day_4")), fill = factor(seurat_norm@meta.data$FinalAssignment, levels = c("FSA006","MBE1006","TOB421","doublet","unassigned","combination_fail_remove")))) +
        geom_bar(stat = "count", position = "stack") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(y="Total Cells", x = "Site Sample", fill = "Cell Line") +
        scale_fill_manual(values = cell_line_colors) 

    ggsave(p_Site_by_CellNumber, filename = paste0(outdir,"Site_by_CellNumberStack.png"))

    p_Site_by_CellProportion <- ggplot(seurat_norm@meta.data, aes(factor(paste(seurat_norm@meta.data$MULTI_ID, seurat_norm@meta.data$Time), 
            levels = c("Brisbane1 Day_0","Brisbane2 Day_0","Brisbane3 Day_0", "Brisbane1 Day_4","Brisbane2 Day_4","Brisbane3 Day_4",
                        "Sydney1 Day_0","Sydney2 Day_0","Sydney3 Day_0", "Sydney1 Day_4","Sydney2 Day_4","Sydney3 Day_4", 
                        "Melbourne1 Day_0","Melbourne2 Day_0","Melbourne3 Day_0",
                        "Melbourne1 Day_4","Melbourne2 Day_4","Melbourne3 Day_4", 
                        "Doublet Day_0","Negative Day_0", "Doublet Day_4","Negative Day_4")), fill = factor(seurat_norm@meta.data$FinalAssignment, levels = c("FSA006","MBE1006","TOB421","doublet","unassigned","combination_fail_remove")))) +
        geom_bar(position = "fill") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(y="Total Cells", x = "Site Sample", fill = "Cell Line") +
        scale_fill_manual(values = cell_line_colors) 

    ggsave(p_Site_by_CellProportion, filename = paste0(outdir,"Site_by_CellProportion.png"))


    print("The seurat object before removing doublets:")
    print(seurat_norm)
    Idents(seurat_norm) <- "ScrubletDoublet"
    seurat_norm <- subset(seurat_norm, idents = FALSE)
    Idents(seurat_norm) <- "FinalAssignment"
    seurat_norm <- subset(seurat_norm,  idents = c("FSA006", "TOB421", "MBE1006"))
    print(seurat_norm)
    Idents(seurat_norm) <- "MULTI_ID"
    seurat_norm <- subset(seurat_norm,  idents = c("Sydney1","Brisbane1", "Melbourne1", "Brisbane2", "Sydney2", "Melbourne2", "Melbourne3", "Sydney3", "Brisbane3"))
    print("The seurat object after removing doublets")
    print(seurat_norm)

    saveRDS(seurat_norm, paste0(outdir,"Seurat_noDoublets_norm_meta.rds"))
    seurat_norm <- readRDS(paste0(outdir,"Seurat_noDoublets_norm_meta.rds"))


    p_Site_by_CellNumber <- ggplot(seurat_norm@meta.data, aes(factor(paste(seurat_norm@meta.data$MULTI_ID, seurat_norm@meta.data$Time), 
            levels = c("Brisbane1 Day_0","Brisbane2 Day_0","Brisbane3 Day_0", "Brisbane1 Day_4","Brisbane2 Day_4","Brisbane3 Day_4",
                        "Sydney1 Day_0","Sydney2 Day_0","Sydney3 Day_0", "Sydney1 Day_4","Sydney2 Day_4","Sydney3 Day_4", 
                        "Melbourne1 Day_0","Melbourne2 Day_0","Melbourne3 Day_0",
                        "Melbourne1 Day_4","Melbourne2 Day_4","Melbourne3 Day_4")), fill = factor(seurat_norm@meta.data$FinalAssignment, levels = c("FSA006","MBE1006","TOB421","doublet")))) +
        geom_bar(stat = "count", position = "dodge") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(y="Total Cells", x = "Site Sample", fill = "Cell Line") +
        scale_fill_manual(values = cell_line_colors) 

    ggsave(p_Site_by_CellNumber, filename = paste0(outdir,"Site_by_CellNumber_singlets.png"))

    p_Site_by_CellNumber <- ggplot(seurat_norm@meta.data, aes(factor(paste(seurat_norm@meta.data$MULTI_ID, seurat_norm@meta.data$Time), 
            levels = c("Brisbane1 Day_0","Brisbane2 Day_0","Brisbane3 Day_0", "Brisbane1 Day_4","Brisbane2 Day_4","Brisbane3 Day_4",
                        "Sydney1 Day_0","Sydney2 Day_0","Sydney3 Day_0", "Sydney1 Day_4","Sydney2 Day_4","Sydney3 Day_4", 
                        "Melbourne1 Day_0","Melbourne2 Day_0","Melbourne3 Day_0",
                        "Melbourne1 Day_4","Melbourne2 Day_4","Melbourne3 Day_4", 
                        "Doublet Day_0","Negative Day_0", "Doublet Day_4","Negative Day_4")), fill = factor(seurat_norm@meta.data$FinalAssignment, levels = c("FSA006","MBE1006","TOB421","doublet","unassigned","combination_fail_remove")))) +
        geom_bar(stat = "count", position = "stack") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(y="Total Cells", x = "Site Sample", fill = "Cell Line") +
        scale_fill_manual(values = cell_line_colors) 

    ggsave(p_Site_by_CellNumber, filename = paste0(outdir,"Site_by_CellNumberStack_singlets.png"))

    p_Site_by_CellProportion <- ggplot(seurat_norm@meta.data, aes(factor(paste(seurat_norm@meta.data$MULTI_ID, seurat_norm@meta.data$Time), 
            levels = c("Brisbane1 Day_0","Brisbane2 Day_0","Brisbane3 Day_0", "Brisbane1 Day_4","Brisbane2 Day_4","Brisbane3 Day_4",
                        "Sydney1 Day_0","Sydney2 Day_0","Sydney3 Day_0", "Sydney1 Day_4","Sydney2 Day_4","Sydney3 Day_4", 
                        "Melbourne1 Day_0","Melbourne2 Day_0","Melbourne3 Day_0",
                        "Melbourne1 Day_4","Melbourne2 Day_4","Melbourne3 Day_4", 
                        "Doublet Day_0","Negative Day_0", "Doublet Day_4","Negative Day_4")), fill = factor(seurat_norm@meta.data$FinalAssignment, levels = c("FSA006","MBE1006","TOB421","doublet")))) +
        geom_bar(position = "fill") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(y="Total Cells", x = "Site Sample", fill = "Cell Line") +
        scale_fill_manual(values = cell_line_colors) 

    ggsave(p_Site_by_CellProportion, filename = paste0(outdir,"Site_by_CellProportion_singlets.png"))



} else if (stage[1]=="QC"){
    pools <- dir("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/", pattern = "DRENEA")
    seurat <- readRDS("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Seurat_w_Hash/Hashtags/Seurat_noDoublets_norm_meta.rds")

    # ## Calculate MADs
    # seurat <- mad_function(seurat, "percent.mt",3)
    # seurat <- mad_function(seurat, "percent.rb",3)
    # seurat <- mad_function(seurat, "nCount_RNA", 3)
    # seurat <- mad_function(seurat, "nFeature_RNA", 3)


    # plot_mt_pct <- VlnPlot(seurat, features = c( "percent.mt"), group.by = "Pool", pt.size = 0, cols = PoolColors) +
    #                 geom_hline(yintercept=0, linetype="dashed", color = "black") +
    #                 geom_hline(yintercept=13.66416, linetype="dashed", color = "black")
    # ggsave(plot_mt_pct, filename = paste0(outdir,"Mt_pct_vln.png"))

    # plot_rb_pct <- VlnPlot(seurat, features = c( "percent.rb"), group.by = "Pool", pt.size = 0, cols = PoolColors) +
    #                 geom_hline(yintercept=15.479573, linetype="dashed", color = "black") +
    #                 geom_hline(yintercept=37.03238, linetype="dashed", color = "black")
    # ggsave(plot_rb_pct, filename = paste0(outdir,"Rb_pct_vln.png"))

    # plot_n_count <- VlnPlot(seurat, features = c( "nCount_RNA"), group.by = "Pool", pt.size = 0, cols = PoolColors) +
    #                 geom_hline(yintercept=0, linetype="dashed", color = "black") +
    #                 geom_hline(yintercept=44381.41, linetype="dashed", color = "black")
    # ggsave(plot_n_count, filename = paste0(outdir,"N_count_vln.png"))

    # plot_nFeature_RNA <- VlnPlot(seurat, features = c( "nFeature_RNA"), group.by = "Pool", pt.size = 0, cols = PoolColors) +
    #                 geom_hline(yintercept=1231.853, linetype="dashed", color = "black") +
    #                 geom_hline(yintercept=7512.147, linetype="dashed", color = "black")
    # ggsave(plot_nFeature_RNA, filename = paste0(outdir,"nFeature_RNA_vln.png"))

    # lib_mt <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Pool", cols = alpha(brewer.pal(n = 6, name = "Dark2"), 0.25)) +
    #                 geom_vline(xintercept=0, linetype="dashed", color = "black") +
    #                 geom_vline(xintercept=44548.25, linetype="dashed", color = "black") +
    #                 geom_hline(yintercept=0, linetype="dashed", color = "black") +
    #                 geom_hline(yintercept=13.732771, linetype="dashed", color = "black") 
    # ggsave(lib_mt, filename = paste0(outdir,"NoDoublets_lib_mt.png"))
    # lib_genes <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Pool", cols = alpha(brewer.pal(n = 6, name = "Dark2"), 0.25)) +
    #                 geom_vline(xintercept=0, linetype="dashed", color = "black") +
    #                 geom_vline(xintercept=44381.41, linetype="dashed", color = "black") +
    #                 geom_hline(yintercept=1231.853, linetype="dashed", color = "black") +
    #                 geom_hline(yintercept=7512.147, linetype="dashed", color = "black")
    # ggsave(lib_genes, filename = paste0(outdir,"NoDoublets_lib_genes.png"))


    # ### Vioplin plots by cell type and site ###
    # p_violin_mt_pct_line_site <- ggplot(seurat@meta.data, aes(factor(Time, levels = c("Day_0", "Day_4")), percent.mt, fill = factor(FinalAssignment, levels = c("FSA006", "MBE1006","TOB421")))) +
    #     geom_violin() +
    #     theme_classic() +
    #     labs(y="Mitochondrial Percent", x = "Sample Day", fill = "Cell Line") +
    #     scale_fill_manual(values = cell_line_colors) +
    #     facet_wrap(vars( MULTI_ID))
    # ggsave(p_violin_mt_pct_line_site, filename = paste0(outdir,"violin_mt_pct_line_site.png"))

    # p_violin_rb_pct_line_site <- ggplot(seurat@meta.data, aes(factor(Time, levels = c("Day_0", "Day_4")), percent.rb, fill = factor(FinalAssignment, levels = c("FSA006", "MBE1006","TOB421")))) +
    #     geom_violin() +
    #     theme_classic() +
    #     labs(y="Ribosomal Percent", x = "Sample Day", fill = "Cell Line") +
    #     scale_fill_manual(values = cell_line_colors) +
    #     facet_wrap(vars( MULTI_ID))
    # ggsave(p_violin_rb_pct_line_site, filename = paste0(outdir,"violin_rb_pct_line_site.png"))

    # p_violin_n_count_line_site <- ggplot(seurat@meta.data, aes(factor(Time, levels = c("Day_0", "Day_4")), nCount_RNA, fill = factor(FinalAssignment, levels = c("FSA006", "MBE1006","TOB421")))) +
    #     geom_violin() +
    #     theme_classic() +
    #     labs(y="Number UMIs", x = "Sample Day", fill = "Cell Line") +
    #     scale_fill_manual(values = cell_line_colors) +
    #     facet_wrap(vars( MULTI_ID))
    # ggsave(p_violin_n_count_line_site, filename = paste0(outdir,"violin_n_count_line_site.png"))

    # p_violin_n_feature_line_site <- ggplot(seurat@meta.data, aes(factor(Time, levels = c("Day_0", "Day_4")), nFeature_RNA, fill = factor(FinalAssignment, levels = c("FSA006", "MBE1006","TOB421")))) +
    #     geom_violin() +
    #     theme_classic() +
    #     labs(y="Number Features", x = "Sample Day", fill = "Cell Line") +
    #     scale_fill_manual(values = cell_line_colors) +
    #     facet_wrap(vars( MULTI_ID))
    # ggsave(p_violin_n_feature_line_site, filename = paste0(outdir,"violin_n_feature_line_site.png"))



    # ##### Remove the outliers #####
    # print(seurat)
    # seurat <- subset(seurat, subset = percent.mt_mad == "NotOutlier") 
    # seurat <- subset(seurat, subset = percent.rb_mad == "NotOutlier")
    # seurat <- subset(seurat, subset = nCount_RNA_mad == "NotOutlier")
    # seurat <- subset(seurat, subset = nFeature_RNA_mad == "NotOutlier")
    # print(seurat)
    # saveRDS(seurat, paste0(outdir,"seurat_norm_NoOutliers.rds"))


    # plot_mt_pct <- VlnPlot(seurat, features = c( "percent.mt"), group.by = "Pool", pt.size = 0, cols = PoolColors)
    # ggsave(plot_mt_pct, filename = paste0(outdir,"Mt_pct_vln_NoOutliers.png"))

    # plot_rb_pct <- VlnPlot(seurat, features = c( "percent.rb"), group.by = "Pool", pt.size = 0, cols = PoolColors)
    # ggsave(plot_rb_pct, filename = paste0(outdir,"Rb_pct_vln_NoOutliers.png"))

    # plot_n_count <- VlnPlot(seurat, features = c( "nCount_RNA"), group.by = "Pool", pt.size = 0, cols = PoolColors)
    # ggsave(plot_n_count, filename = paste0(outdir,"N_count_vln_NoOutliers.png"))

    # plot_nFeature_RNA <- VlnPlot(seurat, features = c( "nFeature_RNA"), group.by = "Pool", pt.size = 0, cols = PoolColors)
    # ggsave(plot_nFeature_RNA, filename = paste0(outdir,"nFeature_RNA_vln_NoOutliers.png"))

    # lib_mt <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Pool", cols = alpha(brewer.pal(n = 6, name = "Dark2"), 0.25))
    # ggsave(lib_mt, filename = paste0(outdir,"NoDoublets_lib_mt_NoOutliers.png"))
    # lib_genes <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Pool", cols = alpha(brewer.pal(n = 6, name = "Dark2"), 0.25))
    # ggsave(lib_genes, filename = paste0(outdir,"NoDoublets_lib_genes_NoOutliers.png"))


    # ### Vioplin plots by cell type and site ###
    # seurat <- readRDS(paste0(outdir,"seurat_norm_NoOutliers.rds"))
    # p_violin_mt_pct_line_site <- ggplot(seurat@meta.data, aes(factor(Time, levels = c("Day_0", "Day_4")), percent.mt, fill = factor(FinalAssignment, levels = c("FSA006", "MBE1006","TOB421")))) +
    #     geom_violin() +
    #     theme_classic() +
    #     labs(y="Mitochondrial Percent", x = "Sample Day", fill = "Cell Line") +
    #     scale_fill_manual(values = cell_line_colors) +
    #     facet_wrap(vars( MULTI_ID))
    # ggsave(p_violin_mt_pct_line_site, filename = paste0(outdir,"violin_mt_pct_line_site_NoOutliers.png"))

    # p_violin_rb_pct_line_site <- ggplot(seurat@meta.data, aes(factor(Time, levels = c("Day_0", "Day_4")), percent.rb, fill = factor(FinalAssignment, levels = c("FSA006", "MBE1006","TOB421")))) +
    #     geom_violin() +
    #     theme_classic() +
    #     labs(y="Ribosomal Percent", x = "Sample Day", fill = "Cell Line") +
    #     scale_fill_manual(values = cell_line_colors) +
    #     facet_wrap(vars( MULTI_ID))
    # ggsave(p_violin_rb_pct_line_site, filename = paste0(outdir,"violin_rb_pct_line_site_NoOutliers.png"))

    # p_violin_n_count_line_site <- ggplot(seurat@meta.data, aes(factor(Time, levels = c("Day_0", "Day_4")), nCount_RNA, fill = factor(FinalAssignment, levels = c("FSA006", "MBE1006","TOB421")))) +
    #     geom_violin() +
    #     theme_classic() +
    #     labs(y="Number UMIs", x = "Sample Day", fill = "Cell Line") +
    #     scale_fill_manual(values = cell_line_colors) +
    #     facet_wrap(vars( MULTI_ID))
    # ggsave(p_violin_n_count_line_site, filename = paste0(outdir,"violin_n_count_line_site_NoOutliers.png"))

    # p_violin_n_feature_line_site <- ggplot(seurat@meta.data, aes(factor(Time, levels = c("Day_0", "Day_4")), nFeature_RNA, fill = factor(FinalAssignment, levels = c("FSA006", "MBE1006","TOB421")))) +
    #     geom_violin() +
    #     theme_classic() +
    #     labs(y="Number Features", x = "Sample Day", fill = "Cell Line") +
    #     scale_fill_manual(values = cell_line_colors) +
    #     facet_wrap(vars( MULTI_ID))
    # ggsave(p_violin_n_feature_line_site, filename = paste0(outdir,"violin_n_feature_line_site_NoOutliers.png"))


} else if (stage == "CellCycle"){
    print("Starting cell cycle step")
    hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
    seurat <- readRDS(paste0(dir,"/output/Seurat_w_Hash/QC/seurat_norm_NoOutliers.rds"))
    GeneConversion <- read.delim("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/DRENEA_1/outs/filtered_feature_bc_matrix/features.tsv.gz", header = F)
    GeneConversion$V2 <- gsub("_","-",GeneConversion$V2)
    rownames(GeneConversion) <- make.unique(GeneConversion$V2)
    RbGeneList <- read.delim(file = "/directflow/SCCGGroupShare/projects/DrewNeavin/References/RibosomalGeneList.txt",header = F)

    countsENSG <- seurat[["RNA"]]@counts # Create a counts matrix that has the ENSG gene clasifiers so can determine cell cycle phase
    print("The size of the counts object in genes x cells is:")
    print(dim(countsENSG))

    ### Test to see if have all the same rownames between the two dataframes ###
    print("Test if the rownames of the counts matrix and the Gene Conversion matrix are the same:")
    print(all(rownames(countsENSG) == rownames(GeneConversion)))

    ### Reassign the row names to be ENSG IDs ###
    row.names(countsENSG) <- GeneConversion$V1

    ### Run the cell cycle identification ###
    print("Starting cell cycle determination")
    assigned <- cyclone(countsENSG, pairs=hs.pairs)  #Note, this takes hours
    table(assigned$phases)
    write.table(assigned, file = paste0(outdir,"CellCycleProportions.txt"), quote = F, sep = "\t") #Save so that can read in and don't have to wait to recompute again


    assigned <- read.table(paste0(outdir,"CellCycleProportions.txt"), sep = "\t")
    assigned <- as.data.frame(assigned)
    rownames(assigned) <- colnames(seurat)
    write.table(assigned, file = paste0(outdir,"CellCycleProportions.txt"), quote = F, sep = "\t") #Save so that can read in and don't have to wait to recompute again

    seurat <- AddMetaData(seurat, assigned)
    saveRDS(seurat, paste0(outdir,"seurat_norm_NoOutliers_wCellCycle.rds"))




} else if ( stage == "Old_scRNAseq_data"){
    # seurat <- readRDS("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Seurat_w_Hash/QC/seurat_norm_NoOutliers.rds")
    # pool_file_list <- list.files(paste0("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Old_scRNAseq_expression/"), pattern = "RawQC_ascendObj_before_filtering_Pool")

    # pool_list <- lapply(pool_file_list, function(x){
    #     readRDS(paste0("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Old_scRNAseq_expression/",x))
    # })
    # names(pool_list) <- gsub("RawQC_ascendObj_before_filtering_","", pool_file_list) %>% gsub(".RDS","", .)

    # FSA <- counts(pool_list[["Pool15"]])[,which(pool_list[["Pool15"]]@colInfo$batch == "22_FSA")]
    # colnames(FSA) <- paste0("FSA006_",colnames(FSA))

    # TOB <- counts(pool_list[["Pool14"]])[,which(pool_list[["Pool14"]]@colInfo$batch == "36_TOB00421_i_E8")]
    # colnames(TOB) <- paste0("TOB421_",colnames(TOB))

    # MBE <- counts(pool_list[["Pool13"]])[,which(pool_list[["Pool13"]]@colInfo$batch == "29_MBE")]
    # colnames(MBE) <- paste0("MBE1006_",colnames(MBE))

    # counts <- cbind(FSA,TOB,MBE)
    # saveRDS(counts, paste0(outdir,"Old_iPSC_counts.RDS"))

    # old_seurat <- CreateSeuratObject(counts = counts)
    # old_seurat$FinalAssignment <- gsub("_.*","",colnames(old_seurat))

    # ## Add QC metrics
    # RbGeneList <- read.delim(file = "/directflow/SCCGGroupShare/projects/DrewNeavin/References/RibosomalGeneList.txt",header = F)
    # print("Calculating Mt %")
    # old_seurat <- PercentageFeatureSet(old_seurat, pattern = "^MT-", col.name = "percent.mt")
    # print("Calculating Rb %")
    # RbGeneList <- RbGeneList[RbGeneList$V1 %in% rownames(old_seurat),]
    # old_seurat <- PercentageFeatureSet(old_seurat, features = RbGeneList, col.name = "percent.rb")
    
    # ## Calculate MADs
    # old_seurat <- mad_function(old_seurat, "percent.mt",3)
    # old_seurat <- mad_function(old_seurat, "percent.rb",3)
    # old_seurat <- mad_function(old_seurat, "nCount_RNA", 3)
    # old_seurat <- mad_function(old_seurat, "nFeature_RNA", 3)


    # ### Plot QC metrics ###
    # plot_mt_pct <- VlnPlot(old_seurat, features = c( "percent.mt"), group.by = "FinalAssignment", pt.size = 0, cols = cell_line_colors) +
    #                 geom_hline(yintercept=0, linetype="dashed", color = "black") +
    #                 geom_hline(yintercept=6.207164, linetype="dashed", color = "black")
    # ggsave(plot_mt_pct, filename = paste0(outdir,"Mt_pct_vln.png"))

    # plot_rb_pct <- VlnPlot(old_seurat, features = c( "percent.rb"), group.by = "FinalAssignment", pt.size = 0, cols = cell_line_colors) +
    #                 geom_hline(yintercept=16.78718, linetype="dashed", color = "black") +
    #                 geom_hline(yintercept=37.06345, linetype="dashed", color = "black")
    # ggsave(plot_rb_pct, filename = paste0(outdir,"Rb_pct_vln.png"))

    # plot_n_count <- VlnPlot(old_seurat, features = c( "nCount_RNA"), group.by = "FinalAssignment", pt.size = 0, cols = cell_line_colors) +
    #                 geom_hline(yintercept=0, linetype="dashed", color = "black") +
    #                 geom_hline(yintercept=12084.99, linetype="dashed", color = "black")
    # ggsave(plot_n_count, filename = paste0(outdir,"N_count_vln.png"))

    # plot_nFeature_RNA <- VlnPlot(old_seurat, features = c( "nFeature_RNA"), group.by = "FinalAssignment", pt.size = 0, cols = cell_line_colors) +
    #                 geom_hline(yintercept=230.8951, linetype="dashed", color = "black") +
    #                 geom_hline(yintercept=3749.105, linetype="dashed", color = "black")
    # ggsave(plot_nFeature_RNA, filename = paste0(outdir,"nFeature_RNA_vln.png"))

    # lib_mt <- FeatureScatter(old_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "FinalAssignment", cols = alpha(brewer.pal(n = 6, name = "Dark2"), 0.25)) +
    #                 geom_vline(xintercept=0, linetype="dashed", color = "black") +
    #                 geom_vline(xintercept=44548.25, linetype="dashed", color = "black") +
    #                 geom_hline(yintercept=0, linetype="dashed", color = "black") +
    #                 geom_hline(yintercept=13.732771, linetype="dashed", color = "black") 
    # ggsave(lib_mt, filename = paste0(outdir,"NoDoublets_lib_mt.png"))
    # lib_genes <- FeatureScatter(old_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "FinalAssignment", cols = alpha(brewer.pal(n = 6, name = "Dark2"), 0.25)) +
    #                 geom_vline(xintercept=0, linetype="dashed", color = "black") +
    #                 geom_vline(xintercept=44381.41, linetype="dashed", color = "black") +
    #                 geom_hline(yintercept=230.8951, linetype="dashed", color = "black") +
    #                 geom_hline(yintercept=3749.105, linetype="dashed", color = "black")
    # ggsave(lib_genes, filename = paste0(outdir,"NoDoublets_lib_genes.png"))


    # ### Remove outliers ###
    # print(old_seurat)
    # old_seurat <- subset(old_seurat, subset = percent.mt_mad == "NotOutlier") 
    # old_seurat <- subset(old_seurat, subset = percent.rb_mad == "NotOutlier")
    # old_seurat <- subset(old_seurat, subset = nCount_RNA_mad == "NotOutlier")
    # old_seurat <- subset(old_seurat, subset = nFeature_RNA_mad == "NotOutlier")
    # print(old_seurat)
    # saveRDS(old_seurat, paste0(outdir,"old_seurat_norm_NoOutliers.rds"))


    # seurat_combined <- merge(seurat, y = old_seurat,project = "Combined_iPSC")

    # seurat_combined@meta.data$MULTI_ID[is.na(seurat_combined@meta.data$MULTI_ID)] <- "Old"
    # seurat_combined@meta.data$Pool[is.na(seurat_combined@meta.data$Pool)] <- "Old"
    # seurat_combined@meta.data$Time[is.na(seurat_combined@meta.data$Time)] <- "Old"
    # seurat_combined@meta.data$Site <- gsub("\\d","", seurat_combined@meta.data$MULTI_ID)

    # seurat_combined <- SCTransform(seurat_combined, verbose = TRUE)
    # seurat_combined <- RunPCA(seurat_combined, npcs = 100)
    # seurat_combined <- RunUMAP(seurat_combined, dims = 1:100, verbose = TRUE)
    # seurat_combined <- FindNeighbors(seurat_combined, dims = 1:100, verbose = TRUE)

    # saveRDS(seurat_combined, paste0(outdir,"seurat_combined_SCT_noCov.rds"))


    # umap <- DimPlot(seurat_combined, reduction = "umap", group.by = "Pool", pt.size = .01, split.by = 'Pool')
    # ggsave(umap, filename = paste0(outdir, "umap_pool_noCov.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat_combined, reduction = "umap", group.by = "Time", pt.size = .01, split.by = 'Time')
    # ggsave(umap, filename = paste0(outdir, "umap_Time_noCov.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat_combined, reduction = "umap", group.by = "Site", pt.size = .01, split.by = 'Site')
    # ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate_noCov.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat_combined, reduction = "umap", group.by = "FinalAssignment", pt.size = .01, split.by = 'FinalAssignment')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine_noCov.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat_combined, reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'Site')
    # ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate_Time_noCov.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat_combined, reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'FinalAssignment')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time_noCov.png"), height = 4, width = 10)

    # seurat_combined$Time_FinalAssignment <- paste0(seurat_combined$Time, seurat_combined$FinalAssignment)
    # umap <- DimPlot(seurat_combined, reduction = "umap", group.by = "Time_FinalAssignment", pt.size = .01, split.by = 'Site')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time_Site_noCov.png"), height = 4, width = 10)

    # UMAPmt <- FeaturePlot(seurat_combined, reduction = "umap", feature = "percent.mt")
    # ggsave(UMAPmt, filename = paste0(outdir,"UMAPmt_noCov.png"))

    # UMAPrb <- FeaturePlot(seurat_combined, reduction = "umap", feature = "percent.rb")
    # ggsave(UMAPrb, filename = paste0(outdir,"UMAPrb_noCov.png"))


    ##### Use SCTransform with covariates #####
    # seurat_combined <- readRDS(paste0(outdir,"seurat_combined_SCT_noCov.rds"))
    # seurat_combined <- SCTransform(seurat_combined, verbose = TRUE, vars.to.regress = c("Site","percent.mt", "percent.rb","nFeature_RNA","nCount_RNA"))
    # seurat_combined <- RunPCA(seurat_combined, npcs = 100)
    # seurat_combined <- RunUMAP(seurat_combined, dims = 1:100, verbose = TRUE)
    # seurat_combined <- FindNeighbors(seurat_combined, dims = 1:100, verbose = TRUE)

    # saveRDS(seurat_combined, paste0(outdir,"seurat_combined_SCT_SiteMtRbNfeaturesNcountCov.rds"))
    # seurat_combined <- readRDS(paste0(outdir,"seurat_combined_SCT_SiteMtRbNfeaturesNcountCov.rds"))

    # umap <- DimPlot(seurat_combined, reduction = "umap", group.by = "Pool", pt.size = .01, split.by = 'Pool')
    # ggsave(umap, filename = paste0(outdir, "umap_pool_SiteMtRbNfeaturesNcountCov.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat_combined, reduction = "umap", group.by = "Time", pt.size = .01, split.by = 'Time')
    # ggsave(umap, filename = paste0(outdir, "umap_Time.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat_combined, reduction = "umap", group.by = "MULTI_ID", pt.size = .01, split.by = 'Site')
    # ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate_SiteMtRbNfeaturesNcountCov.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat_combined, reduction = "umap", group.by = "FinalAssignment", pt.size = .01, split.by = 'FinalAssignment')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine_SiteMtRbNfeaturesNcountCov.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat_combined, reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'Site')
    # ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate_Time_SiteMtRbNfeaturesNcountCov.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat_combined, reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'FinalAssignment')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time_SiteMtRbNfeaturesNcountCov.png"), height = 4, width = 10)

    # seurat_combined$Time_FinalAssignment <- paste0(seurat_combined$Time, seurat_combined$FinalAssignment)
    # umap <- DimPlot(seurat_combined, reduction = "umap", group.by = "Time_FinalAssignment", pt.size = .01, split.by = 'Site')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time_Site_SiteMtRbNfeaturesNcountCov.png"), height = 4, width = 10)

    # UMAPmt <- FeaturePlot(seurat_combined, reduction = "umap", feature = "percent.mt")
    # ggsave(UMAPmt, filename = paste0(outdir,"UMAPmt_SiteMtRbNfeaturesNcountCov.png"))

    # UMAPrb <- FeaturePlot(seurat_combined, reduction = "umap", feature = "percent.rb")
    # ggsave(UMAPrb, filename = paste0(outdir,"UMAPrb_SiteMtRbNfeaturesNcountCov.png"))


    ##### Use Harmony with covariates #####
    # seurat_combined <- readRDS(paste0(outdir,"seurat_combined_SCT_noCov.rds"))
    # DefaultAssay(seurat_combined) <- "RNA"
    # seurat_combined@meta.data$Site <- gsub("\\d","", seurat_combined@meta.data$MULTI_ID)
    # seurat_combined <- NormalizeData(seurat_combined, verbose = FALSE) %>%
    #     FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    #     ScaleData(verbose = FALSE, vars.to.regress = c("percent.mt", "percent.rb","nFeature_RNA","nCount_RNA")) %>% 
    #     RunPCA(pc.genes = seurat_combined@var.genes, npcs = 20, verbose = FALSE)
    # seurat_combined <- RunHarmony(seurat_combined,"Site", assay.use = "RNA")

    # ### Plot the metrics for the results ###
    # p1 <- DimPlot(object = seurat_combined, reduction = "harmony", pt.size = .1, group.by = "Pool", do.return = TRUE)
    # ggsave(p1, filename = paste0(outdir,"harmony_PCA.png"))

    # p2 <- VlnPlot(object = seurat_combined, features = "harmony_1", group.by = "Pool", pt.size = .1)
    # ggsave(p2, filename = paste0(outdir,"harmony_violin.png"))

    # seurat_combined <- RunUMAP(seurat_combined, reduction = "harmony", dims = 1:20)
    # seurat_combined <- FindNeighbors(seurat_combined, reduction = "harmony", dims = 1:20)

    # saveRDS(seurat_combined, paste0(outdir,"seurat_combined_SCT_harmonySite.rds"))

    # umap <- DimPlot(seurat_combined, reduction = "umap", group.by = "Pool", pt.size = .01, split.by = 'Pool')
    # ggsave(umap, filename = paste0(outdir, "umap_pool_harmonySite.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat_combined, reduction = "umap", group.by = "Time", pt.size = .01, split.by = 'Time')
    # ggsave(umap, filename = paste0(outdir, "umap_Time.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat_combined, reduction = "umap", group.by = "MULTI_ID", pt.size = .01, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate_harmonySite.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat_combined, reduction = "umap", group.by = "FinalAssignment", pt.size = .01, split.by = 'FinalAssignment')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine_harmonySite.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat_combined, reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate_Time_harmonySite.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat_combined, reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'FinalAssignment')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time_harmonySite.png"), height = 4, width = 10)

    # seurat_combined$Time_FinalAssignment <- paste0(seurat_combined$Time, seurat_combined$FinalAssignment)
    # umap <- DimPlot(seurat_combined, reduction = "umap", group.by = "Time_FinalAssignment", pt.size = .01, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time_Site_harmonySite.png"), height = 4, width = 10)

    # UMAPmt <- FeaturePlot(seurat_combined, reduction = "umap", feature = "percent.mt")
    # ggsave(UMAPmt, filename = paste0(outdir,"UMAPmt_harmonySite.png"))

    # UMAPrb <- FeaturePlot(seurat_combined, reduction = "umap", feature = "percent.rb")
    # ggsave(UMAPrb, filename = paste0(outdir,"UMAPrb_harmonySite.png"))



    ### Integration for Site
    seurat_combined <- readRDS(paste0(outdir,"seurat_combined_SCT_noCov.rds"))
    seurat_combined@meta.data$Site <- gsub("\\d","", seurat_combined@meta.data$MULTI_ID)
    DefaultAssay(seurat_combined) <- "RNA"
    print(head(seurat_combined@meta.data$Site))
    seurat_list <- SplitObject(seurat_combined, split.by = "Site")

    seurat_list <- lapply(seurat_list, function(x) {
        SCTransform(x, verbose = TRUE, vars.to.regress = c("percent.mt", "percent.rb","nFeature_RNA","nCount_RNA"))
    })

    seurat_features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)
    seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = seurat_features)
    seurat_anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT", 
        anchor.features = seurat_features)
    seurat_integrated <- IntegrateData(anchorset = seurat_anchors, normalization.method = "SCT")

    seurat_integrated <- RunPCA(object = seurat_integrated, verbose = FALSE)
    seurat_integrated <- RunUMAP(object = seurat_integrated, dims = 1:30)

    seurat_integrated <- FindNeighbors(seurat_integrated, reduction = "pca", assay = "integrated", dims = 1:20)

    saveRDS(seurat_integrated,paste0(outdir, "NoOutliers_Site_SeuratIntegrated.rds"))


    plots <- DimPlot(seurat_integrated, group.by = c("Time", "FinalAssignment"), combine = FALSE)
    plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 4, 
        byrow = TRUE, override.aes = list(size = 2.5))))
    ggsave(CombinePlots(plots), filename = paste0(outdir,"umap_integrationSite.png"),width = 12, height = 9)

        umap <- DimPlot(seurat_integrated, reduction = "umap", group.by = "Pool", pt.size = .01, split.by = 'Pool')
    ggsave(umap, filename = paste0(outdir, "umap_pool_integrationSite.png"), height = 4, width = 10)


    umap <- DimPlot(seurat_integrated, reduction = "umap", group.by = "Time", pt.size = .01, split.by = 'Time')
    ggsave(umap, filename = paste0(outdir, "umap_Time_integrationSite.png"), height = 4, width = 10)

    umap <- DimPlot(seurat_integrated, reduction = "umap", group.by = "MULTI_ID", pt.size = .01, split.by = 'MULTI_ID')
    ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate_integrationSite.png"), height = 4, width = 10)

    umap <- DimPlot(seurat_integrated, reduction = "umap", group.by = "FinalAssignment", pt.size = .01, split.by = 'FinalAssignment')
    ggsave(umap, filename = paste0(outdir, "umap_CellLine_integrationSite.png"), height = 4, width = 10)

    umap <- DimPlot(seurat_integrated, reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'MULTI_ID')
    ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate_Time_integrationSite.png"), height = 4, width = 10)

    umap <- DimPlot(seurat_integrated, reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'FinalAssignment')
    ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time_integrationSite.png"), height = 4, width = 10)

    seurat_integrated$Time_FinalAssignment <- paste0(seurat_integrated$Time, seurat_integrated$FinalAssignment)
    umap <- DimPlot(seurat_integrated, reduction = "umap", group.by = "Time_FinalAssignment", pt.size = .001, split.by = 'MULTI_ID')
    ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time_Site_integrationSite.png"), height = 4, width = 10)

    UMAPmt <- FeaturePlot(seurat_integrated, reduction = "umap", feature = "percent.mt")
    ggsave(UMAPmt, filename = paste0(outdir,"UMAPmt_integrationSite.png"))

    UMAPrb <- FeaturePlot(seurat_integrated, reduction = "umap", feature = "percent.rb")
    ggsave(UMAPrb, filename = paste0(outdir,"UMAPrb_integrationSite.png"))



} else if ( stage == "SCTnormalization_no_covatiates"){
    # seurat <- readRDS(paste0(dir,"output/Seurat_w_Hash/QC/seurat_norm_NoOutliers.rds"))
    # seurat <- SCTransform(seurat, verbose = TRUE, return.only.var.genes = FALSE)
    # seurat <- RunPCA(seurat, npcs = 100)
    # seurat <- RunUMAP(seurat, dims = 1:100, verbose = TRUE)
    # seurat <- FindNeighbors(seurat, dims = 1:100, verbose = TRUE)

    # seurat@meta.data$Site <- gsub("\\d","",seurat@meta.data$MULTI_ID)
    # saveRDS(seurat, paste0(outdir, "seurat_SCT_noCov.rds"))
    seurat <- readRDS(paste0(outdir, "seurat_SCT_noCov.rds"))

    assigned <- read.table(paste0(dir,"output/Seurat_w_Hash/CellCycle/CellCycleProportions.txt"), sep = "\t")
    seurat <- AddMetaData(seurat, assigned)

    # plots <- DimPlot(seurat, group.by = c("Time", "FinalAssignment", "Site"), combine = FALSE)
    # plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 4, 
    #     byrow = TRUE, override.aes = list(size = 2.5))))
    # ggsave(CombinePlots(plots, ncol = 3), filename = paste0(outdir,"umap.png"),width = 20, height = 9)


    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Pool", pt.size = .01, split.by = 'Pool')
    # ggsave(umap, filename = paste0(outdir, "umap_pool.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .01, split.by = 'Time')
    # ggsave(umap, filename = paste0(outdir, "umap_Time.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "MULTI_ID", pt.size = .01, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "FinalAssignment", pt.size = .01, split.by = 'FinalAssignment')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate_Time.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'FinalAssignment')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time.png"), height = 4, width = 10)

    # seurat$Time_FinalAssignment <- paste0(seurat$Time, seurat$FinalAssignment)
    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time_FinalAssignment", pt.size = .01, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time_Site.png"), height = 4, width = 10)

    umap <- DimPlot(seurat,  pt.size = .01, group.by = c("phases"))  
    ggsave(umap, filename = paste0(outdir, "umap_CellCycle.png"))



} else if ( stage == "SCTnormalization_SiteRegressed"){
    # seurat <- readRDS(paste0(dir,"output/Seurat_w_Hash/QC/seurat_norm_NoOutliers.rds"))
    # seurat@meta.data$Site <- gsub("\\d","", seurat@meta.data$MULTI_ID)
    # print(head(seurat@meta.data$Site))
    # seurat <- SCTransform(seurat, vars.to.regress = "Site", verbose = TRUE)
    # seurat <- RunPCA(seurat, npcs = 100)
    # seurat <- RunUMAP(seurat, dims = 1:100, verbose = TRUE)
    # seurat <- FindNeighbors(seurat, dims = 1:100, verbose = TRUE)

    # saveRDS(seurat, paste0(outdir, "seurat_SCT_SiteRegressed.rds"))
    seurat <- readRDS(paste0(outdir, "seurat_SCT_SiteRegressed.rds"))

    assigned <- read.table(paste0(dir,"output/Seurat_w_Hash/CellCycle/CellCycleProportions.txt"), sep = "\t")
    seurat <- AddMetaData(seurat, assigned)


    # plots <- DimPlot(seurat, group.by = c("Time", "FinalAssignment", "Site"), combine = FALSE)
    # plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 4, 
    #     byrow = TRUE, override.aes = list(size = 2.5))))
    # ggsave(CombinePlots(plots, ncol = 3), filename = paste0(outdir,"umap.png"),width = 20, height = 9)


    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Pool", pt.size = .01, split.by = 'Pool')
    # ggsave(umap, filename = paste0(outdir, "umap_pool.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .01, split.by = 'Time')
    # ggsave(umap, filename = paste0(outdir, "umap_Time.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "MULTI_ID", pt.size = .01, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "FinalAssignment", pt.size = .01, split.by = 'FinalAssignment')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'Site')
    # ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate_Time.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'FinalAssignment')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time.png"), height = 4, width = 10)

    # seurat$Time_FinalAssignment <- paste0(seurat$Time, seurat$FinalAssignment)
    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time_FinalAssignment", pt.size = .01, split.by = 'Site')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time_Site.png"), height = 4, width = 10)

    # UMAPmt <- FeaturePlot(seurat, reduction = "umap", feature = "percent.mt")
    # ggsave(UMAPmt, filename = paste0(outdir,"UMAPmt.png"))

    # UMAPrb <- FeaturePlot(seurat, reduction = "umap", feature = "percent.rb")
    # ggsave(UMAPrb, filename = paste0(outdir,"UMAPrb.png"))

    umap <- DimPlot(seurat, group.by = c("phases"))  
    ggsave(umap, filename = paste0(outdir, "umap_CellCycle.png"))




} else if ( stage == "SCTnormalization_SiteMtRbRegressed"){
    # seurat <- readRDS(paste0(dir,"output/Seurat_w_Hash/QC/seurat_norm_NoOutliers.rds"))
    # seurat@meta.data$Site <- gsub("\\d","", seurat@meta.data$MULTI_ID)
    # print(head(seurat@meta.data$Site))
    # seurat <- SCTransform(seurat, vars.to.regress = c("Site","percent.mt","percent.rb"), verbose = TRUE, return.only.var.genes = FALSE)
    # seurat <- RunPCA(seurat, npcs = 100)
    # seurat <- RunUMAP(seurat, dims = 1:100, verbose = TRUE)
    # seurat <- FindNeighbors(seurat, dims = 1:100, verbose = TRUE)

    # saveRDS(seurat, paste0(outdir, "seurat_SCT_SiteMtRbRegressed.rds"))
    seurat <- readRDS(paste0(outdir, "seurat_SCT_SiteMtRbRegressed.rds"))
    
    assigned <- read.table(paste0(dir,"output/Seurat_w_Hash/CellCycle/CellCycleProportions.txt"), sep = "\t")
    seurat <- AddMetaData(seurat, assigned)


    # plots <- DimPlot(seurat, group.by = c("Time", "FinalAssignment", "Site"), combine = FALSE)
    # plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 4, 
    #     byrow = TRUE, override.aes = list(size = 2.5))))
    # ggsave(CombinePlots(plots, ncol = 3), filename = paste0(outdir,"umap.png"),width = 20, height = 9)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Pool", pt.size = .01, split.by = 'Pool')
    # ggsave(umap, filename = paste0(outdir, "umap_pool.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .01, split.by = 'Time')
    # ggsave(umap, filename = paste0(outdir, "umap_Time.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "MULTI_ID", pt.size = .01, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "FinalAssignment", pt.size = .01, split.by = 'FinalAssignment')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'Site')
    # ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate_Time.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'FinalAssignment')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time.png"), height = 4, width = 10)

    # seurat$Time_FinalAssignment <- paste0(seurat$Time, seurat$FinalAssignment)
    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time_FinalAssignment", pt.size = .01, split.by = 'Site')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time_Site.png"), height = 4, width = 10)

    # UMAPmt <- FeaturePlot(seurat, reduction = "umap", feature = "percent.mt")
    # ggsave(UMAPmt, filename = paste0(outdir,"UMAPmt.png"))

    # UMAPrb <- FeaturePlot(seurat, reduction = "umap", feature = "percent.rb")
    # ggsave(UMAPrb, filename = paste0(outdir,"UMAPrb.png"))

    umap <- DimPlot(seurat, group.by = c("phases"))  
    ggsave(umap, filename = paste0(outdir, "umap_CellCycle.png"))



} else if ( stage == "SCTnormalization_PoolRegressed" ){
    # seurat <- readRDS(paste0(dir,"output/Seurat_w_Hash/QC/seurat_norm_NoOutliers.rds"))
    # seurat <- SCTransform(seurat, vars.to.regress = "Pool", verbose = TRUE)
    # seurat <- RunPCA(seurat, npcs = 100)
    # seurat <- RunUMAP(seurat, dims = 1:100, verbose = TRUE)
    # seurat <- FindNeighbors(seurat, dims = 1:100, verbose = TRUE)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Pool", pt.size = .01, split.by = 'Pool')
    # ggsave(umap, filename = paste0(outdir, "umap_pool.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .01, split.by = 'Time')
    # ggsave(umap, filename = paste0(outdir, "umap_Time.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "MULTI_ID", pt.size = .01, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "FinalAssignment", pt.size = .01, split.by = 'FinalAssignment')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine.png"), height = 4, width = 10)

    # saveRDS(seurat, paste0(outdir, "seurat_SCT_PoolRegressed.rds"))
    seurat <- readRDS(paste0(outdir, "seurat_SCT_PoolRegressed.rds"))
    
    assigned <- read.table(paste0(dir,"output/Seurat_w_Hash/CellCycle/CellCycleProportions.txt"), sep = "\t")
    seurat <- AddMetaData(seurat, assigned)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .01, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate_Time.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .01, split.by = 'FinalAssignment')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time.png"), height = 4, width = 10)

    # seurat$Time_FinalAssignment <- paste0(seurat$Time, seurat$FinalAssignment)
    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time_FinalAssignment", pt.size = .01, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time_Site.png"), height = 4, width = 10)

    umap <- DimPlot(seurat, group.by = c("phases"))  
    ggsave(umap, filename = paste0(outdir, "umap_CellCycle.png"))


} else if ( stage == "SCTnormalization_PoolDayRegressed"){
    # seurat <- readRDS(paste0(dir,"output/Seurat_w_Hash/QC/seurat_norm_NoOutliers.rds"))
    # seurat <- SCTransform(seurat, vars.to.regress = c("Pool", "Time"), verbose = TRUE)
    # seurat <- RunPCA(seurat, npcs = 100)
    # seurat <- RunUMAP(seurat, dims = 1:100, verbose = TRUE)
    # seurat <- FindNeighbors(seurat, dims = 1:100, verbose = TRUE)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Pool", pt.size = .01, split.by = 'Pool')
    # ggsave(umap, filename = paste0(outdir, "umap_pool.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .01, split.by = 'Time')
    # ggsave(umap, filename = paste0(outdir, "umap_Time.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "MULTI_ID", pt.size = .01, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "FinalAssignment", pt.size = .01, split.by = 'FinalAssignment')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine.png"), height = 4, width = 10)

    # saveRDS(seurat, paste0(outdir, "seurat_SCT_PoolDayRegressed.rds"))
    seurat <- readRDS(paste0(outdir, "seurat_SCT_PoolDayRegressed.rds"))
    
    assigned <- read.table(paste0(dir,"output/Seurat_w_Hash/CellCycle/CellCycleProportions.txt"), sep = "\t")
    seurat <- AddMetaData(seurat, assigned)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .01, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate_Time.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .01, split.by = 'FinalAssignment')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time.png"), height = 4, width = 10)

    # seurat$Time_FinalAssignment <- paste0(seurat$Time, seurat$FinalAssignment)
    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time_FinalAssignment", pt.size = .01, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time_Site.png"), height = 4, width = 10)

    umap <- DimPlot(seurat, group.by = c("phases"))  
    ggsave(umap, filename = paste0(outdir, "umap_CellCycle.png"))




} else if (stage[1]=="Harmony_Pool"){
    # seurat <- readRDS(paste0(dir,"output/Seurat_w_Hash/QC/seurat_norm_NoOutliers.rds"))

    # seurat <- NormalizeData(seurat, verbose = FALSE) %>%
    #     FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    #     ScaleData(verbose = FALSE) %>% 
    #     RunPCA(pc.genes = seurat@var.genes, npcs = 20, verbose = FALSE)
    # seurat <- RunHarmony(seurat,"Pool")

    # ### Plot the metrics for the results ###
    # p1 <- DimPlot(object = seurat, reduction = "harmony", pt.size = .1, group.by = "Pool", do.return = TRUE)
    # ggsave(p1, filename = paste0(outdir,"harmony_PCA.png"))

    # p2 <- VlnPlot(object = seurat, features = "harmony_1", group.by = "Pool", pt.size = .1)
    # ggsave(p2, filename = paste0(outdir,"harmony_violin.png"))

    # seurat <- RunUMAP(seurat, reduction = "harmony", dims = 1:20)
    # seurat <- FindNeighbors(seurat, reduction = "harmony", dims = 1:20)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Pool", pt.size = .01, split.by = 'Pool')
    # ggsave(umap, filename = paste0(outdir, "umap_pool.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .01, split.by = 'Time')
    # ggsave(umap, filename = paste0(outdir, "umap_Time.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "MULTI_ID", pt.size = .01, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "FinalAssignment", pt.size = .01, split.by = 'FinalAssignment')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine.png"), height = 4, width = 10)

    # saveRDS(seurat, paste0(outdir,"NoOutliers_Pool_harmony.rds"))
    seurat <- readRDS(paste0(outdir, "NoOutliers_Pool_harmony.rds"))

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .01, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate_Time.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .01, split.by = 'FinalAssignment')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time.png"), height = 4, width = 10)

    # seurat$Time_FinalAssignment <- paste0(seurat$Time, seurat$FinalAssignment)
    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time_FinalAssignment", pt.size = .01, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time_Site.png"), height = 4, width = 10)

    assigned <- read.table(paste0(dir,"output/Seurat_w_Hash/CellCycle/CellCycleProportions.txt"), sep = "\t")
    seurat <- AddMetaData(seurat, assigned)

    umap <- DimPlot(seurat, group.by = c("phases"))  
    ggsave(umap, filename = paste0(outdir, "umap_CellCycle.png"))


} else if (stage[1]=="Harmony_Site"){
    # seurat <- readRDS(paste0(dir,"output/Seurat_w_Hash/QC/seurat_norm_NoOutliers.rds"))
    # seurat@meta.data$Site <- gsub("\\d","", seurat@meta.data$MULTI_ID)
    # print(head(seurat@meta.data$Site))
    # seurat <- NormalizeData(seurat, verbose = FALSE) %>%
    #     FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    #     ScaleData(verbose = FALSE) %>% 
    #     RunPCA(pc.genes = seurat@var.genes, npcs = 20, verbose = FALSE)
    # seurat <- RunHarmony(seurat,"Site")

    # ### Plot the metrics for the results ###
    # p1 <- DimPlot(object = seurat, reduction = "harmony", pt.size = .1, group.by = "Pool", do.return = TRUE)
    # ggsave(p1, filename = paste0(outdir,"harmony_PCA.png"))

    # p2 <- VlnPlot(object = seurat, features = "harmony_1", group.by = "Pool", pt.size = .1)
    # ggsave(p2, filename = paste0(outdir,"harmony_violin.png"))

    # seurat <- RunUMAP(seurat, reduction = "harmony", dims = 1:20)
    # seurat <- FindNeighbors(seurat, reduction = "harmony", dims = 1:20)


    # saveRDS(seurat, paste0(outdir,"NoOutliers_Site_harmony.rds"))
    seurat <- readRDS(paste0(outdir,"NoOutliers_Site_harmony.rds"))

    # seurat@meta.data$Site <- gsub("\\d","",seurat@meta.data$MULTI_ID)

    # plots <- DimPlot(seurat, group.by = c("Time", "FinalAssignment", "Site"), combine = FALSE)
    # plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 4, 
    #     byrow = TRUE, override.aes = list(size = 2.5))))
    # ggsave(CombinePlots(plots, ncol = 3), filename = paste0(outdir,"umap.png"),width = 20, height = 9)


    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Pool", pt.size = .01, split.by = 'Pool')
    # ggsave(umap, filename = paste0(outdir, "umap_pool.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .01, split.by = 'Time')
    # ggsave(umap, filename = paste0(outdir, "umap_Time.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "MULTI_ID", pt.size = .01, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "FinalAssignment", pt.size = .01, split.by = 'FinalAssignment')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .01, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate_Time.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .01, split.by = 'FinalAssignment')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time.png"), height = 4, width = 10)

    # seurat$Time_FinalAssignment <- paste0(seurat$Time, seurat$FinalAssignment)
    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time_FinalAssignment", pt.size = .01, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time_Site.png"), height = 4, width = 10)
    
    # UMAPmt <- FeaturePlot(seurat, reduction = "umap", feature = "percent.mt")
    # ggsave(UMAPmt, filename = paste0(outdir,"UMAPmt.png"))

    # UMAPrb <- FeaturePlot(seurat, reduction = "umap", feature = "percent.rb")
    # ggsave(UMAPrb, filename = paste0(outdir,"UMAPrb.png"))

    assigned <- read.table(paste0(dir,"output/Seurat_w_Hash/CellCycle/CellCycleProportions.txt"), sep = "\t")
    seurat <- AddMetaData(seurat, assigned)

    umap <- DimPlot(seurat, group.by = c("phases"))  
    ggsave(umap, filename = paste0(outdir, "umap_CellCycle.png"))

} else if (stage[1]=="Harmony_Site_Regress_MtRb"){
    # seurat <- readRDS(paste0(dir,"output/Seurat_w_Hash/QC/seurat_norm_NoOutliers.rds"))
    # seurat@meta.data$Site <- gsub("\\d","", seurat@meta.data$MULTI_ID)
    # print(head(seurat@meta.data$Site))
    # seurat <- NormalizeData(seurat, verbose = FALSE) %>%
    #     FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    #     ScaleData(verbose = FALSE, vars.to.regress = c("percent.mt", "percent.rb")) %>% 
    #     RunPCA(pc.genes = seurat@var.genes, npcs = 20, verbose = FALSE)
    # seurat <- RunHarmony(seurat,"Site")

    # ### Plot the metrics for the results ###
    # p1 <- DimPlot(object = seurat, reduction = "harmony", pt.size = .1, group.by = "Pool", do.return = TRUE)
    # ggsave(p1, filename = paste0(outdir,"harmony_PCA.png"))

    # p2 <- VlnPlot(object = seurat, features = "harmony_1", group.by = "Pool", pt.size = .1)
    # ggsave(p2, filename = paste0(outdir,"harmony_violin.png"))

    # seurat <- RunUMAP(seurat, reduction = "harmony", dims = 1:20)
    # seurat <- FindNeighbors(seurat, reduction = "harmony", dims = 1:20)

    # saveRDS(seurat, paste0(outdir,"NoOutliers_Site_harmony.rds"))
    seurat <- readRDS(paste0(outdir,"NoOutliers_Site_harmony.rds"))

    # seurat@meta.data$Site <- gsub("\\d","",seurat@meta.data$MULTI_ID)

    # plots <- DimPlot(seurat, group.by = c("Time", "FinalAssignment", "Site"), combine = FALSE)
    # plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 4, 
    #     byrow = TRUE, override.aes = list(size = 2.5))))
    # ggsave(CombinePlots(plots, ncol = 3), filename = paste0(outdir,"umap.png"),width = 20, height = 9)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Pool", pt.size = .01, split.by = 'Pool')
    # ggsave(umap, filename = paste0(outdir, "umap_pool.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .01, split.by = 'Time')
    # ggsave(umap, filename = paste0(outdir, "umap_Time.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "MULTI_ID", pt.size = .01, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "FinalAssignment", pt.size = .01, split.by = 'FinalAssignment')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .01, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate_Time.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .01, split.by = 'FinalAssignment')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time.png"), height = 4, width = 10)

    # seurat$Time_FinalAssignment <- paste0(seurat$Time, seurat$FinalAssignment)
    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time_FinalAssignment", pt.size = .01, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time_Site.png"), height = 4, width = 10)
    
    # UMAPmt <- FeaturePlot(seurat, reduction = "umap", feature = "percent.mt")
    # ggsave(UMAPmt, filename = paste0(outdir,"UMAPmt.png"))

    # UMAPrb <- FeaturePlot(seurat, reduction = "umap", feature = "percent.rb")
    # ggsave(UMAPrb, filename = paste0(outdir,"UMAPrb.png"))

    assigned <- read.table(paste0(dir,"output/Seurat_w_Hash/CellCycle/CellCycleProportions.txt"), sep = "\t")
    seurat <- AddMetaData(seurat, assigned)

    umap <- DimPlot(seurat, group.by = c("phases"))  
    ggsave(umap, filename = paste0(outdir, "umap_CellCycle.png"))




} else if (stage[1]=="Harmony_Day"){
    # seurat <- readRDS(paste0(dir,"output/Seurat_w_Hash/QC/seurat_norm_NoOutliers.rds"))

    # seurat <- NormalizeData(seurat, verbose = FALSE) %>%
    #     FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    #     ScaleData(verbose = FALSE) %>% 
    #     RunPCA(pc.genes = seurat@var.genes, npcs = 20, verbose = FALSE)
    # seurat <- RunHarmony(seurat,"Time", plot_convergence = TRUE)

    # ### Plot the metrics for the results ###
    # p1 <- DimPlot(object = seurat, reduction = "harmony", pt.size = .1, group.by = "Time")
    # ggsave(p1, filename = paste0(outdir,"harmony_PCA.png"))

    # p2 <- VlnPlot(object = seurat, features = "harmony_1", group.by = "Time", pt.size = .1)
    # ggsave(p2, filename = paste0(outdir,"harmony_violin.png"))


    # seurat <- RunUMAP(seurat, reduction = "harmony", dims = 1:20)
    # seurat <- FindNeighbors(seurat, reduction = "harmony", dims = 1:20)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Pool", pt.size = .01, split.by = 'Pool')
    # ggsave(umap, filename = paste0(outdir, "umap_pool.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .01, split.by = 'Time')
    # ggsave(umap, filename = paste0(outdir, "umap_Time.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "MULTI_ID", pt.size = .01, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "FinalAssignment", pt.size = .01, split.by = 'FinalAssignment')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine.png"), height = 4, width = 10)

    # saveRDS(seurat, paste0(outdir,"NoOutliers_Day_harmony.rds"))
    seurat <- readRDS(paste0(outdir, "NoOutliers_Day_harmony.rds"))

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .01, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate_Time.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .01, split.by = 'FinalAssignment')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time.png"), height = 4, width = 10)

    # seurat$Time_FinalAssignment <- paste0(seurat$Time, seurat$FinalAssignment)
    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time_FinalAssignment", pt.size = .01, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time_Site.png"), height = 4, width = 10)

    assigned <- read.table(paste0(dir,"output/Seurat_w_Hash/CellCycle/CellCycleProportions.txt"), sep = "\t")
    seurat <- AddMetaData(seurat, assigned)

    umap <- DimPlot(seurat, group.by = c("phases"))  
    ggsave(umap, filename = paste0(outdir, "umap_CellCycle.png"))



} else if (stage[1]=="Harmony_DayPatient"){
    # seurat <- readRDS(paste0(dir,"output/Seurat_w_Hash/QC/seurat_norm_NoOutliers.rds"))

    # seurat <- SCTransform(seurat, vars.to.regress = "Pool", verbose = TRUE)
    # seurat <- RunPCA(seurat, npcs = 100)
    # seurat <- RunHarmony(seurat,c("Time", "FinalAssignment"), plot_convergence = TRUE, assay.use = "SCT")

    # ### Plot the metrics for the results ###
    # p1 <- DimPlot(object = seurat, reduction = "harmony", pt.size = .1, group.by = "Time")
    # ggsave(p1, filename = paste0(outdir,"harmony_PCA.png"))

    # p2 <- VlnPlot(object = seurat, features = "harmony_1", group.by = "Time", pt.size = .1)
    # ggsave(p2, filename = paste0(outdir,"harmony_violin.png"))

    # seurat <- RunUMAP(seurat, reduction = "harmony", dims = 1:20)
    # seurat <- FindNeighbors(seurat, reduction = "harmony", dims = 1:20)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Pool", pt.size = .01, split.by = 'Pool')
    # ggsave(umap, filename = paste0(outdir, "umap_pool.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .01, split.by = 'Time')
    # ggsave(umap, filename = paste0(outdir, "umap_Time.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "MULTI_ID", pt.size = .01, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "FinalAssignment", pt.size = .01, split.by = 'FinalAssignment')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine.png"), height = 4, width = 10)

    # saveRDS(seurat, paste0(outdir,"NoOutliers_DayPatient_harmony.rds"))
    seurat <- readRDS(paste0(outdir, "NoOutliers_DayPatient_harmony.rds"))

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate_Time.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'FinalAssignment')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time.png"), height = 4, width = 10)

    # seurat$Time_FinalAssignment <- paste0(seurat$Time, seurat$FinalAssignment)
    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time_FinalAssignment", pt.size = .001, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time_Site.png"), height = 4, width = 10)

    assigned <- read.table(paste0(dir,"output/Seurat_w_Hash/CellCycle/CellCycleProportions.txt"), sep = "\t")
    seurat <- AddMetaData(seurat, assigned)

    umap <- DimPlot(seurat, group.by = c("phases"))  
    ggsave(umap, filename = paste0(outdir, "umap_CellCycle.png"))



} else if (stage[1]=="Integration_DayPatient"){
    # seurat <- readRDS(paste0(dir,"output/Seurat_w_Hash/QC/seurat_norm_NoOutliers.rds"))
    # seurat@meta.data$Day_Patient <- paste0(seurat@meta.data$Time, "_", seurat@meta.data$FinalAssignment)
    # print(unique( seurat@meta.data$Day_Patient))

    # seurat_list <- SplitObject(seurat, split.by = "Day_Patient")

    # seurat_list <- lapply(seurat_list, function(x) {
    #     SCTransform(x, verbose = TRUE)
    # })

    # seurat_features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)
    # seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = seurat_features)
    # seurat_anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT", 
    #     anchor.features = seurat_features)
    # seurat_integrated <- IntegrateData(anchorset = seurat_anchors, normalization.method = "SCT")

    # seurat_integrated <- RunPCA(object = seurat_integrated, verbose = FALSE)
    # seurat_integrated <- RunUMAP(object = seurat_integrated, dims = 1:30)
    # seurat_integrated <- FindNeighbors(seurat_integrated, reduction = "pca", assay = "integrated", dims = 1:20)


    # plots <- DimPlot(seurat_integrated, group.by = c("Time", "FinalAssignment"), combine = FALSE)
    # plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 4, 
    #     byrow = TRUE, override.aes = list(size = 2.5))))
    # ggsave(CombinePlots(plots), filename = paste0(outdir,"umap.png"),width = 12, height = 9)

    # umap <- DimPlot(seurat_integrated, reduction = "umap", group.by = "Pool", pt.size = .01, split.by = 'Pool')
    # ggsave(umap, filename = paste0(outdir, "umap_pool.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat_integrated, reduction = "umap", group.by = "Day_Patient", pt.size = .01, split.by = 'Day_Patient')
    # ggsave(umap, filename = paste0(outdir, "umap_DayPatient.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat_integrated, reduction = "umap", group.by = "Time", pt.size = .01, split.by = 'Time')
    # ggsave(umap, filename = paste0(outdir, "umap_Time.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat_integrated, reduction = "umap", group.by = "MULTI_ID", pt.size = .01, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat_integrated, reduction = "umap", group.by = "FinalAssignment", pt.size = .01, split.by = 'FinalAssignment')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine.png"), height = 4, width = 10)

    # saveRDS(seurat_integrated, paste0(outdir,"NoOutliers_DayPatient_SeuratIntegratee.rds"))
    seurat <- readRDS(paste0(outdir, "NoOutliers_DayPatient_SeuratIntegratee.rds"))

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate_Time.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'FinalAssignment')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time.png"), height = 4, width = 10)

    # seurat$Time_FinalAssignment <- paste0(seurat$Time, seurat$FinalAssignment)
    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time_FinalAssignment", pt.size = .001, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time_Site.png"), height = 4, width = 10)

    assigned <- read.table(paste0(dir,"output/Seurat_w_Hash/CellCycle/CellCycleProportions.txt"), sep = "\t")
    seurat <- AddMetaData(seurat, assigned)

    umap <- DimPlot(seurat, group.by = c("phases"))  
    ggsave(umap, filename = paste0(outdir, "umap_CellCycle.png"))



} else if (stage[1]=="Integration_Day"){
    # seurat <- readRDS(paste0(dir,"output/Seurat_w_Hash/QC/seurat_norm_NoOutliers.rds"))
    # seurat_list <- SplitObject(seurat, split.by = "Time")

    # seurat_list <- lapply(seurat_list, function(x) {
    #     SCTransform(x, verbose = TRUE)
    # })

    # seurat_features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)
    # seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = seurat_features)
    # seurat_anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT", 
    #     anchor.features = seurat_features)
    # seurat_integrated <- IntegrateData(anchorset = seurat_anchors, normalization.method = "SCT")

    # seurat_integrated <- RunPCA(object = seurat_integrated, verbose = FALSE)
    # seurat_integrated <- RunUMAP(object = seurat_integrated, dims = 1:30)

    # seurat_integrated <- FindNeighbors(seurat_integrated, reduction = "pca", assay = "integrated", dims = 1:20)


    # plots <- DimPlot(seurat_integrated, group.by = c("Time", "FinalAssignment"), combine = FALSE)
    # plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 4, 
    #     byrow = TRUE, override.aes = list(size = 2.5))))
    # ggsave(CombinePlots(plots), filename = paste0(outdir,"umap.png"),width = 12, height = 9)

    #     umap <- DimPlot(seurat_integrated, reduction = "umap", group.by = "Pool", pt.size = .01, split.by = 'Pool')
    # ggsave(umap, filename = paste0(outdir, "umap_pool.png"), height = 4, width = 10)


    # umap <- DimPlot(seurat_integrated, reduction = "umap", group.by = "Time", pt.size = .01, split.by = 'Time')
    # ggsave(umap, filename = paste0(outdir, "umap_Time.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat_integrated, reduction = "umap", group.by = "MULTI_ID", pt.size = .01, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "FinalAssignment", pt.size = .01, split.by = 'FinalAssignment')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine.png"), height = 4, width = 10)

    # saveRDS(seurat_integrated, paste0(outdir,"NoOutliers_Day_SeuratIntegratee.rds"))
    seurat <- readRDS(paste0(outdir, "NoOutliers_Day_SeuratIntegratee.rds"))

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate_Time.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'FinalAssignment')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time.png"), height = 4, width = 10)

    # seurat$Time_FinalAssignment <- paste0(seurat$Time, seurat$FinalAssignment)
    # umap <- DimPlot(seurat, reduction = "umap", group.by = "Time_FinalAssignment", pt.size = .001, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time_Site.png"), height = 4, width = 10)

    assigned <- read.table(paste0(dir,"output/Seurat_w_Hash/CellCycle/CellCycleProportions.txt"), sep = "\t")
    seurat <- AddMetaData(seurat, assigned)

    umap <- DimPlot(seurat, group.by = c("phases"))  
    ggsave(umap, filename = paste0(outdir, "umap_CellCycle.png"))




} else if (stage[1]=="Integration_Site"){
    # seurat <- readRDS(paste0(dir,"output/Seurat_w_Hash/QC/seurat_norm_NoOutliers.rds"))
    # seurat@meta.data$Site <- gsub("\\d","", seurat@meta.data$MULTI_ID)
    # features <- rownames(seurat[["RNA"]][Matrix::rowSums(seurat[["RNA"]]@counts > 0) > (0.01*nrow(seurat))])

    # print(head(seurat@meta.data$Site))
    # seurat_list <- SplitObject(seurat, split.by = "Site")

    # seurat_list <- lapply(seurat_list, function(x) {
    #     SCTransform(x, verbose = TRUE)
    # })

    # seurat_features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)
    # seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = seurat_features)
    # seurat_anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT", 
    #     anchor.features = seurat_features)
    # seurat_integrated <- IntegrateData(anchorset = seurat_anchors, normalization.method = "SCT",features.to.integrate = features)

    # seurat_integrated <- RunPCA(object = seurat_integrated, verbose = FALSE)
    # seurat_integrated <- RunUMAP(object = seurat_integrated, dims = 1:30)

    # seurat_integrated <- FindNeighbors(seurat_integrated, reduction = "pca", assay = "integrated", dims = 1:20)

    # saveRDS(seurat_integrated,paste0(outdir, "NoOutliers_Site_SeuratIntegratee.rds"))
    seurat_integrated <- readRDS(paste0(outdir, "NoOutliers_Site_SeuratIntegratee.rds"))


    # plots <- DimPlot(seurat_integrated, group.by = c("Time", "FinalAssignment", "Site"), combine = FALSE)
    # plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 4, 
    #     byrow = TRUE, override.aes = list(size = 2.5))))
    # ggsave(CombinePlots(plots, ncol = 3), filename = paste0(outdir,"umap.png"),width = 20, height = 9)

    # umap <- DimPlot(seurat_integrated, reduction = "umap", group.by = "Pool", pt.size = .01, split.by = 'Pool')
    # ggsave(umap, filename = paste0(outdir, "umap_pool.png"), height = 4, width = 10)


    # umap <- DimPlot(seurat_integrated, reduction = "umap", group.by = "Time", pt.size = .01, split.by = 'Time')
    # ggsave(umap, filename = paste0(outdir, "umap_Time.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat_integrated, reduction = "umap", group.by = "MULTI_ID", pt.size = .01, split.by = 'Site')
    # ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat_integrated, reduction = "umap", group.by = "FinalAssignment", pt.size = .01, split.by = 'FinalAssignment')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat_integrated, reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'Site')
    # ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate_Time.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat_integrated, reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'FinalAssignment')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time.png"), height = 4, width = 10)

    # seurat_integrated$Time_FinalAssignment <- paste0(seurat_integrated$Time, seurat_integrated$FinalAssignment)
    # umap <- DimPlot(seurat_integrated, reduction = "umap", group.by = "Time_FinalAssignment", pt.size = .001, split.by = 'Site')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time_Site.png"), height = 4, width = 10)

    # UMAPmt <- FeaturePlot(seurat_integrated, reduction = "umap", feature = "percent.mt")
    # ggsave(UMAPmt, filename = paste0(outdir,"UMAPmt.png"))

    # UMAPrb <- FeaturePlot(seurat_integrated, reduction = "umap", feature = "percent.rb")
    # ggsave(UMAPrb, filename = paste0(outdir,"UMAPrb.png"))

    assigned <- read.table(paste0(dir,"output/Seurat_w_Hash/CellCycle/CellCycleProportions.txt"), sep = "\t")
    seurat_integrated <- AddMetaData(seurat_integrated, assigned)

    umap <- DimPlot(seurat_integrated, group.by = c("phases"))  
    ggsave(umap, filename = paste0(outdir, "umap_CellCycle.png"))

} else if (stage[1]=="Integration_Site_Regress_MtRb"){
    # seurat <- readRDS(paste0(dir,"output/Seurat_w_Hash/QC/seurat_norm_NoOutliers.rds"))
    # seurat@meta.data$Site <- gsub("\\d","", seurat@meta.data$MULTI_ID)
    # print(head(seurat@meta.data$Site))
    # seurat_list <- SplitObject(seurat, split.by = "Site")

    # seurat_list <- lapply(seurat_list, function(x) {
    #     SCTransform(x, verbose = TRUE, vars.to.regress = c("percent.mt", "percent.rb"))
    # })

    # seurat_features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)
    # seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = seurat_features)
    # seurat_anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT", 
    #     anchor.features = seurat_features)
    # seurat_integrated <- IntegrateData(anchorset = seurat_anchors, normalization.method = "SCT",features.to.integrate = rownames(seurat[["RNA"]][Matrix::rowSums(seurat[["RNA"]]@counts > 0) > (0.01*nrow(seurat))]))

    # seurat_integrated <- RunPCA(object = seurat_integrated, verbose = FALSE)
    # seurat_integrated <- RunUMAP(object = seurat_integrated, dims = 1:30)

    # seurat_integrated <- FindNeighbors(seurat_integrated, reduction = "pca", assay = "integrated", dims = 1:20)

    # saveRDS(seurat_integrated,paste0(outdir, "NoOutliers_Site_SeuratIntegratee.rds"))


    seurat_integrated <- readRDS(paste0(outdir, "NoOutliers_Site_SeuratIntegratee.rds"))


    # plots <- DimPlot(seurat_integrated, group.by = c("Time", "FinalAssignment", "Site"), combine = FALSE)
    # plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 4, 
    #     byrow = TRUE, override.aes = list(size = 2.5))))
    # ggsave(CombinePlots(plots, ncol = 3), filename = paste0(outdir,"umap.png"),width = 20, height = 9)

    # umap <- DimPlot(seurat_integrated, reduction = "umap", group.by = "Pool", pt.size = .01, split.by = 'Pool')
    # ggsave(umap, filename = paste0(outdir, "umap_pool.png"), height = 4, width = 10)


    # umap <- DimPlot(seurat_integrated, reduction = "umap", group.by = "Time", pt.size = .01, split.by = 'Time')
    # ggsave(umap, filename = paste0(outdir, "umap_Time.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat_integrated, reduction = "umap", group.by = "MULTI_ID", pt.size = .01, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat_integrated, reduction = "umap", group.by = "FinalAssignment", pt.size = .01, split.by = 'FinalAssignment')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat_integrated, reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate_Time.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat_integrated, reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'FinalAssignment')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time.png"), height = 4, width = 10)

    # seurat_integrated$Time_FinalAssignment <- paste0(seurat_integrated$Time, seurat_integrated$FinalAssignment)
    # umap <- DimPlot(seurat_integrated, reduction = "umap", group.by = "Time_FinalAssignment", pt.size = .001, split.by = 'Site')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time_Site.png"), height = 4, width = 10)

    # UMAPmt <- FeaturePlot(seurat_integrated, reduction = "umap", feature = "percent.mt")
    # ggsave(UMAPmt, filename = paste0(outdir,"UMAPmt.png"))

    # UMAPrb <- FeaturePlot(seurat_integrated, reduction = "umap", feature = "percent.rb")
    # ggsave(UMAPrb, filename = paste0(outdir,"UMAPrb.png"))

    assigned <- read.table(paste0(dir,"output/Seurat_w_Hash/CellCycle/CellCycleProportions.txt"), sep = "\t")
    seurat_integrated <- AddMetaData(seurat_integrated, assigned)

    umap <- DimPlot(seurat_integrated, group.by = c("phases"))  
    ggsave(umap, filename = paste0(outdir, "umap_CellCycle.png"))

} else if (stage[1]=="Integration_Day_then_Individual"){
    # seurat <- readRDS(paste0(dir,"output/Seurat_w_Hash/QC/seurat_norm_NoOutliers.rds"))

    # seurat_list <- SplitObject(seurat, split.by = "Time")

    # seurat_list <- lapply(seurat_list, function(x) {
    #     SCTransform(x, verbose = TRUE)
    # })

    # seurat_features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)
    # seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = seurat_features)
    # seurat_anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT", 
    #     anchor.features = seurat_features)
    # seurat_integrated <- IntegrateData(anchorset = seurat_anchors, normalization.method = "SCT")

    # ### Second round of integration ###
    # seurat_integrated_list <- SplitObject(seurat_integrated, split.by = "FinalAssignment")
    # seurat_integrated_list <- lapply(seurat_integrated_list, function(x) {
    #     SCTransform(x, verbose = TRUE)
    # })

    # seurat_features <- SelectIntegrationFeatures(object.list = seurat_integrated_list, nfeatures = 3000, assay = c("integrated","integrated","integrated"))
    # seurat_integrated_list <- PrepSCTIntegration(object.list = seurat_integrated_list, anchor.features = seurat_features, assay = "integrated")
    # seurat_anchors <- FindIntegrationAnchors(object.list = seurat_integrated_list, normalization.method = "SCT", 
    #     anchor.features = seurat_features, assay = c("integrated","integrated","integrated"))
    # seurat_integrated_2 <- IntegrateData(anchorset = seurat_anchors, normalization.method = "SCT", new.assay.name = "integrated_2")


    # seurat_integrated_2 <- RunPCA(object = seurat_integrated_2, verbose = FALSE, assay = "integrated_2")
    # seurat_integrated_2 <- RunUMAP(object = seurat_integrated_2, dims = 1:30, assay = "integrated_2")


    # seurat_integrated_2 <- FindNeighbors(seurat_integrated_2, reduction = "pca", assay = "integrated_2", dims = 1:20)


    # plots <- DimPlot(seurat_integrated_2, group.by = c("Time", "FinalAssignment"), combine = FALSE)
    # plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 4, 
    #     byrow = TRUE, override.aes = list(size = 2.5))))
    # ggsave(CombinePlots(plots), filename = paste0(outdir,"umap.png"),width = 12, height = 9)

    # umap <- DimPlot(seurat_integrated_2, reduction = "umap", group.by = "Pool", pt.size = .01, split.by = 'Pool')
    # ggsave(umap, filename = paste0(outdir, "umap_pool.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat_integrated_2, reduction = "umap", group.by = "Time", pt.size = .01, split.by = 'Time')
    # ggsave(umap, filename = paste0(outdir, "umap_Time.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat_integrated_2, reduction = "umap", group.by = "MULTI_ID", pt.size = .01, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat_integrated_2, reduction = "umap", group.by = "FinalAssignment", pt.size = .01, split.by = 'FinalAssignment')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat_integrated_2, reduction = "umap", group.by = "MULTI_ID", pt.size = .01, split.by = 'Time')
    # ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate.png"), height = 4, width = 10)

    # saveRDS(seurat_integrated_2, paste0(outdir,"NoOutliers_Day_SeuratIntegratee.rds"))
    seurat_integrated_2 <- readRDS(paste0(outdir, "NoOutliers_Day_SeuratIntegratee.rds"))

    # umap <- DimPlot(seurat_integrated_2, reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_SiteSampleReplicate_Time.png"), height = 4, width = 10)

    # umap <- DimPlot(seurat_integrated_2, reduction = "umap", group.by = "Time", pt.size = .001, split.by = 'FinalAssignment')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time.png"), height = 4, width = 10)

    # seurat_integrated_2$Time_FinalAssignment <- paste0(seurat_integrated_2$Time, seurat_integrated_2$FinalAssignment)
    # umap <- DimPlot(seurat_integrated_2, reduction = "umap", group.by = "Time_FinalAssignment", pt.size = .001, split.by = 'MULTI_ID')
    # ggsave(umap, filename = paste0(outdir, "umap_CellLine_Time_Site.png"), height = 4, width = 10)

    assigned <- read.table(paste0(dir,"output/Seurat_w_Hash/CellCycle/CellCycleProportions.txt"), sep = "\t")
    seurat_integrated_2 <- AddMetaData(seurat_integrated_2, assigned)

    umap <- DimPlot(seurat_integrated_2, group.by = c("phases"))  
    ggsave(umap, filename = paste0(outdir, "umap_CellCycle.png"))




} else if (stage[1]=="Multi_Resolution_no_covatiates"){
    seurat <- readRDS(paste0(dir,"output/Seurat/SCTnormalization_no_covatiates/seurat_SCT_noCov.rds"))

    ### Find the clusters for this resolution ###
    res <- sge*0.1-0.1
    print(paste0("Finding Clusters for resolution=",res))
    seurat <- FindClusters(seurat, resolution = res)

    ### Plot the UMAP for this resolution ###
    UMAP <- DimPlot(seurat, reduction = "umap")
    ggsave(file = paste0(outdir,"UMAPresolution",res,".png"), plot = UMAP)

    ### Save the identities for this resolution ###
    identities <- as.data.frame(Idents(seurat))
    colnames(identities) <- c(paste0("Resolution_",res))

    write_delim(identities, paste0(outdir,"Identities_Resolution_",res,".txt"), delim = "\t")

} else if (stage[1]=="Multi_Resolution_PoolRegressed"){
    seurat <- readRDS(paste0(dir,"output/Seurat/SCTnormalization_PoolRegressed/seurat_SCT_PoolRegressed.rds"))

    ### Find the clusters for this resolution ###
    res <- sge*0.1-0.1
    print(paste0("Finding Clusters for resolution=",res))
    seurat <- FindClusters(seurat, resolution = res)

    ### Plot the UMAP for this resolution ###
    UMAP <- DimPlot(seurat, reduction = "umap")
    ggsave(file = paste0(outdir,"UMAPresolution",res,".png"), plot = UMAP)

    ### Save the identities for this resolution ###
    identities <- as.data.frame(Idents(seurat))
    colnames(identities) <- c(paste0("Resolution_",res))

    write_delim(identities, paste0(outdir,"Identities_Resolution_",res,".txt"), delim = "\t")


} else if (stage[1]=="Multi_Resolution_PoolDayRegressed"){
    seurat <- readRDS(paste0(dir,"output/Seurat/SCTnormalization_PoolDayRegressed/seurat_SCT_PoolDayRegressed.rds"))

    ### Find the clusters for this resolution ###
    res <- sge*0.1-0.1
    print(paste0("Finding Clusters for resolution=",res))
    seurat <- FindClusters(seurat, resolution = res)

    ### Plot the UMAP for this resolution ###
    UMAP <- DimPlot(seurat, reduction = "umap")
    ggsave(file = paste0(outdir,"UMAPresolution",res,".png"), plot = UMAP)

    ### Save the identities for this resolution ###
    identities <- as.data.frame(Idents(seurat))
    colnames(identities) <- c(paste0("Resolution_",res))

    write_delim(identities, paste0(outdir,"Identities_Resolution_",res,".txt"), delim = "\t")


} else if (stage[1]=="Multi_Resolution_Harmony_Pool"){
    seurat <- readRDS(paste0(dir,"output/Seurat/Harmony_Pool/NoOutliers_Pool_harmony.rds"))

    ### Find the clusters for this resolution ###
    res <- sge*0.01-0.01
    print(paste0("Finding Clusters for resolution=",res))
    seurat <- FindClusters(seurat, resolution = res)

    ### Plot the UMAP for this resolution ###
    UMAP <- DimPlot(seurat, reduction = "umap")
    ggsave(file = paste0(outdir,"UMAPresolution",res,".png"), plot = UMAP)

    ### Save the identities for this resolution ###
    identities <- as.data.frame(Idents(seurat))
    colnames(identities) <- c(paste0("Resolution_",res))

    write_delim(identities, paste0(outdir,"Identities_Resolution_",res,".txt"), delim = "\t")


} else if (stage[1]=="Multi_Resolution_Harmony_Day"){
    seurat <- readRDS(paste0(dir,"output/Seurat/Harmony_Day/NoOutliers_Day_harmony.rds"))

    ### Find the clusters for this resolution ###
    res <- sge*0.01-0.01
    print(paste0("Finding Clusters for resolution=",res))
    seurat <- FindClusters(seurat, resolution = res)

    ### Plot the UMAP for this resolution ###
    UMAP <- DimPlot(seurat, reduction = "umap")
    ggsave(file = paste0(outdir,"UMAPresolution",res,".png"), plot = UMAP)

    ### Save the identities for this resolution ###
    identities <- as.data.frame(Idents(seurat))
    colnames(identities) <- c(paste0("Resolution_",res))

    write_delim(identities, paste0(outdir,"Identities_Resolution_",res,".txt"), delim = "\t")



} else if (stage[1]=="Multi_Resolution_Harmony_DayPatient"){
    seurat <- readRDS(paste0(dir,"output/Seurat/Harmony_DayPatient/NoOutliers_DayPatient_harmony.rds"))

    ### Find the clusters for this resolution ###
    res <- sge*0.01-0.01
    print(paste0("Finding Clusters for resolution=",res))
    seurat <- FindClusters(seurat, resolution = res)

    ### Plot the UMAP for this resolution ###
    UMAP <- DimPlot(seurat, reduction = "umap")
    ggsave(file = paste0(outdir,"UMAPresolution",res,".png"), plot = UMAP)

    ### Save the identities for this resolution ###
    identities <- as.data.frame(Idents(seurat))
    colnames(identities) <- c(paste0("Resolution_",res))

    write_delim(identities, paste0(outdir,"Identities_Resolution_",res,".txt"), delim = "\t")


} else if (stage[1]=="Multi_Resolution_Integration_Day_then_Individual"){
    seurat <- readRDS(paste0(dir,"output/Seurat/Integration_Day_then_Individual/NoOutliers_Day_SeuratIntegratee.rds"))

    ### Find the clusters for this resolution ###
    res <- sge*0.01-0.01
    print(paste0("Finding Clusters for resolution=",res))
    seurat <- FindClusters(seurat, resolution = res)

    ### Plot the UMAP for this resolution ###
    UMAP <- DimPlot(seurat, reduction = "umap")
    ggsave(file = paste0(outdir,"UMAPresolution",res,".png"), plot = UMAP)

    ### Save the identities for this resolution ###
    identities <- as.data.frame(Idents(seurat))
    colnames(identities) <- c(paste0("Resolution_",res))

    write_delim(identities, paste0(outdir,"Identities_Resolution_",res,".txt"), delim = "\t")



} else if (stage[1]=="Clustree_no_covatiates"){
    ResolutionsFileList <- list.files(paste0(dir,"output/Seurat/Multi_Resolution_no_covatiates/"), pattern = "Identities_Resolution_")

    ResolutionsList <- lapply(ResolutionsFileList, FUN = function(x){
        read_delim(paste0(dir,"output/Seurat/Multi_Resolution_no_covatiates/",x), delim = "\t")
    })
    
    ResolutionsDF <- do.call(cbind, ResolutionsList)

    ##### CLUSTREE #####
    print("Using Clustree to diagram changes with increasing resolutions")
    pClustree <- clustree(ResolutionsDF, prefix = "Resolution_")
    ggsave(file = paste0(outdir, "ClusteringClustree.png"), plot = pClustree)


} else if (stage[1]=="Clustree_PoolRegressed"){
    ResolutionsFileList <- list.files(paste0(dir,"output/Seurat/Multi_Resolution_PoolRegressed/"), pattern = "Identities_Resolution_")

    ResolutionsList <- lapply(ResolutionsFileList, FUN = function(x){
        read_delim(paste0(dir,"output/Seurat/Multi_Resolution_PoolRegressed/",x), delim = "\t")
    })
    
    ResolutionsDF <- do.call(cbind, ResolutionsList)

    ##### CLUSTREE #####
    print("Using Clustree to diagram changes with increasing resolutions")
    pClustree <- clustree(ResolutionsDF, prefix = "Resolution_")
    ggsave(file = paste0(outdir, "ClusteringClustree.png"), plot = pClustree)


} else if (stage[1]=="Clustree_PoolDayRegressed"){
    ResolutionsFileList <- list.files(paste0(dir,"output/Seurat/Multi_Resolution_PoolDayRegressed/"), pattern = "Identities_Resolution_")

    ResolutionsList <- lapply(ResolutionsFileList, FUN = function(x){
        read_delim(paste0(dir,"output/Seurat/Multi_Resolution_PoolDayRegressed/",x), delim = "\t")
    })
    
    ResolutionsDF <- do.call(cbind, ResolutionsList)

    ##### CLUSTREE #####
    print("Using Clustree to diagram changes with increasing resolutions")
    pClustree <- clustree(ResolutionsDF, prefix = "Resolution_")
    ggsave(file = paste0(outdir, "ClusteringClustree.png"), plot = pClustree)


} else if (stage[1]=="Clustree_Harmony_Pool"){
    ResolutionsFileList <- list.files(paste0(dir,"output/Seurat/Multi_Resolution_Harmony_Pool/"), pattern = "Identities_Resolution_")

    ResolutionsList <- lapply(ResolutionsFileList, FUN = function(x){
        read_delim(paste0(dir,"output/Seurat/Multi_Resolution_Harmony_Pool/",x), delim = "\t")
    })
    
    ResolutionsDF <- do.call(cbind, ResolutionsList)

    ##### CLUSTREE #####
    print("Using Clustree to diagram changes with increasing resolutions")
    pClustree <- clustree(ResolutionsDF, prefix = "Resolution_")
    ggsave(file = paste0(outdir, "ClusteringClustree.png"), plot = pClustree)


} else if (stage[1]=="Clustree_Harmony_Day"){
    ResolutionsFileList <- list.files(paste0(dir,"output/Seurat/Multi_Resolution_Harmony_Day/"), pattern = "Identities_Resolution_")

    ResolutionsList <- lapply(ResolutionsFileList, FUN = function(x){
        read_delim(paste0(dir,"output/Seurat/Multi_Resolution_Harmony_Day/",x), delim = "\t")
    })
    
    ResolutionsDF <- do.call(cbind, ResolutionsList)

    ##### CLUSTREE #####
    print("Using Clustree to diagram changes with increasing resolutions")
    pClustree <- clustree(ResolutionsDF, prefix = "Resolution_")
    ggsave(file = paste0(outdir, "ClusteringClustree.png"), plot = pClustree)


} else if (stage[1]=="Clustree_Harmony_DayPatient"){
    ResolutionsFileList <- list.files(paste0(dir,"output/Seurat/Multi_Resolution_Harmony_DayPatient/"), pattern = "Identities_Resolution_")

    ResolutionsList <- lapply(ResolutionsFileList, FUN = function(x){
        read_delim(paste0(dir,"output/Seurat/Multi_Resolution_Harmony_DayPatient/",x), delim = "\t")
    })
    
    ResolutionsDF <- do.call(cbind, ResolutionsList)

    ##### CLUSTREE #####
    print("Using Clustree to diagram changes with increasing resolutions")
    pClustree <- clustree(ResolutionsDF, prefix = "Resolution_")
    ggsave(file = paste0(outdir, "ClusteringClustree.png"), plot = pClustree)



} else if (stage[1]=="Clustree_Integration_DayPatient"){
    ResolutionsFileList <- list.files(paste0(dir,"output/Seurat/Multi_Resolution_Integration_DayPatient/"), pattern = "Identities_Resolution_")

    ResolutionsList <- lapply(ResolutionsFileList, FUN = function(x){
        read_delim(paste0(dir,"output/Seurat/Multi_Resolution_Integration_DayPatient/",x), delim = "\t")
    })
    
    ResolutionsDF <- do.call(cbind, ResolutionsList)

    ##### CLUSTREE #####
    print("Using Clustree to diagram changes with increasing resolutions")
    pClustree <- clustree(ResolutionsDF, prefix = "Resolution_")
    ggsave(file = paste0(outdir, "ClusteringClustree.png"), plot = pClustree)


} else if (stage[1]=="Clustree_Integration_Day_then_Individual"){
    ResolutionsFileList <- list.files(paste0(dir,"output/Seurat/Multi_Resolution_Integration_Day_then_Individual/"), pattern = "Identities_Resolution_")

    ResolutionsList <- lapply(ResolutionsFileList, FUN = function(x){
        read_delim(paste0(dir,"output/Seurat/Multi_Resolution_Integration_Day_then_Individual/",x), delim = "\t")
    })
    
    ResolutionsDF <- do.call(cbind, ResolutionsList)

    ##### CLUSTREE #####
    print("Using Clustree to diagram changes with increasing resolutions")
    pClustree <- clustree(ResolutionsDF, prefix = "Resolution_")
    ggsave(file = paste0(outdir, "ClusteringClustree.png"), plot = pClustree)




} else if (stage[1]=="QCfigurew_no_covatiates"){
    seurat <- readRDS(paste0(dir,"output/Seurat/SCTnormalization_no_covatiates/seurat_SCT_noCov.rds"))
    UMAPpool <- DimPlot(seurat, reduction = "umap", group.by = "Pool")
    ggsave(UMAPpool, filename = paste0(outdir,"UMAPpool.png"))
    UMAPmt <- FeaturePlot(seurat, reduction = "umap", feature = "percent.mt")
    ggsave(UMAPmt, filename = paste0(outdir,"UMAPmt.png"))
    UMAPrb <- FeaturePlot(seurat, reduction = "umap", feature = "percent.rb")
    ggsave(UMAPrb, filename = paste0(outdir,"UMAPrb.png"))
    UMAPlib <- FeaturePlot(seurat, reduction = "umap", feature = "nCount_RNA")
    ggsave(UMAPlib, filename = paste0(outdir,"UMAPlib.png"))
    UMAPgenes <- FeaturePlot(seurat, reduction = "umap", feature = "nFeature_RNA")
    ggsave(UMAPgenes, filename = paste0(outdir,"UMAPgenes.png"))
    UMAPlibSCT <- FeaturePlot(seurat, reduction = "umap", feature = "nCount_SCT")
    ggsave(UMAPlibSCT, filename = paste0(outdir,"UMAPlibSCT.png"))
    UMAPgenesSCT <- FeaturePlot(seurat, reduction = "umap", feature = "nFeature_SCT")
    ggsave(UMAPgenesSCT, filename = paste0(outdir,"UMAPgenesSCT.png"))
    UMAPdonor_souporcell <- DimPlot(seurat, reduction = "umap", group.by = "common_assignment_souporcell")
    ggsave(UMAPdonor_souporcell, filename = paste0(outdir,"UMAPdonor_souporcell.png"))
    UMAPdonor_scSplit <- DimPlot(seurat, reduction = "umap", group.by = "common_assignment_scSplit")
    ggsave(UMAPdonor_scSplit, filename = paste0(outdir,"UMAPdonor_scSplit.png"))
    UMAPdonor_freemuxlet <- DimPlot(seurat, reduction = "umap", group.by = "common_assignment_freemuxlet")
    ggsave(UMAPdonor_freemuxlet, filename = paste0(outdir,"UMAPdonor_freemuxlet.png"))

} else if (stage[1]=="QCfigurew_PoolRegressed"){
    seurat <- readRDS(paste0(dir,"output/Seurat/SCTnormalization_PoolRegressed/seurat_SCT_PoolRegressed.rds"))
    UMAPpool <- DimPlot(seurat, reduction = "umap", group.by = "Pool")
    ggsave(UMAPpool, filename = paste0(outdir,"UMAPpool.png"))
    UMAPmt <- FeaturePlot(seurat, reduction = "umap", feature = "percent.mt")
    ggsave(UMAPmt, filename = paste0(outdir,"UMAPmt.png"))
    UMAPrb <- FeaturePlot(seurat, reduction = "umap", feature = "percent.rb")
    ggsave(UMAPrb, filename = paste0(outdir,"UMAPrb.png"))
    UMAPlib <- FeaturePlot(seurat, reduction = "umap", feature = "nCount_RNA")
    ggsave(UMAPlib, filename = paste0(outdir,"UMAPlib.png"))
    UMAPgenes <- FeaturePlot(seurat, reduction = "umap", feature = "nFeature_RNA")
    ggsave(UMAPgenes, filename = paste0(outdir,"UMAPgenes.png"))
    UMAPlibSCT <- FeaturePlot(seurat, reduction = "umap", feature = "nCount_SCT")
    ggsave(UMAPlibSCT, filename = paste0(outdir,"UMAPlibSCT.png"))
    UMAPgenesSCT <- FeaturePlot(seurat, reduction = "umap", feature = "nFeature_SCT")
    ggsave(UMAPgenesSCT, filename = paste0(outdir,"UMAPgenesSCT.png"))
    UMAPdonor_souporcell <- DimPlot(seurat, reduction = "umap", group.by = "common_assignment_souporcell")
    ggsave(UMAPdonor_souporcell, filename = paste0(outdir,"UMAPdonor_souporcell.png"))
    UMAPdonor_scSplit <- DimPlot(seurat, reduction = "umap", group.by = "common_assignment_scSplit")
    ggsave(UMAPdonor_scSplit, filename = paste0(outdir,"UMAPdonor_scSplit.png"))
    UMAPdonor_freemuxlet <- DimPlot(seurat, reduction = "umap", group.by = "common_assignment_freemuxlet")
    ggsave(UMAPdonor_freemuxlet, filename = paste0(outdir,"UMAPdonor_freemuxlet.png"))


} else if (stage[1]=="QCfigurew_PoolDayRegressed"){
    seurat <- readRDS(paste0(dir,"output/Seurat/SCTnormalization_PoolDayRegressed/seurat_SCT_PoolDayRegressed.rds"))
    UMAPpool <- DimPlot(seurat, reduction = "umap", group.by = "Pool")
    ggsave(UMAPpool, filename = paste0(outdir,"UMAPpool.png"))
    UMAPmt <- FeaturePlot(seurat, reduction = "umap", feature = "percent.mt")
    ggsave(UMAPmt, filename = paste0(outdir,"UMAPmt.png"))
    UMAPrb <- FeaturePlot(seurat, reduction = "umap", feature = "percent.rb")
    ggsave(UMAPrb, filename = paste0(outdir,"UMAPrb.png"))
    UMAPlib <- FeaturePlot(seurat, reduction = "umap", feature = "nCount_RNA")
    ggsave(UMAPlib, filename = paste0(outdir,"UMAPlib.png"))
    UMAPgenes <- FeaturePlot(seurat, reduction = "umap", feature = "nFeature_RNA")
    ggsave(UMAPgenes, filename = paste0(outdir,"UMAPgenes.png"))
    UMAPlibSCT <- FeaturePlot(seurat, reduction = "umap", feature = "nCount_SCT")
    ggsave(UMAPlibSCT, filename = paste0(outdir,"UMAPlibSCT.png"))
    UMAPgenesSCT <- FeaturePlot(seurat, reduction = "umap", feature = "nFeature_SCT")
    ggsave(UMAPgenesSCT, filename = paste0(outdir,"UMAPgenesSCT.png"))
    UMAPdonor_souporcell <- DimPlot(seurat, reduction = "umap", group.by = "common_assignment_souporcell")
    ggsave(UMAPdonor_souporcell, filename = paste0(outdir,"UMAPdonor_souporcell.png"))
    UMAPdonor_scSplit <- DimPlot(seurat, reduction = "umap", group.by = "common_assignment_scSplit")
    ggsave(UMAPdonor_scSplit, filename = paste0(outdir,"UMAPdonor_scSplit.png"))
    UMAPdonor_freemuxlet <- DimPlot(seurat, reduction = "umap", group.by = "common_assignment_freemuxlet")
    ggsave(UMAPdonor_freemuxlet, filename = paste0(outdir,"UMAPdonor_freemuxlet.png"))


} else if (stage[1]=="QCfigures_Harmony_Pool"){
    seurat <- readRDS(paste0(dir,"output/Seurat/Harmony_Pool/NoOutliers_Pool_harmony.rds"))
    UMAPpool <- DimPlot(seurat, reduction = "umap", group.by = "Pool")
    ggsave(UMAPpool, filename = paste0(outdir,"UMAPpool.png"))
    UMAPmt <- FeaturePlot(seurat, reduction = "umap", feature = "percent.mt")
    ggsave(UMAPmt, filename = paste0(outdir,"UMAPmt.png"))
    UMAPrb <- FeaturePlot(seurat, reduction = "umap", feature = "percent.rb")
    ggsave(UMAPrb, filename = paste0(outdir,"UMAPrb.png"))
    UMAPlib <- FeaturePlot(seurat, reduction = "umap", feature = "nCount_RNA")
    ggsave(UMAPlib, filename = paste0(outdir,"UMAPlib.png"))
    UMAPgenes <- FeaturePlot(seurat, reduction = "umap", feature = "nFeature_RNA")
    ggsave(UMAPgenes, filename = paste0(outdir,"UMAPgenes.png"))
    UMAPdonor_souporcell <- DimPlot(seurat, reduction = "umap", group.by = "common_assignment_souporcell")
    ggsave(UMAPdonor_souporcell, filename = paste0(outdir,"UMAPdonor_souporcell.png"))
    UMAPdonor_scSplit <- DimPlot(seurat, reduction = "umap", group.by = "common_assignment_scSplit")
    ggsave(UMAPdonor_scSplit, filename = paste0(outdir,"UMAPdonor_scSplit.png"))
    UMAPdonor_freemuxlet <- DimPlot(seurat, reduction = "umap", group.by = "common_assignment_freemuxlet")
    ggsave(UMAPdonor_freemuxlet, filename = paste0(outdir,"UMAPdonor_freemuxlet.png"))


} else if (stage[1]=="QCfigures_Harmony_Day"){
    seurat <- readRDS(paste0(dir,"output/Seurat/Harmony_Day/NoOutliers_Day_harmony.rds"))
    UMAPpool <- DimPlot(seurat, reduction = "umap", group.by = "Pool")
    ggsave(UMAPpool, filename = paste0(outdir,"UMAPpool.png"))
    UMAPmt <- FeaturePlot(seurat, reduction = "umap", feature = "percent.mt")
    ggsave(UMAPmt, filename = paste0(outdir,"UMAPmt.png"))
    UMAPrb <- FeaturePlot(seurat, reduction = "umap", feature = "percent.rb")
    ggsave(UMAPrb, filename = paste0(outdir,"UMAPrb.png"))
    UMAPlib <- FeaturePlot(seurat, reduction = "umap", feature = "nCount_RNA")
    ggsave(UMAPlib, filename = paste0(outdir,"UMAPlib.png"))
    UMAPgenes <- FeaturePlot(seurat, reduction = "umap", feature = "nFeature_RNA")
    ggsave(UMAPgenes, filename = paste0(outdir,"UMAPgenes.png"))
    UMAPdonor_souporcell <- DimPlot(seurat, reduction = "umap", group.by = "common_assignment_souporcell")
    ggsave(UMAPdonor_souporcell, filename = paste0(outdir,"UMAPdonor_souporcell.png"))
    UMAPdonor_scSplit <- DimPlot(seurat, reduction = "umap", group.by = "common_assignment_scSplit")
    ggsave(UMAPdonor_scSplit, filename = paste0(outdir,"UMAPdonor_scSplit.png"))
    UMAPdonor_freemuxlet <- DimPlot(seurat, reduction = "umap", group.by = "common_assignment_freemuxlet")
    ggsave(UMAPdonor_freemuxlet, filename = paste0(outdir,"UMAPdonor_freemuxlet.png"))


} else if (stage[1]=="QCfigures_Harmony_DayPatient"){
    seurat <- readRDS(paste0(dir,"output/Seurat/Harmony_DayPatient/NoOutliers_DayPatient_harmony.rds"))
    UMAPpool <- DimPlot(seurat, reduction = "umap", group.by = "Pool")
    ggsave(UMAPpool, filename = paste0(outdir,"UMAPpool.png"))
    UMAPmt <- FeaturePlot(seurat, reduction = "umap", feature = "percent.mt")
    ggsave(UMAPmt, filename = paste0(outdir,"UMAPmt.png"))
    UMAPrb <- FeaturePlot(seurat, reduction = "umap", feature = "percent.rb")
    ggsave(UMAPrb, filename = paste0(outdir,"UMAPrb.png"))
    UMAPlib <- FeaturePlot(seurat, reduction = "umap", feature = "nCount_RNA")
    ggsave(UMAPlib, filename = paste0(outdir,"UMAPlib.png"))
    UMAPgenes <- FeaturePlot(seurat, reduction = "umap", feature = "nFeature_RNA")
    ggsave(UMAPgenes, filename = paste0(outdir,"UMAPgenes.png"))
    UMAPdonor_souporcell <- DimPlot(seurat, reduction = "umap", group.by = "FinalAssignment")
    ggsave(UMAPdonor_souporcell, filename = paste0(outdir,"UMAPdonor_souporcell.png"))



} else if (stage[1]=="QCfigures_Integration_DayPatient"){
    seurat <- readRDS(paste0(dir,"output/Seurat/Integration_DayPatient/NoOutliers_DayPatient_SeuratIntegratee.rds"))
    UMAPpool <- DimPlot(seurat, reduction = "umap", group.by = "Pool")
    ggsave(UMAPpool, filename = paste0(outdir,"UMAPpool.png"))
    UMAPmt <- FeaturePlot(seurat, reduction = "umap", feature = "percent.mt")
    ggsave(UMAPmt, filename = paste0(outdir,"UMAPmt.png"))
    UMAPrb <- FeaturePlot(seurat, reduction = "umap", feature = "percent.rb")
    ggsave(UMAPrb, filename = paste0(outdir,"UMAPrb.png"))
    UMAPlib <- FeaturePlot(seurat, reduction = "umap", feature = "nCount_RNA")
    ggsave(UMAPlib, filename = paste0(outdir,"UMAPlib.png"))
    UMAPgenes <- FeaturePlot(seurat, reduction = "umap", feature = "nFeature_RNA")
    ggsave(UMAPgenes, filename = paste0(outdir,"UMAPgenes.png"))
    UMAPdonor_souporcell <- DimPlot(seurat, reduction = "umap", group.by = "FinalAssignment")
    ggsave(UMAPdonor_souporcell, filename = paste0(outdir,"UMAPdonor_souporcell.png"))


} else if (stage[1]=="QCfigures_Integration_Day_then_Individual"){
    seurat <- readRDS(paste0(dir,"output/Seurat/Integration_Day_then_Individual/NoOutliers_Day_SeuratIntegratee.rds"))
    UMAPpool <- DimPlot(seurat, reduction = "umap", group.by = "Pool")
    ggsave(UMAPpool, filename = paste0(outdir,"UMAPpool.png"))
    UMAPmt <- FeaturePlot(seurat, reduction = "umap", feature = "percent.mt")
    ggsave(UMAPmt, filename = paste0(outdir,"UMAPmt.png"))
    UMAPrb <- FeaturePlot(seurat, reduction = "umap", feature = "percent.rb")
    ggsave(UMAPrb, filename = paste0(outdir,"UMAPrb.png"))
    UMAPlib <- FeaturePlot(seurat, reduction = "umap", feature = "nCount_RNA")
    ggsave(UMAPlib, filename = paste0(outdir,"UMAPlib.png"))
    UMAPgenes <- FeaturePlot(seurat, reduction = "umap", feature = "nFeature_RNA")
    ggsave(UMAPgenes, filename = paste0(outdir,"UMAPgenes.png"))
    UMAPdonor_souporcell <- DimPlot(seurat, reduction = "umap", group.by = "FinalAssignment")
    ggsave(UMAPdonor_souporcell, filename = paste0(outdir,"UMAPdonor_souporcell.png"))





} else if (stage[1]=="Compare_Harmony_Seurat_Clusters"){
    dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Seurat/"
    outdir <- paste0(dir,"Harmony_vs_Integration/")
    dir.create(outdir)
    seurat_harmony <- readRDS("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Seurat/Harmony_DayPatient/NoOutliers_DayPatient_harmony.rds")
    seurat_integration <- readRDS("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Seurat/Integration_Day_then_Individual/NoOutliers_Day_SeuratIntegratee.rds")
    

    harmony_resolution_file_list <- list.files(paste0(dir,"/Multi_Resolution_Harmony_DayPatient/"), pattern = "Identities_Resolution_")
    seurat_harmony_res_list <- lapply(harmony_resolution_file_list, FUN = function(x){
        read_delim(paste0(dir,"Multi_Resolution_Harmony_DayPatient/",x), delim = "\t")
    })
    
    harmony_ResolutionsDF <- do.call(cbind, seurat_harmony_res_list)
    colnames(harmony_ResolutionsDF) <- paste0("Harmony_", colnames(harmony_ResolutionsDF))
    rownames(harmony_ResolutionsDF) <- colnames(seurat_harmony)
    seurat_integration <- AddMetaData(seurat_integration, harmony_ResolutionsDF)

    for (res in colnames(seurat_integration@meta.data)[grep("Resolution", colnames(seurat_integration@meta.data))]){
        UMAP <- DimPlot(seurat_integration, reduction = "umap", group.by = res)
        ggsave(UMAP, filename = paste0(outdir,"seuratUMAP_harmony",res,".png"))
    }


    integration_resolution_file_list <- list.files(paste0(dir,"Multi_Resolution_Integration_Day_then_Individual/"), pattern = "Identities_Resolution_")
    seurat_integration_res_list <- lapply(integration_resolution_file_list, FUN = function(x){
        read_delim(paste0(dir,"Multi_Resolution_Integration_Day_then_Individual/",x), delim = "\t")
    })
    
    integration_ResolutionsDF <- do.call(cbind, seurat_integration_res_list)
    colnames(integration_ResolutionsDF) <- paste0("Seurat_",colnames(integration_ResolutionsDF))
    rownames(integration_ResolutionsDF) <- colnames(seurat_integration)
    seurat_harmony <- AddMetaData(seurat_harmony, integration_ResolutionsDF)

    for (res in colnames(seurat_harmony@meta.data)[grep("Resolution", colnames(seurat_harmony@meta.data))]){
        UMAP <- DimPlot(seurat_harmony, reduction = "umap", group.by = res)
        ggsave(UMAP, filename = paste0(outdir,"harmonyUMAP_seurat",res,".png"))
    }

    integration_ResolutionsDF$Barcode <- rownames(integration_ResolutionsDF)
    harmony_ResolutionsDF$Barcode <- rownames(harmony_ResolutionsDF)
    ResolutionsDF <- left_join(integration_ResolutionsDF,harmony_ResolutionsDF)

    plotlist <- list()

for (col in 1:ncol(ResolutionsDF)){
    ResolutionsDF[,col] <- as.factor(ResolutionsDF[,col])
}

for (res in seq(0,0.24,0.01)){
    seurat <- paste0("Seurat_Resolution_",res)
    harmony <- paste0("Harmony_Resolution_",res)
            plotlist[[paste0(seurat,harmony)]] <- ggplot(ResolutionsDF, aes_string(seurat)) +
                        geom_bar(position = "stack", aes_string(fill = harmony)) +
                        theme_classic()

}
    lapply(names(plotlist), function(x){ggsave(plotlist[[x]], filename = paste0(outdir,x,".png"))})
    
    lapply(plotlist1, function(x){ggsave(x, filename = paste0(outdir,"seuratClusters_fillHarmony_",res1,"_",res2,".png"))})
    lapply(plotlist2, function(x){ggsave(x, filename = paste0(outdir,"HarmonyClusters_fillSeurat_",res1,"_",res2,".png"))})




}


