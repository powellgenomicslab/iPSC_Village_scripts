##### Read in Arguments #####
print("Reading and assigning input arguments")

# dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
# outbase <- paste0(dir,"output/Seurat_w_Hash_freeze_thaw/")


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


library(scran)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(Seurat)
library(jcolors)
library(cowplot)
library(RColorBrewer)
library(readr)
library(purrr)
library(clustree)
library(reticulate)
library(BiocParallel)
library(coop)
library(harmony)
library(awtools)
library(Nebulosa)



print(sessionInfo())

### Set colors
cell_line_colors <- c(ppalette[c(2,4,5)], c(gpalette[c(2,3,4)]))
names(cell_line_colors) <- c("FSA0006","MBE1006","TOB0421","doublet","unassigned","combination_fail_remove")

PoolColors <- brewer.pal(n = 8, name = "Dark2")

### Set names ###
pools <- dir("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/")
old_name <- c("22_FSA", "29_MBE", "36_TOB00421_i_E8")
new_name <- c("FSA0006","MBE1006","TOB0421")
cell_key <- data.frame(old_name, new_name)


### Functions ###
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




if (stage == "Hashtags"){
    outdir <- paste0(outbase,stage,"/")
    dirs10x <- paste0("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/", pools, "/outs/filtered_feature_bc_matrix/")
    # pools_hash <- dir("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/", pattern = "DRENEA")
    dirs_hash <- paste0("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Barcodes_200213_A00152_0206_AHK3HYDRXX/Hashing/", pools, "/umi_count/")


    ## Read in expression data
    counts_list <- lapply(dirs10x, function(x){
        Read10X(x, gene.column = 1)
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
    hash_list6 <- lapply(hash_list[1:6], function(x){
        rownames(x) <- gsub("A0252-TGATGGCCTATTGGG","Brisbane2", rownames(x)) %>%
            gsub("A0253-TTCCGCCTCTCTTTG", "Brisbane3", .) %>%
            gsub("A0254-AGTAAGTTCAGCGTA", "Sydney1", .) %>%
            gsub("A0255-AAGTATCGTTTCGCA", "Sydney2", .) %>%
            gsub("A0256-GGTTGCCAGATGTCA", "Sydney3", .) %>%
            gsub("A0257-TGTCTTTCCTGCCAG", "Melbourne1", .) %>%
            gsub("A0258-CTCCTCTGCAATTAC", "Melbourne2", .)
        return(x)
    })

    hash_list <- c(hash_list6, hash_list[7:8])

    rownames(hash_list[["DRENEA_1"]]) <- gsub("A0251-GTCAACTCTTTAGCG","Brisbane1", rownames(hash_list[["DRENEA_1"]]))
    rownames(hash_list[["DRENEA_4"]]) <- gsub("A0251-GTCAACTCTTTAGCG","Brisbane1", rownames(hash_list[["DRENEA_4"]]))
    rownames(hash_list[["DRENEA_3"]]) <- gsub("A0251-GTCAACTCTTTAGCG","Melbourne3", rownames(hash_list[["DRENEA_3"]]))
    rownames(hash_list[["DRENEA_6"]]) <- gsub("A0251-GTCAACTCTTTAGCG","Melbourne3", rownames(hash_list[["DRENEA_6"]]))

    rownames(hash_list[["Village_A_Baseline"]]) <- gsub("A0251-GTCAACTCTTTAGCG","Sydney1", rownames(hash_list[["Village_A_Baseline"]]))
    rownames(hash_list[["Village_A_Baseline"]]) <- gsub("A0252-TGATGGCCTATTGGG","Sydney2", rownames(hash_list[["Village_A_Baseline"]]))
    rownames(hash_list[["Village_A_Baseline"]]) <- gsub("A0253-TTCCGCCTCTCTTTG","Sydney3", rownames(hash_list[["Village_A_Baseline"]]))
    rownames(hash_list[["Village_B_1_week"]]) <- gsub("A0254-AGTAAGTTCAGCGTA","Sydney1", rownames(hash_list[["Village_B_1_week"]]))
    rownames(hash_list[["Village_B_1_week"]]) <- gsub("A0255-AAGTATCGTTTCGCA","Sydney2", rownames(hash_list[["Village_B_1_week"]]))
    rownames(hash_list[["Village_B_1_week"]]) <- gsub("A0256-GGTTGCCAGATGTCA","Sydney3", rownames(hash_list[["Village_B_1_week"]]))

    ## Get just the cells that are also in the hash list ##
    counts_list <- lapply(names(counts_list), function(x){
        print(x)
        colnames(counts_list[[x]]) <- gsub("-1", "", colnames(counts_list[[x]]))
        counts_list[[x]] <- counts_list[[x]][,which(colnames(counts_list[[x]]) %in% colnames(hash_list[[x]]))]
        return(counts_list[[x]])
    })



    ##### Tried analysis with normal normalization as recommended on website #####
    ### Note: very low identification of doublets or negative cells ~3% per pool
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
        prop.table(table(x@meta.data$MULTI_classification))
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
    saveRDS(seurat_norm, paste0(outdir,"seurat_all_cells.rds"))
    # seurat_norm <- readRDS(paste0(outdir,"seurat_all_cells.rds"))

    seurat_norm@meta.data$Pool <- gsub("_[A-Z]{16}","",rownames(seurat_norm@meta.data)) %>% gsub("-1","",.)

    ### Fix barcodes so they are all the same (freeze-thaw have -1 at the end)
    barcode_key <- data.frame("Original_Barcode" = colnames(seurat_norm), "Updated_Barcode" = gsub("-1", "", colnames(seurat_norm)))
    
    seurat_norm <- RenameCells(seurat_norm, new.names = barcode_key$Updated_Barcode)

    seurat_norm@meta.data$Time <- ifelse((seurat_norm@meta.data$Pool == "DRENEA_1" | 
                                            seurat_norm@meta.data$Pool == "DRENEA_2" |
                                            seurat_norm@meta.data$Pool == "DRENEA_3"), "Baseline", 
                                    ifelse((seurat_norm@meta.data$Pool == "DRENEA_4" | 
                                            seurat_norm@meta.data$Pool == "DRENEA_5" |
                                            seurat_norm@meta.data$Pool == "DRENEA_6"),"Village Day 4",
                                    ifelse((seurat_norm@meta.data$Pool == "Village_A_Baseline"), "Thawed Village Day 0",
                                    ifelse((seurat_norm@meta.data$Pool == "Village_B_1_week"), "Thawed Village Day 7", NA))))


    ### Add gene names to feature metadata
    genes_list <- lapply(rev(pools), function(x){
        read.delim(paste0("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/",x,"/outs/filtered_feature_bc_matrix/features.tsv.gz"), header = FALSE, col.names = c("ENSG","Gene_ID", "Descriptor"))
    })

    genes <- do.call(rbind, genes_list)

    genes_unique <- rownames(seurat_norm)
    print(length(genes_unique))
    print(length(unique(genes_unique)))

    genes_unique_df <- genes[match(genes_unique, genes$ENSG),]
    rownames(genes_unique_df) <- genes_unique_df$ENSG
    genes_unique_df$Descriptor <- NULL

    write_delim(genes_unique_df, paste0(outdir,"gene_id_conversion.txt"), delim = "\t")

    seurat_norm[["RNA"]] <- AddMetaData(seurat_norm[["RNA"]], genes_unique_df)

    ## Add QC metrics
    RbGeneList <- read.delim(file = "/directflow/SCCGGroupShare/projects/DrewNeavin/References/RibosomalGeneList_GeneID_ENSG.txt")
    MtGeneList <- read.delim(file = "/directflow/SCCGGroupShare/projects/DrewNeavin/References/MtGeneList_GeneID_ENSG.txt")
    print("Calculating Mt %")
    seurat_norm <- PercentageFeatureSet(seurat_norm, features = MtGeneList$ENSG, col.name = "percent.mt")
    print("Calculating Rb %")
    RbGeneList <- RbGeneList[which(RbGeneList$ENSG %in% rownames(seurat_norm)),]
    seurat_norm <- PercentageFeatureSet(seurat_norm, features = RbGeneList$ENSG, col.name = "percent.rb")

    ### Save final seurat object
    saveRDS(seurat_norm, paste0(outdir,"AllCellsSeurat_norm_meta.rds"))


} else if (stage[1]=="Demultiplex_Doublet"){
    outdir <- paste0(outbase,stage,"/")
    seurat_norm <- readRDS(paste0(dir,"output/Seurat_w_Hash/Hashtags/AllCellsSeurat_norm_meta.rds"))

    ## Add demultiplexing and doublet detection results ##
    demultiplexing <- lapply(pools, function(x){
        read_delim(paste0(dir,"output/CombinedResults/", x, "/CombinedDropletAssignments.tsv"), delim = "\t")
    })
    names(demultiplexing) <- pools

    scrublet <- lapply(pools, function(x){
        read_delim(paste0(dir,"output/Scrublet/default0.85/",x,"/scrublet_results.txt"), delim ="\t")
    })
    names(scrublet) <- pools

    demultiplexing_doublets <- lapply(names(demultiplexing), function(x){
        left_join(demultiplexing[[x]], scrublet[[x]], by = "Barcode")
    })
    names(demultiplexing_doublets) <- pools

    demultiplexing_doublets <- lapply(names(demultiplexing_doublets), function(x){
        demultiplexing_doublets[[x]] <- as.data.frame(demultiplexing_doublets[[x]])
        rownames(demultiplexing_doublets[[x]]) <- paste0(x, "_",gsub("-1","",demultiplexing_doublets[[x]]$Barcode))
        return(demultiplexing_doublets[[x]])
    })
    lapply(demultiplexing_doublets, head)


    ### Update cell line names ###
    demultiplexing_doublets <- lapply(demultiplexing_doublets, function(x){
    for (row in 1:nrow(cell_key)){
        print(cell_key[row,"old_name"])
        print(cell_key[row,"new_name"])
        x[] <- lapply(x, function(y){
        gsub(cell_key[row,"old_name"], cell_key[row,"new_name"], y)
        })
    }
    return(x)
    })

    ### Replace NAs in demuxlet_DropletType and scSplit_Assignment and scSplit_DropletType with unassigned (very few but they exist)
    demultiplexing_doublets <- lapply(demultiplexing_doublets, function(x){
        temp <- x
        temp[c("demuxlet_DropletType", "scSplit_Assignment", "scSplit_DropletType")][is.na(temp[c("demuxlet_DropletType", "scSplit_Assignment", "scSplit_DropletType")])] <- "unassigned"
        temp$souporcell_Assignment <- ifelse(temp$souporcell_DropletType == "unassigned", "unassigned", ifelse(temp$souporcell_DropletType == "doublet", "doublet", temp$souporcell_Assignment))
        temp$demuxlet_Assignment <- ifelse(temp$demuxlet_DropletType == "doublet", "doublet", ifelse(temp$demuxlet_DropletType == "unassigned", "unassigned", temp$demuxlet_Assignment))
        return(temp)
    })


    ##### Rename each of the clusters and then combine into one column for identifying common assignment #####
    demultiplexing_doublets <- lapply(demultiplexing_doublets, function(x){
        colnames(x) <- c("Barcode",colnames(x)[2:ncol(x)])
        x$scSplit_Assignment <- paste0("scSplit_", x$scSplit_Assignment)
        x$souporcell_Assignment <- paste0("souporcell_",x$souporcell_Assignment)
        x$freemuxlet_Assignment <- paste0("freemuxlet_", x$freemuxlet_Assignment)
        x$vireo_Assignment <- paste0("vireo_", x$vireo_Assignment)
        x$demuxlet_Assignment <- paste0("demuxlet_",x$demuxlet_Assignment)
        x$combined_assignments <- paste(x$scSplit_Assignment,x$souporcell_Assignment,x$freemuxlet_Assignment,x$vireo_Assignment, x$demuxlet_Assignment, sep = "-")
        return(x)
    })
    names(demultiplexing_doublets) <- pools

    ### Make a table of hte most common combinations ###
    joined_assignment_counts <- lapply(demultiplexing_doublets, function(x){
        df <- as.data.frame(t(table(x$combined_assignments)))
        df <- df[order(df$Freq, decreasing = TRUE),]
        return(df)
    })
    lapply(joined_assignment_counts, head)


    ### Take just the top three to assign cells to those that were shared as singlets across all to reassign to a common assignment ###
    joined_assignment_counts_top <- lapply(joined_assignment_counts, function(x){
        x$Var1 <- NULL
        x <- x[1:3,]
        return(x)
    })

    joined_assignment_counts_top <- lapply(joined_assignment_counts_top, function(x){
        df <- separate(x,col = Var2, into = c("scSplit","souporcell","freemuxlet","vireo","demuxlet"), sep = "-")
        df$common_assignment <- gsub("demuxlet_","",df$demuxlet)
        return(df)
    })

    joined_assignment_key <- lapply(joined_assignment_counts_top, function(x){
        df <- pivot_longer(x,cols = c("scSplit","souporcell","freemuxlet","vireo","demuxlet"), names_to = "software")
        df$Freq <- NULL
        return(df)
    })
    names(joined_assignment_key) <- pools

    ### Pivot dataframe longer to bind with key ###
    demultiplexing_doublets_long <- lapply(demultiplexing_doublets, function(x){
        temp <- pivot_longer(x, cols = c("freemuxlet_Assignment", "scSplit_Assignment", "souporcell_Assignment", "vireo_Assignment", "demuxlet_Assignment"), names_to = "Software", values_to = "Original_Assignment")
        temp$Software <- gsub("_Assignment","", temp$Software)
        return(temp)
    })
    names(demultiplexing_doublets_long) <- pools

    demultiplexing_doublets_long <- lapply(names(demultiplexing_doublets_long), function(x){
        df <- left_join(demultiplexing_doublets_long[[x]], joined_assignment_key[[x]], by = c("Original_Assignment" = "value", "Software" = "software"))
        df$common_assignment <- ifelse(gsub("[a-zA-Z]+_","",df$Original_Assignment) == "unassigned", "unassigned", ifelse(gsub("[a-zA-Z]+_","",df$Original_Assignment) == "doublet", "doublet", df$common_assignment))
        df$Pool <- x
        return(df)
    })
    names(demultiplexing_doublets_long) <- pools


    ##### Combine the pool dataframes together to do a facet_grid plot
    joined_df4facet <- do.call(rbind,demultiplexing_doublets_long)
    joined_df4facet

    plot <- ggplot(joined_df4facet, aes(common_assignment, fill = Software)) +
        geom_bar(position="dodge") +
        facet_wrap(~Pool, nrow = 2) +
        theme_classic() +
        scale_fill_jcolors(palette = "pal4") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(plot, filename = paste0(outdir,"6Pools_software_bar_plot.png"))
         

    ##### demultiplexing_doublets_long for droplet type #####
    demultiplexing_doublets_droptype_long <- lapply(names(demultiplexing_doublets), function(x){
        temp <- pivot_longer(demultiplexing_doublets[[x]], cols = colnames(demultiplexing_doublets[[x]])[grep("DropletType", colnames(demultiplexing_doublets[[x]]))], names_to = "Software", values_to = "DropletType")
        temp$Software <- gsub("_DropletType","", temp$Software)
        temp$DropletType <- ifelse(temp$DropletType == "unassigned", "unassigned", ifelse(temp$DropletType == "doublet", "doublet", ifelse(temp$DropletType == "FSA0006", "singlet", ifelse(temp$DropletType == "MBE1006", "singlet", ifelse(temp$DropletType == "TOB0421", "singlet", temp$DropletType)))))
        temp$Pool <- x
        return(temp)
    })
    names(demultiplexing_doublets_droptype_long) <- pools


    ##### Combine the pool dataframes together to do a facet_grid plot
    joined_doublet_df4facet <- do.call(rbind, demultiplexing_doublets_droptype_long)
    joined_doublet_df4facet

    colfunc <- jcolors_contin(pal = "pal4")
    jcols   <- colfunc(7)
    n <- length(jcols)

    plot <- ggplot(joined_doublet_df4facet, aes(DropletType, fill = Software)) +
        geom_bar(position="dodge") +
        facet_wrap(~Pool, nrow = 2) +
        theme_classic() +
        scale_fill_manual(values = jcols) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(plot, filename = paste0(outdir,"6Pools_software_bar_plot_droptype.png"))


    ##### Will use majority to call cell type and individual assignment
    ### 1. Pivot wider the dataframes
    droplettype_wide <- pivot_wider(joined_doublet_df4facet[,c("Barcode", "Software", "DropletType", "Pool")], names_from = "Software", values_from = "DropletType", names_prefix = "DropletType_")
    assignment_wide <- pivot_wider(joined_df4facet[,c("Barcode", "Software", "common_assignment", "Pool")], names_from = "Software", values_from = "common_assignment", names_prefix = "Assignment_")

    ### 2. Add up the number of singlets, doublets, FSA0006, MBE10006, TOB0421 (place sums in their own columns)
    for (assignment in c(unique(joined_doublet_df4facet$DropletType))){
        droplettype_wide[,c(assignment)] <- rowSums(droplettype_wide == assignment)
    }

    for (assignment in c(unique(joined_df4facet$common_assignment))){
        assignment_wide[,c(assignment)] <- rowSums(assignment_wide == assignment)
    }

    assignment_wide$doublet <- NULL
    assignment_wide$unassigned <- NULL

    ### 3. Combine DropletType and Assignment Results into a single dataframe with barcodes using Barcode, Software and Pool as the join_by argument
    droplettype_assignments <- left_join(droplettype_wide, assignment_wide, by = c("Barcode", "Pool"))



    ### 4. Use if else statements to call singlets and individuals assigned to those cells
    droplettype_assignments$Final_Assignment <- ifelse(droplettype_assignments$doublet >= 4, "doublet", 
                                                    ifelse(droplettype_assignments$unassigned >= 4, "unassigned", 
                                                    ifelse(droplettype_assignments$FSA0006 >= 3, "FSA0006",
                                                    ifelse(droplettype_assignments$TOB0421 >= 3, "TOB0421",
                                                    ifelse(droplettype_assignments$MBE1006 >= 3, "MBE1006", "unassigned")))))
    table(droplettype_assignments$Final_Assignment)

    ##    doublet    FSA0006    MBE1006    TOB0421 unassigned 
    ##    7877      48109      31576      26692      30746 
    ### Low doublet rate but high unassigned (many because were called as singlets but couldn't be easily assigned to a specific individual) but combined rate is ~26.6% removed which is about right


    ### Plot final assignments ###
    droplettype_assignments$Final_Assignment <- factor(droplettype_assignments$Final_Assignment, levels = c("FSA0006", "MBE1006", "TOB0421", "doublet", "unassigned"))

    plot <- ggplot(droplettype_assignments, aes(Final_Assignment, fill = Final_Assignment)) +
        geom_bar(position="dodge") +
        facet_wrap(~Pool, nrow = 2) +
        theme_classic() +
        scale_fill_manual(values = cell_line_colors) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(plot, filename = paste0(outdir,"6Pools_software_bar_plot_final_assignments.png"))




    droplettype_assignments <- as.data.frame(droplettype_assignments)
    rownames(droplettype_assignments) <- gsub("-1","",paste0(droplettype_assignments$Pool, "_", droplettype_assignments$Barcode))


    ### Add final assignments to seurat object ###
    seurat_norm <- AddMetaData(seurat_norm, droplettype_assignments)


    ### Make pre-QC figures ###
    seurat_norm <- NormalizeData(seurat_norm, verbose = TRUE)
    seurat_norm <- FindVariableFeatures(seurat_norm, selection.method = "mean.var.plot")
    seurat_norm <- ScaleData(seurat_norm, features = VariableFeatures(seurat_norm))
    seurat_norm <- RunPCA(seurat_norm, features = VariableFeatures(object = seurat_norm))
    seurat_norm <- FindNeighbors(seurat_norm, dims = 1:10)
    seurat_norm <- FindClusters(seurat_norm, resolution = 0.5)
    seurat_norm <- RunUMAP(seurat_norm, dims = 1:10)

    ### QC Figures pre filtering ###
    plot_mt_pct <- VlnPlot(seurat_norm, features = c( "percent.mt"), group.by = "Pool", pt.size = 0, cols = PoolColors)
    ggsave(plot_mt_pct, filename = paste0(outdir,"Mt_pct_vln.png"))

    plot_rb_pct <- VlnPlot(seurat_norm, features = c( "percent.rb"), group.by = "Pool", pt.size = 0, cols = PoolColors)
    ggsave(plot_rb_pct, filename = paste0(outdir,"Rb_pct_vln.png"))

    plot_n_count <- VlnPlot(seurat_norm, features = c( "nCount_RNA"), group.by = "Pool", pt.size = 0, cols = PoolColors)
    ggsave(plot_n_count, filename = paste0(outdir,"N_count_vln.png"))

    plot_nFeature_RNA <- VlnPlot(seurat_norm, features = c( "nFeature_RNA"), group.by = "Pool", pt.size = 0, cols = PoolColors)
    ggsave(plot_nFeature_RNA, filename = paste0(outdir,"nFeature_RNA_vln.png"))

    lib_mt <- FeatureScatter(seurat_norm, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Pool", cols = alpha(PoolColors, 0.25))
    ggsave(lib_mt, filename = paste0(outdir,"lib_mt.png"))

    lib_genes <- FeatureScatter(seurat_norm, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Pool", cols = alpha(PoolColors, 0.25))
    ggsave(lib_genes, filename = paste0(outdir,"lib_genes.png"))

    nebulosa_mt_umap <- plot_density(seurat_norm, "percent.mt", pal = "plasma")
    ggsave(nebulosa_mt_umap, filename = paste0(outdir,"mt_percent_umap.png"))

    nebulosa_rb_umap <- plot_density(seurat_norm, "percent.rb", pal = "plasma")
    ggsave(nebulosa_rb_umap, filename = paste0(outdir,"rb_percent_umap.png"))

    umap_Time <- DimPlot(seurat_norm, group.by = "Time")
    ggsave(umap_Time, filename = paste0(outdir,"Time_umap.png"))

    umap_Pool <- DimPlot(seurat_norm, group.by = "Pool")
    ggsave(umap_Pool, filename = paste0(outdir,"pool_umap.png"))

    umap_location <- DimPlot(seurat_norm, group.by = "MULTI_ID")
    ggsave(umap_location, filename = paste0(outdir,"location_umap.png"))

    mt_umap <- FeaturePlot(seurat_norm, features = "percent.mt")
    ggsave(mt_umap, filename = paste0(outdir,"mt_percent_umap_seurat.png"))


    seurat_norm@meta.data$Site_rep <- seurat_norm@meta.data$MULTI_ID
    any(is.na(seurat_norm@meta.data$Site_rep))

    p_Site_by_CellNumber <- ggplot(seurat_norm@meta.data, aes(factor(paste(seurat_norm@meta.data$Site_rep, seurat_norm@meta.data$Time), 
            levels = c("Brisbane1 Baseline","Brisbane2 Baseline","Brisbane3 Baseline", "Brisbane1 Village Day 4","Brisbane2 Village Day 4","Brisbane3 Village Day 4",
                        "Sydney1 Baseline","Sydney2 Baseline","Sydney3 Baseline", "Sydney1 Village Day 4","Sydney2 Village Day 4","Sydney3 Village Day 4", 
                        "Melbourne1 Baseline","Melbourne2 Baseline","Melbourne3 Baseline",
                        "Melbourne1 Village Day 4","Melbourne2 Village Day 4","Melbourne3 Village Day 4", 
                        "Sydney1 Thawed Village Day 0", "Sydney2 Thawed Village Day 0", "Sydney3 Thawed Village Day 0", 
                        "Sydney1 Thawed Village Day 7", "Sydney2 Thawed Village Day 7", "Sydney3 Thawed Village Day 7",
                        "Doublet Baseline","Negative Baseline", "Doublet Village Day 4","Negative Village Day 4", "Doublet Thawed Village Day 0", "Negative Thawed Village Day 0", "Doublet Thawed Village Day 7", "Negative Thawed Village Day 7")), fill = seurat_norm@meta.data$Final_Assignment)) +
        geom_bar(stat = "count", position = "dodge") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(y="Total Cells", x = "Site Sample", fill = "Cell Line") +
        scale_fill_manual(values = cell_line_colors) 

    ggsave(p_Site_by_CellNumber, filename = paste0(outdir,"Site_by_CellNumber.png"))

    p_Site_by_CellNumber <- ggplot(seurat_norm@meta.data, aes(factor(paste(seurat_norm@meta.data$Site_rep, seurat_norm@meta.data$Time), 
            levels = c("Brisbane1 Baseline","Brisbane2 Baseline","Brisbane3 Baseline", "Brisbane1 Village Day 4","Brisbane2 Village Day 4","Brisbane3 Village Day 4",
                        "Sydney1 Baseline","Sydney2 Baseline","Sydney3 Baseline", "Sydney1 Village Day 4","Sydney2 Village Day 4","Sydney3 Village Day 4", 
                        "Melbourne1 Baseline","Melbourne2 Baseline","Melbourne3 Baseline",
                        "Melbourne1 Village Day 4","Melbourne2 Village Day 4","Melbourne3 Village Day 4", 
                        "Sydney1 Thawed Village Day 0", "Sydney2 Thawed Village Day 0", "Sydney3 Thawed Village Day 0", 
                        "Sydney1 Thawed Village Day 7", "Sydney2 Thawed Village Day 7", "Sydney3 Thawed Village Day 7",
                        "Doublet Baseline","Negative Baseline", "Doublet Village Day 4","Negative Village Day 4", "Doublet Thawed Village Day 0", "Negative Thawed Village Day 0", "Doublet Thawed Village Day 7", "Negative Thawed Village Day 7")), fill = factor(seurat_norm@meta.data$Final_Assignment, levels = c("FSA0006","MBE1006","TOB0421","doublet","unassigned","combination_fail_remove")))) +
        geom_bar(stat = "count", position = "stack") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(y="Total Cells", x = "Site Sample", fill = "Cell Line") +
        scale_fill_manual(values = cell_line_colors) 

    ggsave(p_Site_by_CellNumber, filename = paste0(outdir,"Site_by_CellNumberStack.png"))

    p_Site_by_CellProportion <- ggplot(seurat_norm@meta.data, aes(factor(paste(seurat_norm@meta.data$Site_rep, seurat_norm@meta.data$Time), 
            levels = c("Brisbane1 Baseline","Brisbane2 Baseline","Brisbane3 Baseline", "Brisbane1 Village Day 4","Brisbane2 Village Day 4","Brisbane3 Village Day 4",
                        "Sydney1 Baseline","Sydney2 Baseline","Sydney3 Baseline", "Sydney1 Village Day 4","Sydney2 Village Day 4","Sydney3 Village Day 4", 
                        "Melbourne1 Baseline","Melbourne2 Baseline","Melbourne3 Baseline",
                        "Melbourne1 Village Day 4","Melbourne2 Village Day 4","Melbourne3 Village Day 4", 
                        "Sydney1 Thawed Village Day 0", "Sydney2 Thawed Village Day 0", "Sydney3 Thawed Village Day 0", 
                        "Sydney1 Thawed Village Day 7", "Sydney2 Thawed Village Day 7", "Sydney3 Thawed Village Day 7",
                        "Doublet Baseline","Negative Baseline", "Doublet Village Day 4","Negative Village Day 4", "Doublet Thawed Village Day 0", "Negative Thawed Village Day 0", "Doublet Thawed Village Day 7", "Negative Thawed Village Day 7")), fill = factor(seurat_norm@meta.data$Final_Assignment, levels = c("FSA0006","MBE1006","TOB0421","doublet","unassigned","combination_fail_remove")))) +
        geom_bar(position = "fill") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(y="Total Cells", x = "Site Sample", fill = "Cell Line") +
        scale_fill_manual(values = cell_line_colors) 

    ggsave(p_Site_by_CellProportion, filename = paste0(outdir,"Site_by_CellProportion.png"))


    print("The seurat object before removing doublets:")
    print(seurat_norm)

    ### Remove doublets selected ###
    Idents(seurat_norm) <- "Final_Assignment"
    seurat_norm <- subset(seurat_norm,  idents = c("FSA0006", "TOB0421", "MBE1006"))
    print(seurat_norm)
    Idents(seurat_norm) <- "Site_rep"
    seurat_norm <- subset(seurat_norm,  idents = c("Sydney1","Brisbane1", "Melbourne1", "Brisbane2", "Sydney2", "Melbourne2", "Melbourne3", "Sydney3", "Brisbane3"))
    print("The seurat object after removing doublets")
    print(seurat_norm)

    saveRDS(seurat_norm, paste0(outdir,"Seurat_noDoublets_norm_meta.rds"))
    seurat_norm <- readRDS(paste0(outdir,"Seurat_noDoublets_norm_meta.rds"))


    p_Site_by_CellNumber <- ggplot(seurat_norm@meta.data, aes(factor(paste(seurat_norm@meta.data$Site_rep, seurat_norm@meta.data$Time), 
            levels = c("Brisbane1 Baseline","Brisbane2 Baseline","Brisbane3 Baseline", "Brisbane1 Village Day 4","Brisbane2 Village Day 4","Brisbane3 Village Day 4",
                        "Sydney1 Baseline","Sydney2 Baseline","Sydney3 Baseline", "Sydney1 Village Day 4","Sydney2 Village Day 4","Sydney3 Village Day 4", 
                        "Melbourne1 Baseline","Melbourne2 Baseline","Melbourne3 Baseline",
                        "Melbourne1 Village Day 4","Melbourne2 Village Day 4","Melbourne3 Village Day 4", 
                        "Sydney1 Thawed Village Day 0", "Sydney2 Thawed Village Day 0", "Sydney3 Thawed Village Day 0", 
                        "Sydney1 Thawed Village Day 7", "Sydney2 Thawed Village Day 7", "Sydney3 Thawed Village Day 7",
                        "Doublet Baseline","Negative Baseline", "Doublet Village Day 4","Negative Village Day 4", "Doublet Thawed Village Day 0", "Negative Thawed Village Day 0", "Doublet Thawed Village Day 7", "Negative Thawed Village Day 7")), fill = factor(seurat_norm@meta.data$Final_Assignment, levels = c("FSA0006","MBE1006","TOB0421")))) +
        geom_bar(stat = "count", position = "dodge") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(y="Total Cells", x = "Site Sample", fill = "Cell Line") +
        scale_fill_manual(values = cell_line_colors) 

    ggsave(p_Site_by_CellNumber, filename = paste0(outdir,"Site_by_CellNumber_singlets.png"))

    p_Site_by_CellNumber <- ggplot(seurat_norm@meta.data, aes(factor(paste(seurat_norm@meta.data$Site_rep, seurat_norm@meta.data$Time), 
            levels = c("Brisbane1 Baseline","Brisbane2 Baseline","Brisbane3 Baseline", "Brisbane1 Village Day 4","Brisbane2 Village Day 4","Brisbane3 Village Day 4",
                        "Sydney1 Baseline","Sydney2 Baseline","Sydney3 Baseline", "Sydney1 Village Day 4","Sydney2 Village Day 4","Sydney3 Village Day 4", 
                        "Melbourne1 Baseline","Melbourne2 Baseline","Melbourne3 Baseline",
                        "Melbourne1 Village Day 4","Melbourne2 Village Day 4","Melbourne3 Village Day 4", 
                        "Sydney1 Thawed Village Day 0", "Sydney2 Thawed Village Day 0", "Sydney3 Thawed Village Day 0", 
                        "Sydney1 Thawed Village Day 7", "Sydney2 Thawed Village Day 7", "Sydney3 Thawed Village Day 7",
                        "Doublet Baseline","Negative Baseline", "Doublet Village Day 4","Negative Village Day 4", "Doublet Thawed Village Day 0", "Negative Thawed Village Day 0", "Doublet Thawed Village Day 7", "Negative Thawed Village Day 7")), fill = factor(seurat_norm@meta.data$Final_Assignment, levels = c("FSA0006","MBE1006","TOB0421","doublet","unassigned","combination_fail_remove")))) +
        geom_bar(stat = "count", position = "stack") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(y="Total Cells", x = "Site Sample", fill = "Cell Line") +
        scale_fill_manual(values = cell_line_colors) 

    ggsave(p_Site_by_CellNumber, filename = paste0(outdir,"Site_by_CellNumberStack_singlets.png"))

    p_Site_by_CellProportion <- ggplot(seurat_norm@meta.data, aes(factor(paste(seurat_norm@meta.data$Site_rep, seurat_norm@meta.data$Time), 
            levels = c("Brisbane1 Baseline","Brisbane2 Baseline","Brisbane3 Baseline", "Brisbane1 Village Day 4","Brisbane2 Village Day 4","Brisbane3 Village Day 4",
                        "Sydney1 Baseline","Sydney2 Baseline","Sydney3 Baseline", "Sydney1 Village Day 4","Sydney2 Village Day 4","Sydney3 Village Day 4", 
                        "Melbourne1 Baseline","Melbourne2 Baseline","Melbourne3 Baseline",
                        "Melbourne1 Village Day 4","Melbourne2 Village Day 4","Melbourne3 Village Day 4", 
                        "Sydney1 Thawed Village Day 0", "Sydney2 Thawed Village Day 0", "Sydney3 Thawed Village Day 0", 
                        "Sydney1 Thawed Village Day 7", "Sydney2 Thawed Village Day 7", "Sydney3 Thawed Village Day 7",
                        "Doublet Baseline","Negative Baseline", "Doublet Village Day 4","Negative Village Day 4", "Doublet Thawed Village Day 0", "Negative Thawed Village Day 0", "Doublet Thawed Village Day 7", "Negative Thawed Village Day 7")), fill = factor(seurat_norm@meta.data$Final_Assignment, levels = c("FSA0006","MBE1006","TOB0421")))) +
        geom_bar(position = "fill") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(y="Total Cells", x = "Site Sample", fill = "Cell Line") +
        scale_fill_manual(values = cell_line_colors) 

    ggsave(p_Site_by_CellProportion, filename = paste0(outdir,"Site_by_CellProportion_singlets.png"))
 
} else if (stage[1]=="QC"){
    seurat <- readRDS(paste0(outbase,"/Demultiplex_Doublet/Seurat_noDoublets_norm_meta.rds"))
    outdir <- paste0(outbase,"/",stage,"/")
    pools <- dir("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/")

    ## Calculate MADs
    seurat <- mad_function(seurat, "percent.mt",4)
    seurat <- mad_function(seurat, "percent.rb",5)
    seurat <- mad_function(seurat, "nCount_RNA", 5)
    seurat <- mad_function(seurat, "nFeature_RNA", 3.4)


    plot_mt_pct <- VlnPlot(seurat, features = c( "percent.mt"), group.by = "Pool", pt.size = 0, cols = PoolColors) +
                    geom_hline(yintercept=0, linetype="dashed", color = "black") +
                    geom_hline(yintercept=14.10345, linetype="dashed", color = "black")
    ggsave(plot_mt_pct, filename = paste0(outdir,"Mt_pct_vln.png"))

    plot_rb_pct <- VlnPlot(seurat, features = c( "percent.rb"), group.by = "Pool", pt.size = 0, cols = PoolColors) +
                    geom_hline(yintercept=14.20789, linetype="dashed", color = "black") +
                    geom_hline(yintercept=37.16034, linetype="dashed", color = "black")
    ggsave(plot_rb_pct, filename = paste0(outdir,"Rb_pct_vln.png"))

    plot_n_count <- VlnPlot(seurat, features = c( "nCount_RNA"), group.by = "Pool", pt.size = 0, cols = PoolColors) +
                    geom_hline(yintercept=0, linetype="dashed", color = "black") +
                    geom_hline(yintercept=59769.84, linetype="dashed", color = "black")
    ggsave(plot_n_count, filename = paste0(outdir,"N_count_vln.png"))

    plot_nFeature_RNA <- VlnPlot(seurat, features = c( "nFeature_RNA"), group.by = "Pool", pt.size = 0, cols = PoolColors) +
                    geom_hline(yintercept=670.7134, linetype="dashed", color = "black") +
                    geom_hline(yintercept=8205.287, linetype="dashed", color = "black")
    ggsave(plot_nFeature_RNA, filename = paste0(outdir,"nFeature_RNA_vln.png"))

    lib_mt <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Pool", cols = alpha(brewer.pal(n = 8, name = "Dark2"), 0.25)) +
                    geom_vline(xintercept=0, linetype="dashed", color = "black") +
                    geom_vline(xintercept=59769.84, linetype="dashed", color = "black") +
                    geom_hline(yintercept=0, linetype="dashed", color = "black") +
                    geom_hline(yintercept=14.10345, linetype="dashed", color = "black") 
    ggsave(lib_mt, filename = paste0(outdir,"NoDoublets_lib_mt.png"))

    lib_genes <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Pool", cols = alpha(brewer.pal(n = 8, name = "Dark2"), 0.25)) +
                    geom_vline(xintercept=0, linetype="dashed", color = "black") +
                    geom_vline(xintercept=59769.84, linetype="dashed", color = "black") +
                    geom_hline(yintercept=670.7134, linetype="dashed", color = "black") +
                    geom_hline(yintercept=8205.287, linetype="dashed", color = "black")
    ggsave(lib_genes, filename = paste0(outdir,"NoDoublets_lib_genes.png"))


    ### Vioplin plots by cell type and site ###
    p_violin_mt_pct_line_site <- ggplot(seurat@meta.data, aes(factor(Time, levels = c("Day_0", "Day_4")), percent.mt, fill = factor(FinalAssignment, levels = c("FSA006", "MBE1006","TOB421")))) +
        geom_violin() +
        theme_classic() +
        labs(y="Mitochondrial Percent", x = "Sample Day", fill = "Cell Line") +
        scale_fill_manual(values = cell_line_colors) +
        facet_wrap(vars( MULTI_ID))
    ggsave(p_violin_mt_pct_line_site, filename = paste0(outdir,"violin_mt_pct_line_site.png"))

    p_violin_rb_pct_line_site <- ggplot(seurat@meta.data, aes(factor(Time, levels = c("Day_0", "Day_4")), percent.rb, fill = factor(FinalAssignment, levels = c("FSA006", "MBE1006","TOB421")))) +
        geom_violin() +
        theme_classic() +
        labs(y="Ribosomal Percent", x = "Sample Day", fill = "Cell Line") +
        scale_fill_manual(values = cell_line_colors) +
        facet_wrap(vars( MULTI_ID))
    ggsave(p_violin_rb_pct_line_site, filename = paste0(outdir,"violin_rb_pct_line_site.png"))

    p_violin_n_count_line_site <- ggplot(seurat@meta.data, aes(factor(Time, levels = c("Day_0", "Day_4")), nCount_RNA, fill = factor(FinalAssignment, levels = c("FSA006", "MBE1006","TOB421")))) +
        geom_violin() +
        theme_classic() +
        labs(y="Number UMIs", x = "Sample Day", fill = "Cell Line") +
        scale_fill_manual(values = cell_line_colors) +
        facet_wrap(vars( MULTI_ID))
    ggsave(p_violin_n_count_line_site, filename = paste0(outdir,"violin_n_count_line_site.png"))

    p_violin_n_feature_line_site <- ggplot(seurat@meta.data, aes(factor(Time, levels = c("Day_0", "Day_4")), nFeature_RNA, fill = factor(FinalAssignment, levels = c("FSA006", "MBE1006","TOB421")))) +
        geom_violin() +
        theme_classic() +
        labs(y="Number Features", x = "Sample Day", fill = "Cell Line") +
        scale_fill_manual(values = cell_line_colors) +
        facet_wrap(vars( MULTI_ID))
    ggsave(p_violin_n_feature_line_site, filename = paste0(outdir,"violin_n_feature_line_site.png"))



    ##### Remove the outliers #####
    ### Only filter based on mt% since others are pretty well within expected distributions
    print(seurat)
    seurat_filt <- subset(seurat, subset = percent.mt_mad == "NotOutlier") 

    print(seurat_filt)
    saveRDS(seurat_filt, paste0(outdir,"seurat_norm_NoOutliers.rds"))


    plot_mt_pct <- VlnPlot(seurat_filt, features = c( "percent.mt"), group.by = "Pool", pt.size = 0, cols = PoolColors)
    ggsave(plot_mt_pct, filename = paste0(outdir,"Mt_pct_vln_NoOutliers.png"))

    plot_rb_pct <- VlnPlot(seurat_filt, features = c( "percent.rb"), group.by = "Pool", pt.size = 0, cols = PoolColors)
    ggsave(plot_rb_pct, filename = paste0(outdir,"Rb_pct_vln_NoOutliers.png"))

    plot_n_count <- VlnPlot(seurat_filt, features = c( "nCount_RNA"), group.by = "Pool", pt.size = 0, cols = PoolColors)
    ggsave(plot_n_count, filename = paste0(outdir,"N_count_vln_NoOutliers.png"))

    plot_nFeature_RNA <- VlnPlot(seurat_filt, features = c( "nFeature_RNA"), group.by = "Pool", pt.size = 0, cols = PoolColors)
    ggsave(plot_nFeature_RNA, filename = paste0(outdir,"nFeature_RNA_vln_NoOutliers.png"))

    lib_mt <- FeatureScatter(seurat_filt, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Pool", cols = alpha(brewer.pal(n = 8, name = "Dark2"), 0.25))
    ggsave(lib_mt, filename = paste0(outdir,"NoDoublets_lib_mt_NoOutliers.png"))
    
    lib_genes <- FeatureScatter(seurat_filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Pool", cols = alpha(brewer.pal(n = 8, name = "Dark2"), 0.25))
    ggsave(lib_genes, filename = paste0(outdir,"NoDoublets_lib_genes_NoOutliers.png"))


    ### Vioplin plots by cell type and site ###
    p_violin_mt_pct_line_site <- ggplot(seurat_filt@meta.data, aes(factor(Time, levels = c("Baseline", "Village Day 4", "Thawed Village Day 0", "Thawed Village Day 7")), percent.mt, fill = factor(Final_Assignment, levels = c("FSA0006", "MBE1006","TOB0421")))) +
        geom_violin() +
        theme_classic() +
        labs(y="Mitochondrial Percent", x = "Sample Day", fill = "Cell Line") +
        scale_fill_manual(values = cell_line_colors) +
        facet_wrap(vars( MULTI_ID))
    ggsave(p_violin_mt_pct_line_site, filename = paste0(outdir,"violin_mt_pct_line_site_NoOutliers.png"))

    p_violin_rb_pct_line_site <- ggplot(seurat_filt@meta.data, aes(factor(Time, levels = c("Baseline", "Village Day 4", "Thawed Village Day 0", "Thawed Village Day 7")), percent.rb, fill = factor(Final_Assignment, levels = c("FSA0006", "MBE1006","TOB0421")))) +
        geom_violin() +
        theme_classic() +
        labs(y="Ribosomal Percent", x = "Sample Day", fill = "Cell Line") +
        scale_fill_manual(values = cell_line_colors) +
        facet_wrap(vars( MULTI_ID))

    ggsave(p_violin_rb_pct_line_site, filename = paste0(outdir,"violin_rb_pct_line_site_NoOutliers.png"))

    p_violin_n_count_line_site <- ggplot(seurat_filt@meta.data, aes(factor(Time, levels = c("Baseline", "Village Day 4", "Thawed Village Day 0", "Thawed Village Day 7")), nCount_RNA, fill = factor(Final_Assignment, levels = c("FSA0006", "MBE1006","TOB0421")))) +
        geom_violin() +
        theme_classic() +
        labs(y="Number UMIs", x = "Sample Day", fill = "Cell Line") +
        scale_fill_manual(values = cell_line_colors) +
        facet_wrap(vars( MULTI_ID))
    ggsave(p_violin_n_count_line_site, filename = paste0(outdir,"violin_n_count_line_site_NoOutliers.png"))

    p_violin_n_feature_line_site <- ggplot(seurat_filt@meta.data, aes(factor(Time, levels = c("Baseline", "Village Day 4", "Thawed Village Day 0", "Thawed Village Day 7")), nFeature_RNA, fill = factor(Final_Assignment, levels = c("FSA0006", "MBE1006","TOB0421")))) +
        geom_violin() +
        theme_classic() +
        labs(y="Number Features", x = "Sample Day", fill = "Cell Line") +
        scale_fill_manual(values = cell_line_colors) +
        facet_wrap(vars( MULTI_ID))
    ggsave(p_violin_n_feature_line_site, filename = paste0(outdir,"violin_n_feature_line_site_NoOutliers.png"))



} else if (stage == "CellCycle"){

    print("Starting cell cycle step")
    seurat <- readRDS(paste0(dir,"/output/Seurat_w_Hash_freeze_thaw/QC/seurat_norm_NoOutliers.rds"))

    seurat <- NormalizeData(seurat, verbose = TRUE)
    seurat <- FindVariableFeatures(seurat, selection.method = "mean.var.plot")
    seurat <- ScaleData(seurat, features = VariableFeatures(seurat))
    seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))
    seurat <- FindNeighbors(seurat, dims = 1:10)
    seurat <- FindClusters(seurat, resolution = 0.5)
    seurat <- RunUMAP(seurat, dims = 1:10)

    s.genes <- cc.genes$s.genes
    s.genes_ENSG <- seurat[["RNA"]][[]][which(seurat[["RNA"]][[]]$Gene_ID %in% s.genes),"ENSG"]
    g2m.genes <- cc.genes$g2m.genes
    g2m.genes_ENSG <- seurat[["RNA"]][[]][which(seurat[["RNA"]][[]]$Gene_ID %in% g2m.genes),"ENSG"]


    seurat <- CellCycleScoring(seurat, s.features = s.genes_ENSG, g2m.features = g2m.genes_ENSG, set.ident = TRUE)

    umap_cellcycle <- DimPlot(seurat, group.by = "Phase")
    ggsave(umap_cellcycle, filename = paste0(outdir,"pool_umap_cell_cycle.png"))

    lib_mt <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Pool", cols = alpha(PoolColors, 0.25))
    ggsave(lib_mt, filename = paste0(outdir,"lib_mt.png"))

    lib_genes <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Pool", cols = alpha(PoolColors, 0.25))
    ggsave(lib_genes, filename = paste0(outdir,"lib_genes.png"))

    nebulosa_mt_umap <- plot_density(seurat, "percent.mt", pal = "plasma")
    ggsave(nebulosa_mt_umap, filename = paste0(outdir,"mt_percent_umap.png"))

    nebulosa_rb_umap <- plot_density(seurat, "percent.rb", pal = "plasma")
    ggsave(nebulosa_rb_umap, filename = paste0(outdir,"rb_percent_umap.png"))

    umap_Time <- DimPlot(seurat, group.by = "Time")
    ggsave(umap_Time, filename = paste0(outdir,"Time_umap.png"))

    umap_Pool <- DimPlot(seurat, group.by = "Pool")
    ggsave(umap_Pool, filename = paste0(outdir,"pool_umap.png"))

    umap_location <- DimPlot(seurat, group.by = "MULTI_ID")
    ggsave(umap_location, filename = paste0(outdir,"location_umap.png"))

    mt_umap <- FeaturePlot(seurat, features = "percent.mt")
    ggsave(mt_umap, filename = paste0(outdir,"mt_percent_umap_seurat.png"))

    assignment_umap <- DimPlot(seurat, group.by = "Final_Assignment")
    ggsave(assignment_umap, filename = paste0(outdir,"Final_Assignment_umap_seurat.png"))


    ### Separate the seurat image by site and time and save separately
    seurat@meta.data$Location_Time <- gsub(" ", "_", paste0(gsub("\\d","",seurat@meta.data$MULTI_ID), "_", seurat@meta.data$Time))

    for (location in unique(seurat@meta.data$Location_Time)){
        saveRDS(subset(seurat, subset = Location_Time == location), paste0(outdir,location,"_seurat.rds"))
    }
    
}
