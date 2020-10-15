.libPaths( c("/share/ClusterShare/software/contrib/sccg/miniconda3/envs/EdgeR/lib/R/library" , .libPaths()))
print(.libPaths())

library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(ggpubr)
library(edgeR)
library(Seurat)

##### Read in the arguments #####
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


if (stage[1] == "PrePost"){
    ##### Read in Data #####
    data <- readRDS("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Seurat/Harmony_DayPatient/NoOutliers_DayPatient_harmony.rds")
    clusters <- read_delim("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Seurat/Multi_Resolution_Harmony_DayPatient/Identities_Resolution_0.02.txt", delim = "\t")
    data@meta.data[,c(colnames(clusters))] <- clusters

##### Make counts DGE lists and metadata matrices for each sub cell type #####
    cluster <- unique(data@meta.data$Resolution_0.02)[sge]
    print(cluster)
    SUB <- subset(data ,subset = Resolution_0.02 == cluster)
    counts <- DGEList(counts=SUB[["RNA"]]@counts, genes=rownames(SUB[["RNA"]]@counts))
    metadata <- cbind(row.names(SUB@meta.data),
                            as.numeric(SUB@meta.data$percent.mt), 
                            as.numeric(SUB@meta.data$percent.rb),
                            as.character(SUB@meta.data$Day),
                            as.character(SUB@meta.data$FinalAssignment))

    print(head(metadata))
    colnames(metadata) <- c("Barcode",
                                        "PercentMito",
                                        "PercentRibo",
                                        "Day",
                                        "Individual")
    print(head(metadata))

    print("Metadata dim 1")
    print(dim(metadata))
    print("Counts dim 1")
    print(dim(counts$counts))

    ##### Remove genes with less than 1 RPM in 25% of cells #####
    print("Starting RPM Analysis")
    counts$genes <- counts$genes[which(rowSums((counts$counts/(rowSums(counts$counts)/1000000))>1) > ncol(counts$counts)*0.25),]
    counts$counts <- counts$counts[which(rowSums((counts$counts/(rowSums(counts$counts)/1000000))>1) > ncol(counts$counts)*0.25),]
    
    print("RPM analysis done")

    print("Counts dim 2")
    print(dim(counts$counts))

    ##### Calculate Detection Rate #####
    print("starting calculate detection rate")
    metadata <- as.data.frame(cbind(metadata,DetRate = (colSums(counts$counts > 0)/nrow(counts$counts))))

    ##### Make bins for detection rate #####
    print("Calculating bins for detection rate")
    metadata$MtPercentBins <- cut(as.numeric(as.character(metadata$PercentMito)),breaks = 50)
    metadata$DetRateBins <- cut(as.numeric(as.character(metadata$DetRate)),breaks = 10)


    saveRDS(metadata, file = paste0(outdir,cluster,"_CellSubtypeMetadata.rds"))
    print("Metadata dim 2")
    dim(metadata)


    print("Finished calculating detection rate")


    ##### Normalize by Cell Size #####
    counts <- calcNormFactors(counts)
    print("Counts dim 3")
    print(dim(counts$counts))

    saveRDS(counts, file = paste0(outdir,cluster,"_CellSubtypeFilteredCounts.rds"))


} else if (stage[1] == "PrePost_DE"){
    metadataFileList <- list.files(paste0(dir,"output/EdgeR/PrePost/"), pattern = "Metadata.rds")
    countsFileList <- list.files(paste0(dir,"output/EdgeR/PrePost/"), pattern = "FilteredCounts.rds")

    metadata <- readRDS(file = paste0(dir,"output/EdgeR/PrePost/",metadataFileList[sge]))
    counts <- readRDS(file = paste0(dir,"output/EdgeR/PrePost/",countsFileList[sge]))
    counts$genes <- data.frame(Symbol=counts$genes)

    clusters <- gsub("_CellSubtypeFilteredCounts.rds","", countsFileList[sge])

    print(head(metadata))
    metadata <- as.data.frame(metadata)
    print(head(metadata))

    #### Use EdgeR to run differential expression
    ### Pre-Post DE ###
    ## Create Group variables
    print("DE genes with detection rate as covariate")
    
    individual <- factor(metadata$Individual)
    detrate <- metadata$DetRateBins
    day <- factor(metadata$Day, levels = c("Day0","Day4"))
    designTimeDet <- model.matrix(~individual + detrate + day)

    saveRDS(designTimeDet, file = paste0(outdir,clusters,"_designDay_DetIndCovar.rds"))

    counts <- estimateDisp(counts, designTimeDet)

    print(counts$common.dispersion)
    saveRDS(counts, file = paste0(outdir,clusters,"_countsDet.rds"))

    fitTime <- glmQLFit(counts, designTimeDet)
    saveRDS(fitTime, file = paste0(outdir,clusters,"_fitTimeDet.rds"))

    qlfTime <- glmQLFTest(fitTime)
    saveRDS(qlfTime, file = paste0(outdir,clusters,"_qlfTimeDet.rds"))

    write.table(topTags(qlfTime,n=Inf, p=0.05), file = paste0(outdir,clusters,"_WithDetRate.txt"), sep = "\t", quote = FALSE)


} else if (stage[1] == "PrePost_All"){
    ##### Read in Data #####
    data <- readRDS("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Seurat_w_Hash/SCTnormalization_no_covatiates/seurat_SCT_noCov.rds")
    data@meta.data$Site <- gsub("\\d","",data@meta.data$MULTI_ID)

##### Make counts DGE lists and metadata matrices for each sub cell type #####
    counts <- DGEList(counts=data[["SCT"]]@counts, genes=rownames(data[["SCT"]]@counts))
    metadata <- cbind(row.names(data@meta.data),
                            as.numeric(data@meta.data$percent.mt), 
                            as.numeric(data@meta.data$percent.rb),
                            as.character(data@meta.data$Site),
                            as.character(data@meta.data$FinalAssignment),
                            as.character(data@meta.data$Time))

    print(head(metadata))
    colnames(metadata) <- c("Barcode",
                            "PercentMito",
                            "PercentRibo",
                            "Site",
                            "CellLine",
                            "Time")
    print(head(metadata))

    print("Metadata dim 1")
    print(dim(metadata))
    print("Counts dim 1")
    print(dim(counts$counts))

    ##### Remove genes with less than 1 RPM in 25% of cells #####
    print("Starting RPM Analysis")
    counts$genes <- counts$genes[which(rowSums((counts$counts/(rowSums(counts$counts)/1000000))>1) > ncol(counts$counts)*0.25),]
    counts$counts <- counts$counts[which(rowSums((counts$counts/(rowSums(counts$counts)/1000000))>1) > ncol(counts$counts)*0.25),]
    
    print("RPM analysis done")

    print("Counts dim 2")
    print(dim(counts$counts))

    ##### Calculate Detection Rate #####
    print("starting calculate detection rate")
    metadata <- as.data.frame(cbind(metadata,DetRate = (colSums(counts$counts > 0)/nrow(counts$counts))))

    ##### Make bins for detection rate #####
    print("Calculating bins for detection rate")
    metadata$MtPercentBins <- cut(as.numeric(as.character(metadata$PercentMito)),breaks = 50)
    metadata$RbPercentBins <- cut(as.numeric(as.character(metadata$PercentRibo)), breaks = 50)
    metadata$DetRateBins <- cut(as.numeric(as.character(metadata$DetRate)),breaks = 10)


    saveRDS(metadata, file = paste0(outdir,"AllCells_metadata.rds"))
    print("Metadata dim 2")
    dim(metadata)


    print("Finished calculating detection rate")


    ##### Normalize by Cell Size #####
    counts <- calcNormFactors(counts)
    print("Counts dim 3")
    print(dim(counts$counts))

    saveRDS(counts, file = paste0(outdir,"AllCells_counts.rds"))


} else if (stage[1] == "PrePost_All_DE"){
    metadata <- readRDS(file = paste0(dir,"output/EdgeR/PrePost_All/AllCells_metadata.rds"))
    counts <- readRDS(file = paste0(dir,"output/EdgeR/PrePost_All/AllCells_counts.rds"))
    counts$genes <- data.frame(Symbol=counts$genes)

    print(head(metadata))
    metadata <- as.data.frame(metadata)
    print(head(metadata))

    #### Use EdgeR to run differential expression
    ### Pre-Post DE ###
    ## Create Group variables
    print("DE genes with detection rate as covariate")
    
    site <- factor(metadata$Site)
    detrate <- metadata$DetRateBins
    mt_pct <- metadata$MtPercentBins
    rb_pct <- metadata$RbPercentBins
    cell_line <- factor(metadata$CellLine)
    day <- factor(metadata$Time, levels = c("Day_0","Day_4"))
    designTime <- model.matrix(~ site + detrate + mt_pct + rb_pct + cell_line + day)

    saveRDS(designTime, file = paste0(outdir,"AllCells_designDay_SiteDetIndMtRbCovar.rds"))

    counts <- estimateDisp(counts, designTime)

    print(counts$common.dispersion)
    saveRDS(counts, file = paste0(outdir,"AllCells_counts.rds"))

    fitTime <- glmQLFit(counts, designTime)
    saveRDS(fitTime, file = paste0(outdir,"AllCells_fitTime.rds"))

    qlfTime <- glmQLFTest(fitTime)
    saveRDS(qlfTime, file = paste0(outdir,"AllCells_qlfTime.rds"))

    write.table(topTags(qlfTime,n=Inf, p=0.05), file = paste0(outdir,"AllCells_DE.txt"), sep = "\t", quote = FALSE)





} else if (stage[1] == "PrePost_CellLines"){
    ##### Read in Data #####
    data <- readRDS("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Seurat_w_Hash/SCTnormalization_no_covatiates/seurat_SCT_noCov.rds")
    data@meta.data$Site <- gsub("\\d","",data@meta.data$MULTI_ID)

##### Make counts DGE lists and metadata matrices for each sub cell type #####
    line <- unique(data@meta.data$FinalAssignment)[sge]
    SUB <- subset(data ,subset = FinalAssignment == line)
    counts <- DGEList(counts=SUB[["SCT"]]@counts, genes=rownames(SUB[["SCT"]]@counts))
    metadata <- cbind(row.names(SUB@meta.data),
                            as.numeric(SUB@meta.data$percent.mt), 
                            as.numeric(SUB@meta.data$percent.rb),
                            as.character(SUB@meta.data$Time),
                            as.character(SUB@meta.data$Site),
                            as.character(SUB@meta.data$FinalAssignment))

    print(head(metadata))
    colnames(metadata) <- c("Barcode",
                                        "PercentMito",
                                        "PercentRibo",
                                        "Day",
                                        "Site",
                                        "Individual")
    print(head(metadata))

    print("Metadata dim 1")
    print(dim(metadata))
    print("Counts dim 1")
    print(dim(counts$counts))

    ##### Remove genes with less than 1 RPM in 25% of cells #####
    print("Starting RPM Analysis")
    counts$genes <- counts$genes[which(rowSums((counts$counts/(rowSums(counts$counts)/1000000))>1) > ncol(counts$counts)*0.25),]
    counts$counts <- counts$counts[which(rowSums((counts$counts/(rowSums(counts$counts)/1000000))>1) > ncol(counts$counts)*0.25),]
    
    print("RPM analysis done")

    print("Counts dim 2")
    print(dim(counts$counts))

    ##### Calculate Detection Rate #####
    print("starting calculate detection rate")
    metadata <- as.data.frame(cbind(metadata,DetRate = (colSums(counts$counts > 0)/nrow(counts$counts))))

    ##### Make bins for detection rate #####
    print("Calculating bins for detection rate")
    metadata$MtPercentBins <- cut(as.numeric(as.character(metadata$PercentMito)),breaks = 50)
    metadata$RbPercentBins <- cut(as.numeric(as.character(metadata$PercentRibo)), breaks = 50)
    metadata$DetRateBins <- cut(as.numeric(as.character(metadata$DetRate)),breaks = 10)


    saveRDS(metadata, file = paste0(outdir,line,"_metadata.rds"))
    print("Metadata dim 2")
    dim(metadata)


    print("Finished calculating detection rate")


    ##### Normalize by Cell Size #####
    counts <- calcNormFactors(counts)
    print("Counts dim 3")
    print(dim(counts$counts))

    saveRDS(counts, file = paste0(outdir,line,"_counts.rds"))



} else if (stage[1] == "PrePost_CellLines_DE"){
    metadataFileList <- list.files(paste0(dir,"output/EdgeR/PrePost_CellLines/"), pattern = "metadata.rds")
    countsFileList <- list.files(paste0(dir,"output/EdgeR/PrePost_CellLines/"), pattern = "counts.rds")

    metadata <- readRDS(file = paste0(dir,"output/EdgeR/PrePost_CellLines/",metadataFileList[sge]))
    counts <- readRDS(file = paste0(dir,"output/EdgeR/PrePost_CellLines/",countsFileList[sge]))
    counts$genes <- data.frame(Symbol=counts$genes)

    line <- gsub("_counts.rds","", countsFileList[sge])

    print(head(metadata))
    metadata <- as.data.frame(metadata)
    print(head(metadata))

    #### Use EdgeR to run differential expression
    ### Pre-Post DE ###
    ## Create Group variables
    print("DE genes with detection rate as covariate")
    
    site <- factor(metadata$Site)
    detrate <- metadata$DetRateBins
    mt_pct <- metadata$MtPercentBins
    rb_pct <- metadata$RbPercentBins
    day <- factor(metadata$Day, levels = c("Day_0","Day_4"))
    designTime <- model.matrix(~ site + detrate + mt_pct + rb_pct + day)


    saveRDS(designTime, file = paste0(outdir,line,"_designDay_DetIndCovar.rds"))

    counts <- estimateDisp(counts, designTime)

    print(counts$common.dispersion)
    saveRDS(counts, file = paste0(outdir,line,"_countsDet.rds"))

    fitTime <- glmQLFit(counts, designTime)
    saveRDS(fitTime, file = paste0(outdir,line,"_fitTimeDet.rds"))

    qlfTime <- glmQLFTest(fitTime)
    saveRDS(qlfTime, file = paste0(outdir,line,"_qlfTimeDet.rds"))

    write.table(topTags(qlfTime,n=Inf, p=0.05), file = paste0(outdir,line,"_WithDetRate.txt"), sep = "\t", quote = FALSE)




}