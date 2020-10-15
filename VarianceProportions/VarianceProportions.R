.libPaths("/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/generalR/lib/R/library")
library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(limma, lib.loc = "/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/generalR/lib/R/library")
library(variancePartition)
library(SingleCellExperiment)
library(scater)
library(Seurat)

##### Read in Arguments #####
print("Reading and assigning input arguments")
args <- commandArgs(TRUE)
arguments <- read.table(args, header = F)
stage <- arguments[1,]
outdir <- arguments[2,]
dir <- arguments[3,]
sge <- arguments[4,]
sge <- as.numeric(as.character(sge))
print(sge)
print(outdir)


if (stage[1]=="SCTnormalizedReplicate"){
    seurat <- readRDS("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Seurat_w_Hash/SCTnormalization_no_covatiates/seurat_SCT_noCov.rds")
    seurat@meta.data$Site <- gsub("\\d","",seurat@meta.data$MULTI_ID)

    ### Pull out the expression and metadata ###
    expression <- seurat[["SCT"]]@scale.data
    meta <- seurat@meta.data[,c("FinalAssignment","Site","MULTI_ID","Time","Pool","percent.mt","percent.rb")]

    ### Identify the variables to test for explaining gene variance ###
    form <- ~ (1|FinalAssignment) + (1|Site) + (1|MULTI_ID) + (1|Time) + (1|Pool) + percent.mt + percent.rb

    ### Test the variance explained for each of the genes by the covariates ###
    varPart <- fitExtractVarPartModel(expression, form, meta)

    ### Sort variables (i.e. columns) by median fraction of variance explained
    vp <- sortCols(varPart)
    vp <- vp[order(vp[,1]), decreasing=TRUE]
    write_delim(vp, paste0(outdir,"VarianceContributionsResults.txt"), delim = "\t")

    ### Bar plot of variance fractions for the first 10 genes ###
    TopVarianceExplainedGenes <- plotPercentBars(vp[1:10,])
    ggsave(TopVarianceExplainedGenes, filename = paste0(outdir,"TopVarianceExplainedGenes.png"))
    
    ### violin plot of contribution of each variable to total variance ###
    p_violin <-mplotVarPart(vp)
    ggsave(p_violin, filename = paste0(outdir,"ViolingVariance.png"))


} else if (stage[1]=="SCTnormalized") {
    seurat <- readRDS("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Seurat_w_Hash/SCTnormalization_no_covatiates/seurat_SCT_noCov.rds")
    seurat@meta.data$Site <- gsub("\\d","",seurat@meta.data$MULTI_ID)

    ### Pull out the expression and metadata ###
    expression <- seurat[["SCT"]]@scale.data
    meta <- seurat@meta.data[,c("FinalAssignment","Site","Time","Pool","percent.mt","percent.rb")]

    ### Identify the variables to test for explaining gene variance ###
    form <- ~ (1|FinalAssignment) + (1|Site)  + (1|Time) + (1|Pool) + percent.mt + percent.rb

    ### Test the variance explained for each of the genes by the covariates ###
    varPart <- fitExtractVarPartModel(expression, form, meta)

    ### Sort variables (i.e. columns) by median fraction of variance explained
    vp <- sortCols(varPart)
    vp <- vp[order(vp[,1]), decreasing=TRUE]
    write_delim(vp, paste0(outdir,"VarianceContributionsResults.txt"), delim = "\t")

    ### Bar plot of variance fractions for the first 10 genes ###
    TopVarianceExplainedGenes <- plotPercentBars(vp[1:10,])
    ggsave(TopVarianceExplainedGenes, filename = paste0(outdir,"TopVarianceExplainedGenes.png"))
    
    ### violin plot of contribution of each variable to total variance ###
    p_violin <-mplotVarPart(vp)
    ggsave(p_violin, filename = paste0(outdir,"ViolingVariance.png"))


} else if (stage[1] == "scaterVariance" ) {
    seurat <- readRDS("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Seurat_w_Hash/SCTnormalization_no_covatiates/seurat_SCT_noCov.rds")
    seurat@meta.data$Site <- gsub("\\d","",seurat@meta.data$MULTI_ID)
    sce <- SingleCellExperiment(assays = list(counts = seurat[["SCT"]]@counts, logcounts = seurat[["SCT"]]@data, normcounts = seurat[["SCT"]]@scale.data), colData = seurat@meta.data)


    vars <- getVarianceExplained(sce, variables=c("FinalAssignment","Site","Time","percent.mt","percent.rb"), exprs_values = "logcounts")
    head(vars)

    vars <- sortCols(vars)
    vars <- vars[order(vars[,1], decreasing=TRUE),]
    write.table(as.data.frame(vars), paste0(outdir,"VarianceContributionsResults.txt"), sep = "\t")

    ### Bar plot of variance fractions for the first 10 genes ###
    TopVarianceExplainedGenes <- plotExplanatoryVariables(vars)
    ggsave(TopVarianceExplainedGenes, filename = paste0(outdir,"TopVarianceExplainedGenes.png"))

    ### violin plot of contribution of each variable to total variance ###
    vars <- read.table(paste0(outdir,"VarianceContributionsResults.txt"), sep = "\t", header = TRUE)
    vars$gene <- rownames(vars)
    vars_long <- pivot_longer(vars, cols = c("percent.rb", "Site", "FinalAssignment", "Time", "percent.mt"), names_to = "covariate", values_to = "percent_variance")
    
    p_density <- ggplot(vars_long, aes(x = log10(percent_variance), color = covariate)) +
        geom_density() +
        theme_classic()
    ggsave(p_density, filename = paste0(outdir,"Density.png"))

    p_violin <-ggplot(vars_long, aes(covariate, percent_variance)) +
        geom_violin()
    ggsave(p_violin, filename = paste0(outdir,"ViolingVariance.png"))

    p_violin <-ggplot(vars_long, aes(covariate, log(percent_variance))) +
        geom_violin()
    ggsave(p_violin, filename = paste0(outdir,"ViolingVarianceLog.png"))


} else if (stage[1] == "scaterVarianceALL" ) {
    seurat <- readRDS("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Seurat_w_Hash/SCTnormalization_no_covatiates/seurat_SCT_noCov.rds")
    seurat@meta.data$Site <- gsub("\\d","",seurat@meta.data$MULTI_ID)
    sce <- SingleCellExperiment(assays = list(counts = seurat[["SCT"]]@counts, logcounts = seurat[["SCT"]]@data, normcounts = seurat[["SCT"]]@scale.data), colData = seurat@meta.data)


    vars <- getVarianceExplained(sce, variables=c("FinalAssignment","Site","Time","Pool","MULTI_ID","percent.mt","percent.rb"), exprs_values = "normcounts")
    head(vars)

    vars <- sortCols(vars)
    vars <- vars[order(vars[,1], decreasing=TRUE),]
    write.table(as.data.frame(vars), paste0(outdir,"VarianceContributionsResults.txt"), sep = "\t")

    ### Bar plot of variance fractions for the first 10 genes ###
    TopVarianceExplainedGenes <- plotExplanatoryVariables(vars)
    ggsave(TopVarianceExplainedGenes, filename = paste0(outdir,"TopVarianceExplainedGenes.png"))

    ### violin plot of contribution of each variable to total variance ###
    vars <- read.table(paste0(outdir,"VarianceContributionsResults.txt"), sep = "\t", header = TRUE)
    vars$gene <- rownames(vars)
    vars_long <- pivot_longer(vars, cols = c("FinalAssignment","Site","Time","Pool","MULTI_ID","percent.mt","percent.rb"), names_to = "covariate", values_to = "percent_variance")
    
    p_violin <-ggplot(vars_long, aes(covariate, percent_variance)) +
        geom_violin()
    ggsave(p_violin, filename = paste0(outdir,"ViolingVariance.png"))

    p_violin <-ggplot(vars_long, aes(covariate, log(percent_variance))) +
        geom_violin()
    ggsave(p_violin, filename = paste0(outdir,"ViolingVarianceLog.png"))



} else if (stage[1] == "scaterVariance_RbMtSiteCov" ) {
    seurat <- readRDS("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Seurat_w_Hash/SCTnormalization_SiteMtRbRegressed/seurat_SCT_SiteMtRbRegressed.rds")
    seurat@meta.data$Site <- gsub("\\d","",seurat@meta.data$MULTI_ID)
    sce <- as.SingleCellExperiment(seurat)

    vars <- getVarianceExplained(sce, 
        variables=c("FinalAssignment","Site","Time","percent.mt","percent.rb"), exprs_values = "logcounts")
    head(vars)

    vars <- sortCols(vars)
    vars <- vars[order(vars[,1], decreasing=TRUE),]
    write.table(as.data.frame(vars), paste0(outdir,"VarianceContributionsResults.txt"), sep = "\t")

    ### Bar plot of variance fractions for the first 10 genes ###
    TopVarianceExplainedGenes <- plotExplanatoryVariables(vars)
    ggsave(TopVarianceExplainedGenes, filename = paste0(outdir,"TopVarianceExplainedGenes.png"))

    ### violin plot of contribution of each variable to total variance ###
    vars <- read.table(paste0(outdir,"VarianceContributionsResults.txt"), sep = "\t", header = TRUE)
    vars$gene <- rownames(vars)
    vars_long <- pivot_longer(vars, cols = c("percent.rb", "FinalAssignment", "Time", "percent.mt"), names_to = "covariate", values_to = "percent_variance")
    
    p_violin <-ggplot(vars_long, aes(covariate, percent_variance)) +
        geom_violin()
    ggsave(p_violin, filename = paste0(outdir,"ViolingVariance.png"))

    p_violin <-ggplot(vars_long, aes(covariate, log(percent_variance))) +
        geom_violin()
    ggsave(p_violin, filename = paste0(outdir,"ViolingVariance.png"))


} else if (stage[1] == "scaterVariance_IntegratedSite" ) {
    seurat <- readRDS("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Seurat_w_Hash/Integration_Site/NoOutliers_Site_SeuratIntegratee.rds")
    seurat@meta.data$Site <- gsub("\\d","",seurat@meta.data$MULTI_ID)
    sce <- as.SingleCellExperiment(seurat)

    vars <- getVarianceExplained(sce, 
        variables=c("FinalAssignment","Site","Time","percent.mt","percent.rb"), exprs_values = "logcounts")
    head(vars)

    vars <- sortCols(vars)
    vars <- vars[order(vars[,1], decreasing=TRUE),]
    write.table(as.data.frame(vars), paste0(outdir,"VarianceContributionsResults.txt"), sep = "\t")

    ### Bar plot of variance fractions for the first 10 genes ###
    TopVarianceExplainedGenes <- plotExplanatoryVariables(vars)
    ggsave(TopVarianceExplainedGenes, filename = paste0(outdir,"TopVarianceExplainedGenes.png"))

    ### violin plot of contribution of each variable to total variance ###
    vars <- read.table(paste0(outdir,"VarianceContributionsResults.txt"), sep = "\t", header = TRUE)
    vars$gene <- rownames(vars)
    vars_long <- pivot_longer(vars, cols = c("percent.rb", "FinalAssignment", "Time", "percent.mt"), names_to = "covariate", values_to = "percent_variance")
    
    p_violin <-ggplot(vars_long, aes(covariate, percent_variance)) +
        geom_violin()
    ggsave(p_violin, filename = paste0(outdir,"ViolingVariance.png"))

    p_violin <-ggplot(vars_long, aes(covariate, log(percent_variance))) +
        geom_violin()
    ggsave(p_violin, filename = paste0(outdir,"ViolingVariance.png"))



} else if (stage[1] == "scaterVariance_SeparatedSite" ) {
    seurat <- readRDS("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Seurat_w_Hash/QC/seurat_norm_NoOutliers.rds")
    seurat@meta.data$Site <- gsub("\\d","", seurat@meta.data$MULTI_ID)
    print(head(seurat@meta.data$Site))
    site <- unique(seurat@meta.data$Site)[sge]
    print(paste0("The site is: ", site))
    seurat <- subset(seurat, subset = Site == site)
    seurat <- SCTransform(seurat, vars.to.regress = c("percent.mt","percent.rb"), verbose = TRUE, return.only.var.genes = FALSE)
    seurat <- RunPCA(seurat, npcs = 100)
    seurat <- RunUMAP(seurat, dims = 1:100, verbose = TRUE)
    seurat <- FindNeighbors(seurat, dims = 1:100, verbose = TRUE)

    sce <- SingleCellExperiment(assays = list(counts = seurat[["SCT"]]@counts, logcounts = seurat[["SCT"]]@data, normcounts = seurat[["SCT"]]@scale.data), colData = seurat@meta.data)

    ### Use the scater package to get the variance that is explained for each gene by each of the variables.
    vars <- getVarianceExplained(sce, 
        variables=c("FinalAssignment","Time","percent.mt","percent.rb","MULTI_ID"), exprs_values = "logcounts")
    head(vars)

    ### Sort the columns of the resulting dataframe (ie the covariate with the highest average contribution to gene variation across all genes is the first column and the one that has the lowest average contribution to the gene variation is the last column)
    vars <- sortCols(vars)

    ### Order the rows fo the dataframe so that the top row is the gene where the first has the highest amount of its variance described by the covariate in the first column
    vars <- vars[order(vars[,1], decreasing=TRUE),]

    ### Write the dataframe to file for future reference
    write.table(as.data.frame(vars), paste0(outdir,site,"_VarianceContributionsResults.txt"), sep = "\t")

    ### Density plot of the variance explained by each covariate for each gene ###
    TopVarianceExplainedGenes <- plotExplanatoryVariables(vars)
    ggsave(TopVarianceExplainedGenes, filename = paste0(outdir,site,"_TopVarianceExplainedGenes.png"))

    ### violin plot of contribution of each variable to total variance ###
    vars <- read.table(paste0(outdir,site,"_VarianceContributionsResults.txt"), sep = "\t", header = TRUE)
    vars$gene <- rownames(vars)
    vars_long <- pivot_longer(vars, cols = c("percent.rb", "FinalAssignment", "Time", "percent.mt", "MULTI_ID"), names_to = "covariate", values_to = "percent_variance")
    
    p_violin <-ggplot(vars_long, aes(covariate, percent_variance)) +
        geom_violin() +
        geom_boxplot(width=0.1)
    ggsave(p_violin, filename = paste0(outdir,site,"_ViolingVariance.png"))

    p_violin <-ggplot(vars_long, aes(covariate, log(percent_variance))) +
        geom_violin() +
        geom_boxplot(width=0.1)
    ggsave(p_violin, filename = paste0(outdir,site,"_Log_ViolingVariance.png"))



} else if (stage[1] == "scaterVariance_SeparatedSite_MtRb_regressed" ) {
    # seurat <- readRDS("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Seurat_w_Hash/QC/seurat_norm_NoOutliers.rds")
    # seurat@meta.data$Site <- gsub("\\d","", seurat@meta.data$MULTI_ID)
    # print(head(seurat@meta.data$Site))
    # site <- unique(seurat@meta.data$Site)[sge]
    # print(paste0("The site is: ", site))
    # seurat <- subset(seurat, subset = Site == site)
    # seurat <- SCTransform(seurat, vars.to.regress = c("percent.mt","percent.rb"), verbose = TRUE, return.only.var.genes = FALSE)
    # seurat <- RunPCA(seurat, npcs = 100)
    # seurat <- RunUMAP(seurat, dims = 1:100, verbose = TRUE)
    # seurat <- FindNeighbors(seurat, dims = 1:100, verbose = TRUE)

    # sce <- SingleCellExperiment(assays = list(counts = seurat[["SCT"]]@counts, logcounts = seurat[["SCT"]]@data, normcounts = seurat[["SCT"]]@scale.data), colData = seurat@meta.data)

    # vars <- getVarianceExplained(sce, 
    #     variables=c("FinalAssignment","Time","percent.mt","percent.rb","MULTI_ID"), exprs_values = "normcounts")
    # head(vars)

    # vars <- sortCols(vars)
    # vars <- vars[order(vars[,1], decreasing=TRUE),]
    # write.table(as.data.frame(vars), paste0(outdir,site,"_VarianceContributionsResults.txt"), sep = "\t")

    # ### Bar plot of variance fractions for the first 10 genes ###
    # TopVarianceExplainedGenes <- plotExplanatoryVariables(vars)
    # ggsave(TopVarianceExplainedGenes, filename = paste0(outdir,site,"_TopVarianceExplainedGenes.png"))

    ### violin plot of contribution of each variable to total variance ###
    site <- c("Melbourne","Sydney","Brisbane")[sge]
    vars <- read.table(paste0(outdir,site,"_VarianceContributionsResults.txt"), sep = "\t", header = TRUE)
    vars$gene <- rownames(vars)
    vars_long <- pivot_longer(vars, cols = c("percent.rb", "FinalAssignment", "Time", "percent.mt", "MULTI_ID"), names_to = "covariate", values_to = "percent_variance")
    
    p_violin <-ggplot(vars_long, aes(covariate, percent_variance)) +
        geom_violin() +
        geom_boxplot(width=0.1) +
        theme_classic()
    ggsave(p_violin, filename = paste0(outdir,site,"_ViolingVariance.png"))

    p_violin <-ggplot(vars_long, aes(covariate, log(percent_variance))) +
        geom_violin() +
        geom_boxplot(width=0.1) +
        theme_classic()
    ggsave(p_violin, filename = paste0(outdir,site,"_Log_ViolingVariance.png"))




} else if (stage[1] == "scaterVariance_SeparatedReplicate" ) {
    seurat <- readRDS("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Seurat_w_Hash/QC/seurat_norm_NoOutliers.rds")
    seurat@meta.data$Replicate <- paste0(seurat@meta.data$MULTI_ID, "_", seurat@meta.data$Time)
    print(head(seurat@meta.data$Replicate))
    rep <- unique(seurat@meta.data$Replicate)[sge]
    print(paste0("The sample replicate is ", rep))
    Idents(seurat) <- seurat@meta.data$Replicate
    seurat <- subset(seurat, idents = rep)
    seurat <- SCTransform(seurat, verbose = TRUE, return.only.var.genes = FALSE)
    seurat <- RunPCA(seurat, npcs = 100)
    seurat <- RunUMAP(seurat, dims = 1:100, verbose = TRUE)
    seurat <- FindNeighbors(seurat, dims = 1:100, verbose = TRUE)

    sce <- SingleCellExperiment(assays = list(counts = seurat[["SCT"]]@counts, logcounts = seurat[["SCT"]]@data, normcounts = seurat[["SCT"]]@scale.data), colData = seurat@meta.data)

    vars <- getVarianceExplained(sce, 
        variables=c("FinalAssignment","Time","percent.mt","percent.rb","nCount_RNA", "nFeature_RNA", "nCount_HTO", "nFeature_HTO"), exprs_values = "normcounts")
    head(vars)

    vars <- sortCols(vars)
    vars <- vars[order(vars[,1], decreasing=TRUE),]
    write.table(as.data.frame(vars), paste0(outdir,rep,"_VarianceContributionsResults.txt"), sep = "\t")

    ### Bar plot of variance fractions for the first 10 genes ###
    TopVarianceExplainedGenes <- plotExplanatoryVariables(vars)
    ggsave(TopVarianceExplainedGenes, filename = paste0(outdir,rep,"_TopVarianceExplainedGenes.png"))

    ### violin plot of contribution of each variable to total variance ###
    vars <- read.table(paste0(outdir,rep,"_VarianceContributionsResults.txt"), sep = "\t", header = TRUE)
    vars$gene <- rownames(vars)
    vars_long <- pivot_longer(vars, cols = c("FinalAssignment","Time","percent.mt","percent.rb","nCount_RNA", "nFeature_RNA", "nCount_HTO", "nFeature_HTO"), names_to = "covariate", values_to = "percent_variance")
    
    p_violin <-ggplot(vars_long, aes(covariate, percent_variance)) +
        geom_violin() +
        geom_boxplot(width=0.1)
    ggsave(p_violin, filename = paste0(outdir,rep,"_ViolingVariance.png"))

    p_violin <-ggplot(vars_long, aes(covariate, log(percent_variance))) +
        geom_violin() +
        geom_boxplot(width=0.1)
    ggsave(p_violin, filename = paste0(outdir,rep,"_Log_ViolingVariance.png"))




} else if (stage[1] == "scaterVariance_SeparatedSiteReplicate" ) {
    seurat <- readRDS("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Seurat_w_Hash/QC/seurat_norm_NoOutliers.rds")
    rep <- unique(seurat@meta.data$MULTI_ID)[sge]
    print(paste0("The sample replicate is ", rep))
    Idents(seurat) <- seurat@meta.data$MULTI_ID
    seurat <- subset(seurat, idents = rep)
    seurat <- SCTransform(seurat, verbose = TRUE, return.only.var.genes = FALSE)
    seurat <- RunPCA(seurat, npcs = 100)
    seurat <- RunUMAP(seurat, dims = 1:100, verbose = TRUE)
    seurat <- FindNeighbors(seurat, dims = 1:100, verbose = TRUE)

    sce <- SingleCellExperiment(assays = list(counts = seurat[["SCT"]]@counts, logcounts = seurat[["SCT"]]@data, normcounts = seurat[["SCT"]]@scale.data), colData = seurat@meta.data)

    vars <- getVarianceExplained(sce, 
        variables=c("FinalAssignment","Time","percent.mt","percent.rb","nCount_RNA", "nFeature_RNA", "nCount_HTO", "nFeature_HTO"), exprs_values = "normcounts")
    head(vars)

    vars <- sortCols(vars)
    vars <- vars[order(vars[,1], decreasing=TRUE),]
    write.table(as.data.frame(vars), paste0(outdir,rep,"_VarianceContributionsResults.txt"), sep = "\t")

    ### Bar plot of variance fractions for the first 10 genes ###
    TopVarianceExplainedGenes <- plotExplanatoryVariables(vars)
    ggsave(TopVarianceExplainedGenes, filename = paste0(outdir,rep,"_TopVarianceExplainedGenes.png"))

    ### violin plot of contribution of each variable to total variance ###
    # vars <- read.table(paste0(outdir,rep,"_VarianceContributionsResults.txt"), sep = "\t", header = TRUE)
    vars$gene <- rownames(vars)
    vars_long <- pivot_longer(vars, cols = c("FinalAssignment","Time","percent.mt","percent.rb","nCount_RNA", "nFeature_RNA", "nCount_HTO", "nFeature_HTO"), names_to = "covariate", values_to = "percent_variance")
    
    p_violin <-ggplot(vars_long, aes(covariate, percent_variance)) +
        geom_violin() +
        geom_boxplot(width=0.1)
    ggsave(p_violin, filename = paste0(outdir,rep,"_ViolingVariance.png"))

    p_violin <-ggplot(vars_long, aes(covariate, log(percent_variance))) +
        geom_violin() +
        geom_boxplot(width=0.1)
    ggsave(p_violin, filename = paste0(outdir,rep,"_Log_ViolingVariance.png"))


} else if (stage[1] == "scaterVariance_SeparatedSiteReplicate_RbMtRegressed" ) {
    seurat <- readRDS("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Seurat_w_Hash/QC/seurat_norm_NoOutliers.rds")
    rep <- unique(seurat@meta.data$MULTI_ID)[sge]
    print(paste0("The sample replicate is ", rep))
    Idents(seurat) <- seurat@meta.data$MULTI_ID
    seurat <- subset(seurat, idents = rep)
    seurat <- SCTransform(seurat, verbose = TRUE, vars.to.regress = c("percent.mt","percent.rb"), return.only.var.genes = FALSE)
    seurat <- RunPCA(seurat, npcs = 100)
    seurat <- RunUMAP(seurat, dims = 1:100, verbose = TRUE)
    seurat <- FindNeighbors(seurat, dims = 1:100, verbose = TRUE)

    sce <- SingleCellExperiment(assays = list(counts = seurat[["SCT"]]@counts, logcounts = seurat[["SCT"]]@data, normcounts = seurat[["SCT"]]@scale.data), colData = seurat@meta.data)

    vars <- getVarianceExplained(sce, 
        variables=c("FinalAssignment","Time","percent.mt","percent.rb","nCount_RNA", "nFeature_RNA", "nCount_HTO", "nFeature_HTO"), exprs_values = "normcounts")
    head(vars)

    vars <- sortCols(vars)
    vars <- vars[order(vars[,1], decreasing=TRUE),]
    write.table(as.data.frame(vars), paste0(outdir,rep,"_VarianceContributionsResults.txt"), sep = "\t")

    ### Bar plot of variance fractions for the first 10 genes ###
    TopVarianceExplainedGenes <- plotExplanatoryVariables(vars)
    ggsave(TopVarianceExplainedGenes, filename = paste0(outdir,rep,"_TopVarianceExplainedGenes.png"))

    ### violin plot of contribution of each variable to total variance ###
    # vars <- read.table(paste0(outdir,rep,"_VarianceContributionsResults.txt"), sep = "\t", header = TRUE)
    vars$gene <- rownames(vars)
    vars_long <- pivot_longer(vars, cols = c("FinalAssignment","Time","percent.mt","percent.rb","nCount_RNA", "nFeature_RNA", "nCount_HTO", "nFeature_HTO"), names_to = "covariate", values_to = "percent_variance")
    
    p_violin <-ggplot(vars_long, aes(covariate, percent_variance)) +
        geom_violin() +
        geom_boxplot(width=0.1)
    ggsave(p_violin, filename = paste0(outdir,rep,"_ViolingVariance.png"))

    p_violin <-ggplot(vars_long, aes(covariate, log(percent_variance))) +
        geom_violin() +
        geom_boxplot(width=0.1)
    ggsave(p_violin, filename = paste0(outdir,rep,"_Log_ViolingVariance.png"))




}