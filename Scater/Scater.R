##### Analysis of Immunotherapy 10X scRNA-seq patients #####
## Read in Libraries
library(dplyr)
library(scater) 
library(ggplot2)
library(gridExtra,
        tidyr)
library(Seurat)
library(mvoutlier)
library(scran)
library(SingleCellExperiment)

source("/directflow/SCCGGroupShare/projects/DrewNeavin/tools/Rfunctions/Scater/ScaterFigures.R")
source("/directflow/SCCGGroupShare/projects/DrewNeavin/tools/Rfunctions/Scater/NormalizeScaterIdentifyOutliers.R")

## Read in arguments from data submission

args <- commandArgs(TRUE)
arguments <- read.table(args, header = F)

datadir <- arguments[1,]
outdir <- arguments[2,]
stage <- arguments[3,]
numsamp <- arguments[4,]

dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"


pools <- dir("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/200128_A00152_0196_BH3HNFDSXY/GE/", pattern = "DRENEA")
dataDIRS <- paste0("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/200128_A00152_0196_BH3HNFDSXY/GE/", pools, "/outs/filtered_feature_bc_matrix/")
GeneConversion <- read.delim("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/200128_A00152_0196_BH3HNFDSXY/GE/DRENEA_1/outs/filtered_feature_bc_matrix/features.tsv.gz", header = F)
RbGeneList <- read.delim(file = "/directflow/SCCGGroupShare/projects/DrewNeavin/References/RibosomalGeneList.txt",header = F)


print("Loading in data for cell cycle determination")
hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))

if (stage[1] == "CellCycle"){
    print("Starting CellCycle analysis")##### Determine Cell Cycle stage for each cell #####
    ### Set up the directory
    dataDIR <- as.character(dataDIRS[as.numeric(as.character(numsamp))])
    print(dataDIR)
    pool <- pools[as.numeric(as.character(numsamp))]
    print(pool)
    outdir <- paste0(outdir,pool,"/")
    dir.create(outdir)

    ##### Read in the Data #####
    counts <- Read10X(dataDIR) ## files must be unzipped for v2.0 of cellranger and gzipped for v3.0 of cellranger
    counts <- as.matrix(counts)
    sce <- SingleCellExperiment( assays = list(counts = counts))
    colnames(sce) <- paste0(pool,"_",colnames(sce))
    saveRDS(sce, paste0(outdir, pool,"_sce.rds"))

    print("Starting cell cycle determination")
    countsENSG <- counts(sce) # Create a counts matrix that has the ENSG gene clasifiers so can determine cell cycle phase
    print(dim(countsENSG))

    ### Reassign the row names to be ENSG IDs ###
    row.names(countsENSG) <- GeneConversion$V1
    ### Run the cell cycle identification ###
    assigned <- cyclone(countsENSG, pairs=hs.pairs)  #Note, this takes hours
    table(assigned$phases)
    write.table(assigned, file = paste0(outdir,"CellCycleProportions.txt"), quote = F, sep = "\t") #Save so that can read in and don't have to wait to recompute again

} else if (stage[1] == "merged"){

    print("Starting merged analysis")
    ##### Read in Data #####
    # pool <- pools[as.numeric(as.character(numsamp))]
    # print(pool)
    CellCycleDirList <- dir(path = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/ScaterQC/CellCycle/", pattern = "DRENEA")
    CellCycleList <- lapply(CellCycleDirList, FUN = function(x){read.table(paste0("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/ScaterQC/CellCycle/",x,"/CellCycleProportions.txt"))})
    CellCycles <- do.call(rbind, CellCycleList)

    SCEobjectList <- lapply(CellCycleDirList, FUN = function(x){readRDS(paste0("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/ScaterQC/CellCycle/",x,"/",x,"_sce.rds"))})
    sce <- do.call(cbind, SCEobjectList)

    sce$CellCycle <- CellCycles$phases
    saveRDS(sce, paste0(outdir,"/AllCellsSCEobject.rds"))

    ##### Make a list of control gene #####
    MtGenes <- rownames(sce)[grep("^MT-",rownames(sce))] # Make a list of mitochondrial genes
    RbGenes <- rownames(sce)[rownames(sce) %in% RbGeneList$V1]
    feature_controls <- list(mito = MtGenes, ribo = RbGenes)

    ##### Normalize data and identify outliers #####
    print("Starting normalization and 'outlier' identification")
    SCEallGenes <- NormalizeScaterIdentifyOutliers(ScaterObj = sce, SubType = "AllExpGene", feature_controls = feature_controls)
    saveRDS(SCEallGenes, file = paste0(outdir,"sce_AllGenesNormalized.rds"))
    SCEnoRbMt <- NormalizeScaterIdentifyOutliers(ScaterObj = sce, SubType = "NoRbMt", feature_controls = feature_controls)
    saveRDS(SCEnoRbMt, file = paste0(outdir,"sce_NoRbMtNormalized.rds"))

    # SCEallGenes <- readRDS(file = paste0(outdir,"sce_AllGenesNormalized.rds"))
    # SCEnoRbMt <- readRDS(file = paste0(outdir,"sce_NoRbMtNormalized.rds"))


    dir.create(paste0(outdir,"AllGenesRawCounts/"))
    dir.create(paste0(outdir,"AllGenesLogCounts/"))
    dir.create(paste0(outdir,"NoRbMtRawCounts/"))
    dir.create(paste0(outdir,"NoRbMtLogCounts/"))

    ScaterFigures(ScaterObj = SCEallGenes, 
                Assay = "counts",
                RbMtPresent = TRUE,
                outdir = paste0(outdir,"AllGenesRawCounts/AllGenesRawCounts"),
                colour1 = "CellCycle", 
                colour2 = "outlier",
                colour3 = "Pool")

    ScaterFigures(ScaterObj = SCEallGenes, 
                    Assay = "logcounts",
                    RbMtPresent = TRUE,
                    outdir = paste0(outdir,"AllGenesLogCounts/AllGenesLogCounts"),
                    colour1 = "CellCycle", 
                    colour2 = "outlier",
                    colour3 = "Pool")

    ScaterFigures(ScaterObj = SCEnoRbMt, 
                    Assay = "counts",
                    RbMtPresent = FALSE,
                    outdir = paste0(outdir,"NoRbMtRawCounts/NoRbMtRawCounts"),
                    colour1 = "CellCycle", 
                    colour2 = "outlier",
                    colour3 = "Pool")

    ScaterFigures(ScaterObj = SCEnoRbMt, 
                    Assay = "logcounts",
                    RbMtPresent = FALSE,
                    outdir = paste0(outdir,"NoRbMtLogCounts/NoRbMtLogCounts"),
                    colour1 = "CellCycle", 
                    colour2 = "outlier",
                    colour3 = "Pool")


    ### Remove greater than 3 median absolute deviations for the mitochondiral content
    mt_mad <- mad(SCEallGenes@meta.data$percent.mt)
    mt_low <- median(SCEallGenes@meta.data$percent.mt) - 3*mt_mad
    mt_high <- median(SCEallGenes@meta.data$percent.mt) + 3*mt_mad
    print("The mt mad is:")
    print(mt_mad)
    print("The lower bound for mt percent is:")
    print(mt_low)
    print("The upper bound for mt percent is:")
    print(mt_high)
    sce$mt_mad_outlier <- ifelse((mt_mad > mt_low & mt_mad < mt_high),FALSE, TRUE)

    percent.mito <- plotColData(object = SCEallGenes,
                    x = "Pool",
                    y = "pct_counts_mito",
                    colour_by = "outlier") +
                    geom_hline(yintercept=mt_low, linetype="dashed", color = "red") +
                    geom_hline(yintercept=mt_high, linetype="dashed", color = "red")
    ggsave(filename = paste0(outdir,"pctMitoViolinOutliersRemoved.png"), plot = percent.mito)


    rb_mad <- mad(SCEallGenes@meta.data$percent.rb)
    rb_low <- median(SCEallGenes@meta.data$percent.rb) - 3*rb_mad
    rb_high <- median(SCEallGenes@meta.data$percent.rb) + 3*rb_mad
    print("The rb mad is:")
    print(rb_mad)
    print("The lower bound for rb percent is:")
    print(rb_low)
    print("The upper bound for rb percent is:")
    print(rb_high)
    sce$rb_mad_outlier <- ifelse((rb_mad > rb_low & rb_mad < rb_high),FALSE, TRUE)

    percent.ribo <- plotColData(object = SCEallGenes,
                    x = "Pool",
                    y = "pct_counts_mito",
                    colour_by = "outlier") +
                    geom_hline(yintercept=rb_low, linetype="dashed", color = "red") +
                    geom_hline(yintercept=rb_high, linetype="dashed", color = "red")
    ggsave(filename = paste0(outdir,"pctMitoViolinOutliersRemoved.png"), plot = percent.mito)

} else if (stage[1] == "merged_post_scrublet"){
    


}