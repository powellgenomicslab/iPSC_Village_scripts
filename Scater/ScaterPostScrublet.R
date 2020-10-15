#### Analysis of Immunotherapy 10X scRNA-seq patients #####

## Read in Libraries
library(dplyr)
library(scater)
library(ggplot2)
library(gridExtra,
        tidyr)
library(Seurat)
library(mvoutlier)
library(scran)
library(rmutil)

source("directflow/SCCGGroupShare/projects/DrewNeavin/tools/Rfunctions/Scater/ScaterFigures.R")
source("directflow/SCCGGroupShare/projects/DrewNeavin/tools/Rfunctions/Scater/NormalizeScaterIdentifyOutliers.R")

## Read in arguments from data submission

args <- commandArgs(TRUE)
arguments <- read.table(args, header = F)

DataDIR <- arguments[1,]
outdir <- arguments[2,]
stage <- arguments[3,]
numsamp <- arguments[4,]


SCRUBLETdataDIR <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Scrublet/default0.85/"
names <- dir(SCRUBLETdataDIR, pattern = "DRENEA")

barcodeDIR <- dir("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/200128_A00152_0196_BH3HNFDSXY/GE/", pattern = "Pool_", full.names = TRUE)
print(names)
GeneConversion <- read.delim("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/200128_A00152_0196_BH3HNFDSXY/GE/DRENEA_1/outs/filtered_feature_bc_matrix/features.tsv", header = F)
RbGeneList <- read.delim(file = "/share/ScratchGeneral/drenea/References/RibosomalGeneList.txt",header = F)

##### Read in prepared SCE object and the predicted_doublet_mask.txt from Scrublet output and 
print("Loading in data")
### Read in SCE object that was originally used and low MT only (< 10%) as well as outliers removed - check for overlap with the doublets -> the cells that you get back will be those that DON't overlap!!!
pbmcAllGenesNormalizedNoOutliersLowMt <- readRDS("/share/ScratchGeneral/drenea/SarcomaIOproject/analysis/ScaterQC/mergedRemovedOutliers/pbmcAllGenesNormalizedNoOutliersLowMt.rds")

# ### Read in the SCE object that has low mitochondrial content only and has been QC metrics done - will be used for eventually feeding in to Scater after removal of doublets
pbmcLowMito <- readRDS(paste0(DataDIR,"/FilteredCounts/pbmcLowMito.rds"))
# pbmc <- readRDS(paste0("/share/ScratchGeneral/drenea/SarcomaIOproject/analysis/ScaterQC/merged/pbmc.rds"))

# ### Read in the Scrublet output data
SCRUBLETdata <- list()
SCRUBLETdata <- lapply(paste0(ScrubletDIRS,"/predicted_doublet_mask.txt"), read.table, sep = "\t")
names(SCRUBLETdata) <- names
SCRUBLETbarcodes <- list()
SCRUBLETbarcodes <- lapply(paste0(barcodeDIR,"/outs/filtered_feature_bc_matrix/barcodes.tsv"), read.table, sep = "\t")
names(SCRUBLETbarcodes) <- names

# ##### Merge Scrublet data and barcodes and then identify the number that were remaining in the 
SCRUBLETmerged <- mapply(cbind, SCRUBLETbarcodes, SCRUBLETdata, SIMPLIFY=F)

SCRUBLETmergedTable <- do.call(rbind, SCRUBLETmerged)


doublet <- lapply(SCRUBLETmerged, function(x) x[which(x[,2] == "True"),])
NOTdoublet <- lapply(SCRUBLETmerged, function(x) x[which(x[,2] == "False"),])

# doubletsNOToriginallyRMV <- lapply(doublet, function(x) x[x[,1] %in% colnames(pbmcAllGenesNormalizedNoOutliersLowMt),])

# print("The number of doublets detected for each batch are:")
# print(lapply(doublet, dim))

# print("The number of doublets that were not removed in the original analysis but were detected by Scrublet are:")
# print(lapply(doubletsNOToriginallyRMV, dim))

# saveRDS(doubletsNOToriginallyRMV, file = paste0(outdir,"/doubletsNOToriginallyRMV.rds"))

batch <- data.frame(matrix(ncol = 0, nrow = nrow((pbmcLowMito@colData))))
batch$barcode <- as.character(colnames(pbmcLowMito))
batch$batch <- pbmcLowMito@colData$BatchRun
SCRUBLETmergedTable$batch <- row.names(SCRUBLETmergedTable)
SCRUBLETmergedTable$batch <- gsub("\\.\\d+", "", SCRUBLETmergedTable$batch)
colnames(SCRUBLETmergedTable) <- c("barcode", "doublet", "batch")
SCRUBLETmerged$barcode <- as.character(SCRUBLETmergedTable$barcode)
print(head(SCRUBLETmergedTable))
print(anyDuplicated(SCRUBLETmergedTable$barcode))
print(head(batch))
batch$barcode <- substr(batch$barcode, 1, 16)
batch$barcode <- paste0(batch$barcode,"-1")
print(head(batch))
print(anyDuplicated(batch$batch))
batchmerged <- inner_join(batch,SCRUBLETmergedTable, by = c("barcode", "batch"))
print(dim(batch))
print(dim(SCRUBLETmergedTable))
print(dim(batchmerged))

print(head(batch))
print(head(batchmerged))

pbmcLowMito@colData$doublet <- batchmerged[,3]
print(head(colData(pbmcLowMito)))

# pbmcLowMito <- pbmc[,colData(pbmc)$pct_counts_mito < 10]
# saveRDS(pbmcLowMito, file = paste0(outdir,"/pbmcLowMitoWithDoubletInfo.rds"))

##### Make QC plots for doublet data #####
feature_controls <- readRDS(file = "/share/ScratchGeneral/drenea/SarcomaIOproject/analysis/ScaterQC/merged/RbMtFeatureControls.rds")

print("Starting normalization and 'outlier' identification")
pbmcLowMito <- NormalizeScaterIdentifyOutliers(ScaterObj = pbmcLowMito, SubType = "AllExpGene", feature_controls = feature_controls)
saveRDS(pbmcLowMito, file = paste0(outdir,"/pbmcLowMitoWithDoubletNormalized.rds"))

dir.create(paste0(outdir,"/LowMito"))

ScaterFigures(ScaterObj = pbmcLowMito, 
            Assay = "counts",
            RbMtPresent = TRUE,
            outdir = paste0(outdir,"/LowMito/LowMito"),
            colour1 = "doublet", 
            colour2 = "BatchRun",
            colour3 = "outlier")

##### Remove doublets from data #####
# pbmcLowMito <- readRDS(file = paste0(outdir,"/pbmcLowMitoWithDoubletNormalized.rds"))

print(head(colData(pbmcLowMito)))

pbmcLowMitoNoDoublet <- filter(pbmcLowMito, doublet == "False")

print("The dimension of the dataset before removing doublets:")
print(pbmcLowMito)

print("The dimension of the dataset after removing doublets:")
print(pbmcLowMitoNoDoublet)

##### Rerun Scater to detect outliers and generate new plots #####
feature_controls <- readRDS(file = "/share/ScratchGeneral/drenea/SarcomaIOproject/analysis/ScaterQC/merged/RbMtFeatureControls.rds")

print("Starting normalization and 'outlier' identification")
pbmcLowMitoNoDoublet <- NormalizeScaterIdentifyOutliers(ScaterObj = pbmcLowMitoNoDoublet, SubType = "AllExpGene", feature_controls = feature_controls)
saveRDS(pbmcLowMitoNoDoublet, file = paste0(outdir,"/pbmcLowMitoNoDoubletNormalized.rds"))

pbmcLowMitoNoDoublet <- readRDS(file = paste0(outdir,"/pbmcLowMitoNoDoubletNormalized.rds"))
dir.create(paste0(outdir,"/LowMitoNoDoublet/"))

ScaterFigures(ScaterObj = pbmcLowMitoNoDoublet, 
            Assay = "counts",
            RbMtPresent = TRUE,
            outdir = paste0(outdir,"/LowMitoNoDoublet/LowMitoNoDoublet"),
            colour1 = "CellCycle", 
            colour2 = "BatchRun",
            colour3 = "outlier")


