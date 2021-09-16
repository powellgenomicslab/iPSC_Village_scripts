library(Seurat)
library(data.table)
library(tradeSeq)

BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- 5 # use 2 cores

##### Setting up Directories
dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
datadir <- paste0(dir,"output/All_data_integrated_remove_bad/")
scvelo_dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/scVelo/velocyto/scvelo_umap/"
outdir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/scVelo/tradeseq_DE/"

dir.create(outdir, recursive = TRUE)


line_colors = c(FSA0006 = "#F79E29", MBE1006 = "#9B2C99", TOB0421 = "#35369C")


##### Read in data #####
seurat <- readRDS(paste0(datadir,"seurat_integrated_all_times_clustered.rds"))

latent_files <- list.files(scvelo_dir, pattern = ".csv")

meta_list <- lapply(latent_files, function(x){
    fread(paste0(scvelo_dir,x), sep = ",")
})


meta_df <- do.call(rbind, meta_list)


meta_latent <- data.table(meta_df[,latent_time])
colnames(meta_latent) <- "latent_time"
rownames(meta_latent) <- meta_df[,V1]

seurat <- AddMetaData(seurat, meta_latent)


seurat_sub <- seurat[rowSums(seurat[["SCT"]]@counts) > ncol(seurat)*0.1,]


icMat <- evaluateK(counts = seurat_sub[["SCT"]]@counts, pseudotime = data.frame(seurat@meta.data[,"latent_time"]), k = 3:10, 
                   nGenes = 200, verbose = T,  parallel=TRUE, BPPARAM = BPPARAM)



batch <- 

sce <- fitGAM(counts = seurat[["SCT"]]@counts, pseudotime = seurat@meta.data$latent_time, 
                 nknots = 6, verbose = FALSE, parallel=TRUE, BPPARAM = BPPARAM)
