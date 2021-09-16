library(Seurat)
library(SeuratDisk)


##### Setting up Directories
dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
datadir <- paste0(dir,"output/All_data_integrated_remove_bad/")
phatedir <- paste0(dir,"output/LineageTracing/Phate/")
outdir <- paste0(dir,"output/scVelo/preprocess/seurat/")

dir.create(outdir, recursive = TRUE)


##### Read in seurat object
seurat <- readRDS(paste0(datadir,"seurat_integrated_all_times_clustered.rds"))

write.csv(Cells(seurat), file = paste0(outdir, "cellID_obs.csv"), row.names = FALSE)
write.csv(Embeddings(seurat, reduction = "umap"), file = paste0(outdir,"cell_embeddings.csv"))
write.csv(seurat@meta.data, file = paste0(outdir,"metadata.csv"))

dir.create(paste0(outdir,"h5ad/"), recursive = TRUE)
SaveH5Seurat(seurat, filename = paste0(outdir,"seurat_integrated_all_times_clustered.h5Seurat"))

Convert(paste0(outdir,"seurat_integrated_all_times_clustered.h5Seurat"), dest = "h5ad")


##### Read in phate object
phate <- readRDS(paste0(phatedir,"phate_knn20.rds"))


write.csv(rownames(phate$embedding), file = paste0(outdir,"cellID_obs_phate.csv"), row.names = FALSE)
write.csv(phate$embedding, file = paste0(outdir,"cell_embeddings_phate.csv"))
