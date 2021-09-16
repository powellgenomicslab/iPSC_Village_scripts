library(tidyverse)
library(slingshot)
library(Seurat)


##### Set up directories #####
dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
datadir <- paste0(dir,"output/LineageTracing/Phate/")
outdir <- paste0(dir,"output/LineageTracing/Slingshot/")
dir.create(outdir, recursive = TRUE)



##### Read in Phate Data #####
phate <- readRDS(paste0(datadir,"phate_knn5.rds"))
seurat <- readRDS(paste0(dir,"output/Distribution_tests/seurat_integrated_all_times_clustered_1pct_expressing.rds"))



##### Slingshot using Clusters from PCA on PHATE embeddings #####
sling <- slingshot(phate$embedding, Idents(seurat))
saveRDS(sling, paste0(outdir,"slingshot_5knn.rds"))


png(paste0(outdir,"slingshot.png"))
plot(as.data.frame(phate$embedding), col = Idents(seurat), pch=20, asp = 1,
    xlab="PHATE 1", 
    ylab="PHATE 2",
    cex.main=1.25, 
    cex.lab=1, 
    cex.axis=1)
lines(sling, lwd=2, type = 'lineages', col = 'grey')
dev.off()   



png(paste0(outdir,"slingshot_curves.png"))
plot(as.data.frame(phate$embedding), col = Idents(seurat), pch=20, asp = 1,
    xlab="PHATE 1", 
    ylab="PHATE 2",
    cex.main=1.25, 
    cex.lab=1, 
    cex.axis=1)
lines(sling, lwd=2, col = 'grey')
dev.off()   



##### Slingshot using Clusters from PCA on PHATE embeddings knn 20 #####
phate20 <- readRDS(paste0(datadir,"phate_knn20.rds"))
sling20 <- slingshot(phate20$embedding, Idents(seurat), start.clus = 1)
sling20_0 <- slingshot(phate20$embedding, Idents(seurat), start.clus = 0)
saveRDS(sling, paste0(outdir,"slingshot_20knn.rds"))


png(paste0(outdir,"slingshot.png"))
plot(as.data.frame(phate$embedding), col = Idents(seurat), pch=20, asp = 1,
    xlab="PHATE 1", 
    ylab="PHATE 2",
    cex.main=1.25, 
    cex.lab=1, 
    cex.axis=1)
lines(sling, lwd=2, type = 'lineages', col = 'grey')
dev.off()   



png(paste0(outdir,"slingshot_curves.png"))
plot(as.data.frame(phate$embedding), col = Idents(seurat), pch=20, asp = 1,
    xlab="PHATE 1", 
    ylab="PHATE 2",
    cex.main=1.25, 
    cex.lab=1, 
    cex.axis=1)
lines(sling, lwd=2, col = 'grey')
dev.off()   