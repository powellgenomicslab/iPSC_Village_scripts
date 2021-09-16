library(tidyverse)
library(phateR)
library(ggplot2)
library(Rmagic)
library(Seurat)
library(viridis)


dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
datadir <- paste0(dir,"output/All_data_integrated_remove_bad/")
outdir <- paste0(dir,"output/LineageTracing/Phate/")
dir.create(outdir, recursive = TRUE)


### Set up colors ###
line_colors = c(FSA0006 = "#F79E29", MBE1006 = "#9B2C99", TOB0421 = "#35369C")
village_colors = c(Baseline = "#7d57a0", Village = "#853786")
location_colors = c(Brisbane = "#C9D8EA", Melbourne = "#5D59AB", Sydney = "#179085", Sydney_Cryopreserved = "#A7C9A9")



seurat <- readRDS(paste0(datadir,"seurat_integrated_all_times_clustered.rds"))
seurat@meta.data$Location <- gsub("_.+", "", seurat@meta.data$Location)
seurat@meta.data$Location <- ifelse(grepl("Thawed", seurat@meta.data$Location_Time), "Sydney_Cryopreserved", seurat@meta.data$Location)


### Double check that library size filtering was OK
pLibSizeDist <- ggplot() +
  geom_histogram(aes(x=rowSums(seurat[["RNA"]]@counts)), bins=50) +
  geom_vline(xintercept = 1000, color='red') +
  theme_classic()
ggsave(pLibSizeDist, filename = paste0(outdir, "lib_size_dist.png"))

pLibSizeDist <- ggplot() +
  geom_histogram(aes(x=rowSums(seurat[["SCT"]]@counts)), bins=50) +
  geom_vline(xintercept = 1000, color='red') +
  theme_classic()
ggsave(pLibSizeDist, filename = paste0(outdir, "lib_size_dist_SCT.png"))

min(colSums(seurat[["RNA"]]@counts))

### Run PCA
bmmsc_PCA <- as.data.frame(prcomp(t(seurat[["integrated"]]@scale.data))$x)


pPCA_POU5F1 <- ggplot(bmmsc_PCA) +
  geom_point(aes(PC1, PC2, color=seurat[["SCT"]]@scale["ENSG00000204531",])) +
  labs(color="POU5F1") +
  scale_color_viridis(option="B") +
  theme_classic()
ggsave(pPCA_POU5F1, filename = paste0(outdir, "PCA_POU5F1.png"))


if(!file.exists(paste0(outdir,"phate_knn5.rds"))){
	bmmsc_PHATE <- phate(t(seurat[["integrated"]]@scale.data))
	saveRDS(bmmsc_PHATE, paste0(outdir,"phate_knn5.rds"))
} else {
	bmmsc_PHATE <- readRDS(paste0(outdir,"phate_knn5.rds"))
}


pPHATE <- ggplot(bmmsc_PHATE) +
  geom_point(aes(PHATE1, PHATE2, color=seurat[["SCT"]]@scale.data["ENSG00000204531",]), size = 0.5) +
  labs(color="POU5F1") +
  scale_color_viridis(option="B") +
  theme_classic()
ggsave(pPHATE, filename = paste0(outdir, "PHATE_POU5F1.png"))


pPHATE_village <- ggplot(bmmsc_PHATE) +
  geom_point(aes(PHATE1, PHATE2, color=ifelse(seurat@meta.data$Time == "Baseline", "Baseline", "Village")), size = 0.5) +
  labs(color="Village") +
  scale_color_manual(values = village_colors) +
  theme_classic()
ggsave(pPHATE_village, filename = paste0(outdir, "PHATE_village.png"))


pPHATE_line <- ggplot(bmmsc_PHATE) +
  geom_point(aes(PHATE1, PHATE2, color = seurat@meta.data$Final_Assignment), size = 0.5) +
  labs(color="Village") +
  scale_color_manual(values = line_colors) +
  theme_classic()
ggsave(pPHATE_line, filename = paste0(outdir, "PHATE_cell_liune.png"))



pPHATE_clusters <- ggplot(bmmsc_PHATE) +
  geom_point(aes(PHATE1, PHATE2, color=Idents(seurat)), size = 0.5) +
  labs(color="Clusters") +
  theme_classic()
ggsave(pPHATE_clusters, filename = paste0(outdir, "PHATE_clusters.png"))


pPHATE_clusters_facet <- ggplot(bmmsc_PHATE) +
  geom_point(aes(PHATE1, PHATE2, color=Idents(seurat)), size = 0.5) +
  labs(color="Clusters") +
  theme_classic() +
  facet_wrap(vars(Idents(seurat)))
ggsave(pPHATE_clusters_facet, filename = paste0(outdir, "PHATE_clusters_facet.png"))




pPHATE_location <- ggplot(bmmsc_PHATE) +
  geom_point(aes(PHATE1, PHATE2, color = seurat@meta.data$Location), size = 0.5) +
  labs(color="Location") +
  scale_color_manual(values = location_colors) +
  theme_classic()
ggsave(pPHATE_location, filename = paste0(outdir, "PHATE_location.png"))


pPHATE_location <- ggplot(bmmsc_PHATE) +
  geom_point(aes(PHATE1, PHATE2, color = seurat@meta.data$Location), size = 0.5) +
  labs(color="Location") +
  scale_color_manual(values = location_colors) +
  theme_classic()
ggsave(pPHATE_location, filename = paste0(outdir, "PHATE_location.png"))



clusters <- cluster_phate(bmmsc_PHATE)


pPHATE_clusters_phate <- ggplot(bmmsc_PHATE) +
  geom_point(aes(PHATE1, PHATE2, color=factor(clusters)), size = 0.5) +
  labs(color="POU5F1")  +
  theme_classic()
ggsave(pPHATE_clusters_phate, filename = paste0(outdir, "PHATE_clusters_phate.png"))





### Try phate with 10 knn instead of 5 ###

if(!file.exists(paste0(outdir,"phate_knn10.rds"))){
	bmmsc_PHATE2 <- phate(t(seurat[["integrated"]]@scale.data), knn= 10, init=bmmsc_PHATE)
	saveRDS(bmmsc_PHATE2, paste0(outdir,"phate_knn10.rds"))
} else {
	bmmsc_PHATE2 <- readRDS(paste0(outdir,"phate_knn10.rds"))
}



pPHATE2 <- ggplot(bmmsc_PHATE2) +
  geom_point(aes(PHATE1, PHATE2, color=seurat[["SCT"]]@data["ENSG00000204531",]), size = 0.5) +
  labs(color="POU5F1") +
  scale_color_viridis(option="B") +
  theme_classic()
ggsave(pPHATE2, filename = paste0(outdir, "PHATE_POU5F1_knn10.png"))


pPHATE_village2 <- ggplot(bmmsc_PHATE2) +
  geom_point(aes(PHATE1, PHATE2, color=ifelse(seurat@meta.data$Time == "Baseline", "Baseline", "Village")), size = 0.5) +
  labs(color="Village") +
  scale_color_manual(values = village_colors) +
  theme_classic()
ggsave(pPHATE_village2, filename = paste0(outdir, "PHATE_village_knn10.png"))


pPHATE_clusters2 <- ggplot(bmmsc_PHATE2) +
  geom_point(aes(PHATE1, PHATE2, color=Idents(seurat)), size = 0.5) +
  labs(color="Clusters") +
  theme_classic()
ggsave(pPHATE_clusters2, filename = paste0(outdir, "PHATE_clusters_knn10.png"))



pPHATE_location2 <- ggplot(bmmsc_PHATE2) +
  geom_point(aes(PHATE1, PHATE2, color = seurat@meta.data$Location), size = 0.5) +
  labs(color="Location") +
  scale_color_manual(values = location_colors) +
  theme_classic()
ggsave(pPHATE_location2, filename = paste0(outdir, "PHATE_location_knn10.png"))


pPHATE2_ID1 <- ggplot(bmmsc_PHATE2) +
  geom_point(aes(PHATE1, PHATE2, color=seurat[["SCT"]]@data["ENSG00000125968",]), size = 0.5) +
  labs(color="POU5F1") +
  scale_color_viridis(option="B")
ggsave(pPHATE2_ID1, filename = paste0(outdir, "PHATE_ID1_knn10.png"))


##### Cluster #####
clusters2 <- cluster_phate(bmmsc_PHATE2)

pPHATE2_clusters <- ggplot(bmmsc_PHATE2) +
  geom_point(aes(PHATE1, PHATE2, color=factor(clusters)), size = 0.5) +
  labs(color="POU5F1")  +
  theme_classic()
ggsave(pPHATE2_clusters, filename = paste0(outdir, "PHATE_clusters_phate_knn10.png"))


##### Try with knn = 20 #####
if(!file.exists(paste0(outdir,"phate_knn20.rds"))){
	bmmsc_PHATE3 <- phate(t(seurat[["integrated"]]@scale.data), knn= 20, init=bmmsc_PHATE)
	saveRDS(bmmsc_PHATE3, paste0(outdir,"phate_knn20.rds"))
} else {
	bmmsc_PHATE3 <- readRDS(paste0(outdir,"phate_knn20.rds"))
}


pPHATE3 <- ggplot(bmmsc_PHATE3) +
  geom_point(aes(PHATE1, PHATE2, color=seurat[["SCT"]]@data["ENSG00000204531",]), size = 0.5) +
  labs(color="POU5F1") +
  scale_color_viridis(option="B") +
  theme_classic()
ggsave(pPHATE3, filename = paste0(outdir, "PHATE_POU5F1_knn20.png"))


pPHATE_village3 <- ggplot(bmmsc_PHATE3) +
  geom_point(aes(PHATE1, PHATE2, color=ifelse(seurat@meta.data$Time == "Baseline", "Baseline", "Village")), size = 0.5) +
  labs(color="Village") +
  scale_color_manual(values = village_colors) +
  theme_classic()
ggsave(pPHATE_village3, filename = paste0(outdir, "PHATE_village_knn20.png"))


pPHATE_clusters3 <- ggplot(bmmsc_PHATE3) +
  geom_point(aes(PHATE1, PHATE2, color=Idents(seurat)), size = 0.5) +
  labs(color="Clusters") +
  theme_classic()
ggsave(pPHATE_clusters3, filename = paste0(outdir, "PHATE_clusters_knn20.png"))



pPHATE_location3 <- ggplot(bmmsc_PHATE3) +
  geom_point(aes(PHATE1, PHATE2, color = seurat@meta.data$Location), size = 0.5) +
  labs(color="Location") +
  scale_color_manual(values = location_colors) +
  theme_classic()
ggsave(pPHATE_location3, filename = paste0(outdir, "PHATE_location_knn20.png"))


pPHATE3_ID1 <- ggplot(bmmsc_PHATE3) +
  geom_point(aes(PHATE1, PHATE2, color=seurat[["SCT"]]@data["ENSG00000125968",]), size = 0.5) +
  labs(color="ID1") +
  scale_color_viridis(option="B")
ggsave(pPHATE3_ID1, filename = paste0(outdir, "PHATE_ID1_knn20.png"))


##### Cluster #####
clusters3 <- cluster_phate(bmmsc_PHATE3)

pPHATE3_clusters <- ggplot(bmmsc_PHATE3) +
  geom_point(aes(PHATE1, PHATE2, color=factor(clusters)), size = 0.5) +
  labs(color="POU5F1")  +
  theme_classic()
ggsave(pPHATE3_clusters, filename = paste0(outdir, "PHATE_clusters_phate_knn20.png"))
