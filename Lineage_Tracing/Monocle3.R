library(tidyverse)
library(monocle3)
library(Seurat)


##### Set up directories #####
dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
datadir <- paste0(dir,"output/LineageTracing/Phate/")
outdir <- paste0(dir,"output/LineageTracing/Monocle3/")
dir.create(outdir, recursive = TRUE)



##### Read in Phate #####
phate20 <- readRDS(paste0(datadir,"phate_knn20.rds"))
seurat <- readRDS(paste0(dir,"output/Distribution_tests/seurat_integrated_all_times_clustered_1pct_expressing.rds"))


##### Reformat data for monocle 3 #####
cds <- new_cell_data_set(seurat@,
                          cell_metadata = cell_metadata,
                          gene_metadata = gene_annotation)



##### Add Clusters #####
cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition")




##### Learn Trajectory #####
cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "cell.type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)