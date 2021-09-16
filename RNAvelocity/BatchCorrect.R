library(anndata)
library(data.table)
library(tidyverse)
library(reticulate)
scv <- import("scvelo")


##### Directories #####
datadir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/scVelo/preprocess/seurat/"


##### Get Data #####
samples = fread('/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/scVelo/preprocess/velocyto_files.tsv', sep="\t")

lfile_list <- lapply(samples$Directory, function(file){
	scv$read(filename = file)
})

lapply(lfile_list, function(x){
	x$var_names_make_unique(join = ".")
})

adata_concat = concat(lfile_list, axis = 0L)

### Update rownames to match those in metadata files ###
adata_concat$obs_names <- gsub(":", "_", adata_concat$obs_names) %>% gsub("x","",.)


meta = fread(paste0(datadir, "metadata.csv"))



### Filter anndata object for droplet ids
adata_concat = adata_concat[rownames(adata_concat$obs) %in% meta$V1]


### Add metadata to object
rownames(meta) <- meta$V1
meta <- meta[match(adata_concat$obs_names, rownames(meta)),]
rownames(meta) <- meta$V1
all(adata_concat$obs_names == rownames(meta))


### Update necessary columns ###
meta$Location <- gsub("_.+", "", meta$Location_Time)
meta$Time <- gsub("Village Day 4","Village", meta$Time) %>% gsub("Thawed Village Day 7", "Village", .) %>% gsub("Thawed Village Day 0","Baseline", .)



adata_concat$obs$clusters = meta$integrated_snn_res.0.28
adata_concat$obs$individual = meta$Final_Assignment
adata_concat$obs$cell_cycle = meta$phases
adata_concat$obs$site_rep = meta$Site_rep
adata_concat$obs$Location = meta$Location
adata_concat$obs$Time = meta$Time


##### Get just the Brisbane FSA0006 cells #####
adata_sub <- adata_concat[adata_concat$obs$Location == "Brisbane"]

M <- matrix(adata_sub$layers$spliced) + (adata_sub$layers$unspliced)
R <- matrix(adata_sub$layers$spliced)/M
