
##### Setting up Directories
dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
outdir <- paste0(dir,"output/Cluster_Pathway_Analysis/")
datadir <- paste0(dir,"output/All_data_integrated_remove_bad/")
datadir_DEGs <- paste0(datadir,"/DEGs/")


dir.create(outdir)


library(Seurat)
library(tidyverse)
library(ggplot2)
library(jcolors)
library(cowplot)
library(RColorBrewer)
library(readr)
library(purrr)
library(reticulate)
library(awtools)
library(Nebulosa)
library(plyr)
library(clusterProfiler)
library(enrichplot)
library(magrittr)


organism = "org.Hs.eg.db"
# BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)


wp2gene <- read.gmt("/directflow/SCCGGroupShare/projects/DrewNeavin/References/WikiPathways/wikipathways-20210510-gmt-Homo_sapiens.gmt")
wp2gene <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME

seurat <- readRDS(paste0(datadir,"seurat_integrated_all_times_clustered.rds"))


# degs <- FindAllMarkers(seurat, assay = "RNA", logfc.threshold = 0)

# saveRDS(degs, paste0(outdir,"DEG_clusters.rds"))
degs <- readRDS(paste0(outdir,"DEG_clusters.rds"))

degs_sig <- degs[which(degs$p_val_adj < 0.05),]



##### Add gene IDs for easy identification downstream #####
GeneConversion1 <- read_delim("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/DRENEA_1/outs/filtered_feature_bc_matrix/features.tsv.gz", col_names = F, delim = "\t")
GeneConversion2 <- read_delim("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/Village_B_1_week/outs/filtered_feature_bc_matrix/features.tsv.gz", col_names = F, delim = "\t")

GeneConversion <- unique(rbind(GeneConversion1, GeneConversion2))
GeneConversion <- GeneConversion[!duplicated(GeneConversion$X1),]
GeneConversion <- left_join(data.frame(rownames(seurat[["RNA"]])), GeneConversion, by= c("rownames.seurat...RNA...." = "X1"))
GeneConversion$X3 <- NULL
colnames(GeneConversion) <- c("gene", "Gene_ID")


degs_sig <- left_join(degs_sig, GeneConversion, by = c("gene"))

top_degs_sig <- degs_sig %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)




seurat_list <- list()
for (cluster in unique(Idents(seurat))){
	seurat_list[[cluster]] <- subset(seurat, idents = cluster)
}




cluster_deg_list <- list()
for (cluster in unique(degs$cluster)){
	cluster_deg_list[[cluster]] <- degs[which(degs$cluster == cluster),]
}


GeneList <- lapply(cluster_deg_list, function(x){
	tmp <- x$avg_log2FC
	names(tmp) <- x$gene
	tmp <- sort(tmp, decreasing = TRUE)
	return(tmp)
})


Genes <- lapply(GeneList, function(x){
	names(x)[abs(x) > 0.25]
})

lapply(Genes, length)

### Convert to entrz IDs ###
Genes_ENTREZ <- lapply(Genes, function(x){
	select(org.Hs.eg.db, keys=x, columns="ENTREZID", keytype="ENSEMBL")$ENTREZID
})


GeneList_ENTREZ <- lapply(GeneList, function(x){
	names(x) <- select(org.Hs.eg.db, keys=names(x), columns="ENTREZID", keytype="ENSEMBL")$ENTREZID
	return(x)
})




ewp <- lapply(Genes_ENTREZ, function(x){
	enricher(x, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
})

lapply(ewp, head)

ewp <- lapply(Genes_ENTREZ, function(x){
	GSEA(x, TERM2GENE = wpid2gene, TERM2NAME = wpid2name, verbose=TRUE)
})

# # we want the log2 fold change 
# original_gene_list <- lapply(cluster_deg_list, function(x){
# 	temp <- x$avg_log2FC

# 	# name the vector
# 	names(temp) <- x$gene

# 	# omit any NA values 
# 	temp <- na.omit(temp)

# 	# sort the list in decreasing order (required for clusterProfiler)
# 	temp = sort(temp, decreasing = TRUE)

# 	return(temp)
# })


# keytypes(org.Hs.eg.db)






# gseGO(geneList = original_gene_list[[6]],
#             #  ont ="ALL", 
#              keyType = "ENSEMBL", 
#              verbose = TRUE, 
#              OrgDb = organism)