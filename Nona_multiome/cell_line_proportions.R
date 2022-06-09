##### Author: Drew Neavin
##### Date: 2 December, 2021
##### Reason: Look at the proportions of each line at each time of Nona's multi-ome experiment




##### Load in libraries #####
library(data.table)
library(Seurat)
library(tidyverse)
library("ggpomological")


##### Set up directories #####
datadir <- "/directflow/SCCGGroupShare/projects/himaro/imputing_snp/demultiplexing/demultiplex_Nona/processed_data_demultiplex2/log_dir/"
non_dir <- "/directflow/SCCGGroupShare/projects/nonfar/analysis/cardiac_multiome_directflow/demux_obj/"
outdir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Nona_multiome/"

dir.create(outdir, recursive = TRUE)


##### Get a list of the village pools #####
villages <- list.files(non_dir, pattern = "Village")
villages <- grep("DemuxALL", villages, value = TRUE)



##### Get the singlets from the file #####
village_id_list <- lapply(villages, function(x){
    print(x)
    # tmp <- fread(paste0(datadir,x,"/CombinedResults/Final_Assignments_demultiplexing_doublets_new_edit.txt"), sep = "\t")
    tmp <- readRDS(paste0(non_dir,x))
    dt <- data.table(tmp@meta.data)
    dt$Pool_ID <- gsub("_DemuxALL.rds", "",x))
    dt$Day <- as.numeric(as.character(gsub("Village_Day", "", dt$Pool_ID)))
    return(dt)
})

village_id <- do.call(rbind, village_id_list)

village_id$Pool_ID_updated <- gsub("Day7$", "Day5b", village_id$Pool_ID) %>%
                                gsub("Day5$", "Day7b", .) %>% 
                                    gsub("Day15$", "Day4b", .) %>%
                                        gsub("Day4$", "Day15b", .) %>%
                                            gsub("b", "", .)

village_id$Day_updated <- as.numeric(as.character(gsub("Village_Day", "", village_id$Pool_ID_updated)))

### Update columnes ###
village_id$Day <- factor(village_id$Day, levels = c(0,1,2,3,4,5,7,15))
village_id$Assignment <- gsub("^0_", "", village_id$Assignment) %>%
                                gsub("D-", "", .) %>%
                                    gsub("\\.\\d\\.", "", .) %>%
                                        gsub("N-", "", .) %>%
                                            gsub("-P36", "", .)  %>%
                                                gsub("-", "", .)

village_id$DropletType <- ifelse(village_id$DropletType != "singlet", "doublet", village_id$DropletType)

data.table(prop.table(table(village_id[,c("DropletType", "Day_updated")]), margin = 2))

village_summary <- data.table(prop.table(table(village_id[,c("Assignment", "Day_updated")]), margin = 2))

village_summary_singlets <- data.table(prop.table(table(village_id[Assignment != "unassigned" & Assignment != "doublet",c("Assignment", "Day_updated")]), margin = 2))

village_summary_singlets$Assignment <- factor(village_summary_singlets$Assignment, levels = rev(village_summary_singlets[Day_updated == 15]$Assignment[order(village_summary_singlets[Day_updated == 15]$N)]))


##### Make proportion plots (area plot) #####
p_stacked_area <- ggplot(village_summary_singlets, aes(x = as.numeric(as.character(Day_updated)), y = N, fill = factor(Assignment), group = Assignment)) +
    geom_area(alpha=0.6 , size=0.5, colour="black") +
    theme_classic() +
    scale_fill_manual(values = c("#f44336", "#e81f63", "#9c27b0", "#673ab7", "#3f51b5", "#2096f3","#2096f3", "#009688", "#4caf50", "#8bc34a", "#cddc39", "#ffeb3b", "#ffc108", "#ff9801", "#ff5723" ,"#795548", "#9e9e9e", "#607d8b")) +
    xlab("Days") +
    ylab("Proportion of Cells")
ggsave(p_stacked_area, filename = paste0(outdir,"stacked_area.png"), width = 7, height = 4)
ggsave(p_stacked_area, filename = paste0(outdir,"stacked_area.pdf"), width = 7, height = 4)




##### Make line plot of propotion over time #####
p_line <- ggplot(village_summary_singlets, aes(x = as.numeric(as.character(Day_updated)), y = N, color = Assignment)) +
    geom_point() +
    theme_classic() +
    geom_line() +
    scale_color_manual(values = c("#f44336", "#e81f63", "#9c27b0", "#673ab7", "#3f51b5", "#2096f3","#2096f3", "#009688", "#4caf50", "#8bc34a", "#cddc39", "#ffeb3b", "#ffc108", "#ff9801", "#ff5723" ,"#795548", "#9e9e9e", "#607d8b")) +
    xlab("Days") +
    ylab("Proportion of Cells")

ggsave(p_line, filename = paste0(outdir,"line_proportions.png"), width = 7, height = 4)
ggsave(p_line, filename = paste0(outdir,"line_proportions.pdf"), width = 7, height = 4)



##### Check QC metrics #####
### Load in Data ###
tenxdir <- "/directflow/SCCGGroupShare/projects/annsen/ATACseq/"

tenx_list <- lapply(villages, function(x){
    Read10X(paste0(tenxdir,x, "/outs/filtered_feature_bc_matrix"))
})


seurat_list <- lapply(tenx_list, function(x){
    CreateSeuratObject(counts = x)
})

seurat <- merge(seurat_list[[1]], y = seurat_list[2:length(seurat_list)], add.cell.ids = villages)

seurat$Pool <- gsub("_[GTAC]+-1", "", rownames(seurat@meta.data))
seurat$Day <- as.numeric(as.character(gsub("Village_Day", "", seurat$Pool)))

seurat$mt_percent <- PercentageFeatureSet(seurat, pattern = "^MT-")
seurat$rb_percent <- PercentageFeatureSet(seurat, pattern = "^RP[SL]")

### Generate some plots by pool to show the QC metrics ###
### N UMI ###
p_UMI_vnl <- VlnPlot(seurat, features = "nCount_RNA", group.by = "Day", split.by = "Pool", pt.size = 0) +
                    scale_fill_manual(values = c("#c03728", "#919c4c", "#f18721", "#f5c049", "#e68c7c", "#828585", "#c3c377", "#4f5157")) +
                    NoLegend()

ggsave(p_UMI_vnl, filename = paste0(outdir,"umi_vln.png"))


### N Genes ###
p_gene_vnl <- VlnPlot(seurat, features = "nFeature_RNA", group.by = "Day", split.by = "Day", pt.size = 0) +
                    scale_fill_manual(values = c("#c03728", "#919c4c", "#f18721", "#f5c049", "#e68c7c", "#828585", "#c3c377", "#4f5157")) +
                    NoLegend()

ggsave(p_gene_vnl, filename = paste0(outdir,"gene_vln.png"))


### Mt % ###
p_mt_vnl <- VlnPlot(seurat, features = "mt_percent", group.by = "Day", split.by = "Day", pt.size = 0) +
                    scale_fill_manual(values = c("#c03728", "#919c4c", "#f18721", "#f5c049", "#e68c7c", "#828585", "#c3c377", "#4f5157")) +
                    NoLegend()

ggsave(p_mt_vnl, filename = paste0(outdir,"mt_vln.png"))

### Rb % ###
p_rb_vnl <- VlnPlot(seurat, features = "rb_percent", group.by = "Day", split.by = "Day", pt.size = 0) +
                    scale_fill_manual(values = c("#c03728", "#919c4c", "#f18721", "#f5c049", "#e68c7c", "#828585", "#c3c377", "#4f5157")) +
                    NoLegend()

ggsave(p_rb_vnl, filename = paste0(outdir,"rb_vln.png"))

### Mt % vs N UMIs ###
plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "mt_percent") +
            scale_color_manual()


### N UMI vs N Genes ###
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

### Data was already pretty clean since just used intronic reads, probably don't need to filter further for the high quality cells