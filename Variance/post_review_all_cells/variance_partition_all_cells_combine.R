##### Reason: combine the results for the variance explained by different factors for each gene
##### Author: Drew Neavin
##### Date: 14 March, 2022


##### Load in libraries #####
library(data.table)
library(tidyverse)
library(ggridges)
library(raincloudplots)
library(ggdist)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOSemSim)
library(RColorBrewer)



##### Set up directories #####
dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
icc_dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/variance_partitioning_all_cells/gene_separated/icc/"
icc_interaction_dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/variance_partitioning_all_cells/gene_separated/icc_interaction/"
outdir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/variance_partitioning_all_cells/combined/"

dir.create(outdir, recursive = TRUE)


vars <- c("Line", "Village", "Cryopreserved", "Site",  "Replicate","Line:Village", "Line:Cryopreserved", "Line:Site", "Village:Cryopreserved","Village:Site",  "Replicate:Village", "Replicate:Line","Replicate:Cryopreserved",  "Replicate:Site", "Residual")
var_colors <- c("#4734a9", rev(c("#115cc7", "#a4c3c8", "#499090", "#405940", "#685a54", "#f7d312", "#f8bf33", "#e4910e", "#f65d19", "#931519", "#f2c1ce", "#e17aab", "#a186aa")), "gray90")
names(var_colors) <- vars




##### Get list of icc files #####
icc_files <- list.files(icc_dir)



##### Read in icc results #####
icc_results_list <- lapply(icc_files, function(x){
    readRDS(paste0(icc_dir,x))
})
names(icc_results_list) <- icc_files



##### Get list of icc interaction files #####
icc_interaction_files <- list.files(icc_interaction_dir, pattern = "_icc.rds")



##### Read in icc results #####
icc_interaction_results_list <- lapply(icc_interaction_files, function(x){
    readRDS(paste0(icc_interaction_dir,x))
})
names(icc_interaction_results_list) <- icc_interaction_files



##### Merge icc results into a single data.table #####
icc_dt <- do.call(rbind, icc_results_list)

icc_dt$percent_round <- round(icc_dt$percent)

icc_dt$grp <- factor(icc_dt$grp, levels= rev(vars))

group_size  <- data.table(table(icc_dt$grp))
colnames(group_size) <- c("grp", "size")
group_size$grp_size <- paste0(group_size$grp, "\nN = ", group_size$size)

icc_dt <- group_size[icc_dt, on = "grp"]
icc_dt$grp_size <- factor(icc_dt$grp_size, levels = unique(group_size$grp_size))


##### Merge icc_interaction results into a single data.table #####
icc_interaction_dt <- do.call(rbind, icc_interaction_results_list)

icc_interaction_dt$percent_round <- round(icc_interaction_dt$percent)

icc_interaction_dt$grp <- factor(icc_interaction_dt$grp, levels= rev(vars))

group_size  <- data.table(table(icc_interaction_dt$grp))
colnames(group_size) <- c("grp", "size")
group_size$grp_size <- paste0(group_size$grp, "\nN = ", group_size$size)

icc_interaction_dt <- group_size[icc_interaction_dt, on = "grp"]
icc_interaction_dt$grp_size <- factor(icc_interaction_dt$grp_size, levels = unique(group_size$grp_size))

### *** Need to add individual effects without interaction in to interaction dt *** ###
icc_interaction_plus_dt <- rbind(icc_interaction_dt, icc_dt[!(gene %in% icc_interaction_dt$gene)])


group_size  <- data.table(table(icc_interaction_plus_dt$grp))
colnames(group_size) <- c("grp", "size")
group_size$grp_size <- paste0(group_size$grp, "\nN = ", group_size$size)

icc_interaction_plus_dt <- group_size[icc_interaction_plus_dt, on = "grp"]
icc_interaction_plus_dt$grp_size <- factor(icc_interaction_plus_dt$grp_size, levels = unique(group_size$grp_size))



##### Check difference in percent explained with and without interactions #####
icc_interaction_dt_joined <- icc_dt[icc_interaction_dt, on = c("grp", "gene")]

icc_interaction_dt_joined$difference <- icc_interaction_dt_joined$percent - icc_interaction_dt_joined$i.percent


pRaincloud_dif <- ggplot(icc_interaction_dt_joined, aes(x = difference, y = factor(grp_size, levels = rev(levels(grp_size))), fill = factor(grp, levels = rev(vars)))) + 
                geom_density_ridges(stat = "binline", bins = 90, scale = 0.7, draw_baseline = FALSE, aes(height =..ndensity..), alpha = 0.75) +
                geom_boxplot(size = 0.5,width = .15, outlier.size = 0.25, position = position_nudge(y=-0.12), alpha = 0.75) +
                # coord_cartesian(xlim = c(1.2, NA), clip = "off") +
                theme_classic() +
                theme(axis.title.y=element_blank()) +
                xlab("Percent Variance Explained") +
                scale_y_discrete(expand = c(0.03, 0)) +
                scale_fill_manual(values = var_colors)

ggsave(pRaincloud_dif, filename = paste0(outdir, "variance_explained_interaction_difference_raincloud.png"), height = 8, width = 7)
ggsave(pRaincloud_dif, filename = paste0(outdir, "variance_explained_interaction_difference_raincloud.pdf"), height = 8, width = 7)



##### Make a figure of stacked variance explained #####
### Order based on line variance explained ###
genes_list <- list()

for (group in vars){
    genes_list[[group]] <- icc_dt[grp == group][rev(order(percent_round))]$gene
}

genes <- unique(unlist(genes_list))

icc_dt$gene <- factor(icc_dt$gene, levels = genes)

icc_dt$grp <- factor(icc_dt$grp, levels= rev(vars))



pRaincloud <- ggplot(icc_dt, aes(x = percent, y = factor(grp_size, levels = rev(levels(grp_size))), fill = factor(grp, levels = rev(vars)))) + 
                geom_density_ridges(stat = "binline", bins = 90, scale = 0.7, draw_baseline = FALSE, aes(height =..ndensity..), alpha = 0.75) +
                geom_boxplot(size = 0.5,width = .15, outlier.size = 0.25, position = position_nudge(y=-0.12), alpha = 0.75) +
                coord_cartesian(xlim = c(1.2, NA), clip = "off") +
                theme_classic() +
                theme(axis.title.y=element_blank()) +
                xlab("Percent Variance Explained") +
                scale_y_discrete(expand = c(0.03, 0)) +
                scale_fill_manual(values = var_colors)

ggsave(pRaincloud, filename = paste0(outdir, "variance_explained_raincloud.png"), height = 8, width = 7)
ggsave(pRaincloud, filename = paste0(outdir, "variance_explained_raincloud.pdf"), height = 8, width = 7)



icc_interaction_plus_dt$grp_size <- factor(icc_interaction_plus_dt$grp_size, levels = c(unique(grep("^Line\nN = ", icc_interaction_plus_dt$grp_size, value = TRUE)), 
                                                                                       unique(grep("^Village\nN = ",  icc_interaction_plus_dt$grp_size, value = TRUE)), 
                                                                                       unique(grep("^Site\nN = ",  icc_interaction_plus_dt$grp_size, value = TRUE)), 
                                                                                       unique(grep("^Cryopreserved\nN =",  icc_interaction_plus_dt$grp_size, value = TRUE)), 
                                                                                       unique(grep("^Replicate\nN =",  icc_interaction_plus_dt$grp_size, value = TRUE)), 
                                                                                       unique(grep("Line:Village\nN =",  icc_interaction_plus_dt$grp_size, value = TRUE)), 
                                                                                       unique(grep("Line:Site\nN =",  icc_interaction_plus_dt$grp_size, value = TRUE)), 
                                                                                       unique(grep("Line:Cryopreserved\nN =",  icc_interaction_plus_dt$grp_size, value = TRUE)), 
                                                                                       unique(grep("Village:Site\nN =",  icc_interaction_plus_dt$grp_size, value = TRUE)), 
                                                                                       unique(grep("Village:Cryopreserved\nN =",  icc_interaction_plus_dt$grp_size, value = TRUE)), 
                                                                                       unique(grep("Replicate:Village\nN ",  icc_interaction_plus_dt$grp_size, value = TRUE)), 
                                                                                       unique(grep("Replicate:Line\nN ",  icc_interaction_plus_dt$grp_size, value = TRUE)), 
                                                                                       unique(grep("Replicate:Site\nN ",  icc_interaction_plus_dt$grp_size, value = TRUE)), 
                                                                                       unique(grep("Replicate:Cryopreserved\nN ",  icc_interaction_plus_dt$grp_size, value = TRUE)), 
                                                                                       unique(grep("Residual\nN = ", icc_interaction_plus_dt$grp_size, value = TRUE))))


pRaincloud_interaction <- ggplot(icc_interaction_plus_dt, aes(x = percent, y = grp_size, fill = factor(grp, levels = rev(vars)))) + 
                geom_density_ridges(stat = "binline", bins = 100, scale = 0.7, draw_baseline = FALSE, aes(height =..ndensity..), alpha = 0.75) +
                geom_boxplot(size = 0.5,width = .15, outlier.size = 0.25, position = position_nudge(y=-0.12), alpha = 0.75) +
                coord_cartesian(xlim = c(1.2, NA), clip = "off") +
                theme_classic() +
                theme(axis.title.y=element_blank()) +
                xlab("Percent Variance Explained") +
                scale_y_discrete(expand = c(0.03, 0)) +
                scale_fill_manual(values = var_colors) +
                geom_vline(xintercept = 1, linetype = "dashed", color = "red")

ggsave(pRaincloud_interaction, filename = paste0(outdir, "variance_explained_raincloud_interaction.png"), height = 8, width = 7)
ggsave(pRaincloud_interaction, filename = paste0(outdir, "variance_explained_raincloud_interaction.pdf"), height = 8, width = 7)



icc_interaction_sig_list <- list()

for (ensg in unique(icc_interaction_plus_dt$gene)){
    icc_interaction_sig_list[[ensg]] <- icc_interaction_plus_dt[gene == ensg][P < 0.05/(nrow(icc_interaction_plus_dt[gene == ensg])-1)]
    icc_interaction_sig_list[[ensg]] <- rbind(icc_interaction_sig_list[[ensg]], icc_interaction_plus_dt[gene == ensg & grp == "Residual"])
}

icc_interaction_sig_dt <- do.call(rbind, icc_interaction_sig_list) 



group_size  <- data.table(table(icc_interaction_sig_dt$grp))
colnames(group_size) <- c("grp", "size")
group_size$grp_size <- paste0(group_size$grp, "\nN = ", group_size$size)

icc_interaction_sig_dt <- group_size[icc_interaction_sig_dt, on = "grp"]
icc_interaction_sig_dt$grp_size <- factor(icc_interaction_sig_dt$grp_size, levels = c(unique(grep("^Line\nN = ", icc_interaction_sig_dt$grp_size, value = TRUE)), 
                                                                                       unique(grep("^Village\nN = ",  icc_interaction_sig_dt$grp_size, value = TRUE)), 
                                                                                       unique(grep("^Site\nN = ",  icc_interaction_sig_dt$grp_size, value = TRUE)), 
                                                                                       unique(grep("^Cryopreserved\nN =",  icc_interaction_sig_dt$grp_size, value = TRUE)), 
                                                                                       unique(grep("^Replicate\nN =",  icc_interaction_sig_dt$grp_size, value = TRUE)), 
                                                                                       unique(grep("Line:Village\nN =",  icc_interaction_sig_dt$grp_size, value = TRUE)), 
                                                                                       unique(grep("Line:Site\nN =",  icc_interaction_sig_dt$grp_size, value = TRUE)), 
                                                                                       unique(grep("Line:Cryopreserved\nN =",  icc_interaction_sig_dt$grp_size, value = TRUE)), 
                                                                                       unique(grep("Village:Site\nN =",  icc_interaction_sig_dt$grp_size, value = TRUE)), 
                                                                                       unique(grep("Village:Cryopreserved\nN =",  icc_interaction_sig_dt$grp_size, value = TRUE)), 
                                                                                       unique(grep("Replicate:Village\nN ",  icc_interaction_sig_dt$grp_size, value = TRUE)), 
                                                                                       unique(grep("Replicate:Line\nN ",  icc_interaction_sig_dt$grp_size, value = TRUE)), 
                                                                                       unique(grep("Replicate:Site\nN ",  icc_interaction_sig_dt$grp_size, value = TRUE)), 
                                                                                       unique(grep("Replicate:Cryopreserved\nN ",  icc_interaction_sig_dt$grp_size, value = TRUE)), 
                                                                                       unique(grep("Residual\nN = ", icc_interaction_sig_dt$grp_size, value = TRUE))))




pRaincloud_interaction_sig <- ggplot(icc_interaction_sig_dt, aes(x = percent, y = grp_size, fill = factor(grp, levels = rev(vars)))) + 
                geom_density_ridges(stat = "binline", bins = 100, scale = 0.7, draw_baseline = FALSE, aes(height =..ndensity..), alpha = 0.75) +
                geom_boxplot(size = 0.5,width = .15, outlier.size = 0.25, position = position_nudge(y=-0.12), alpha = 0.75) +
                coord_cartesian(xlim = c(1.2, NA), clip = "off") +
                theme_classic() +
                theme(axis.title.y=element_blank()) +
                xlab("Percent Variance Explained") +
                scale_y_discrete(expand = c(0.03, 0)) +
                scale_fill_manual(values = var_colors) +
                geom_vline(xintercept = 1, linetype = "dashed", color = "firebrick3")

ggsave(pRaincloud_interaction_sig, filename = paste0(outdir, "variance_explained_raincloud_interaction_significant.png"), height = 8, width = 7)
ggsave(pRaincloud_interaction_sig, filename = paste0(outdir, "variance_explained_raincloud_interaction_significant.pdf"), height = 8, width = 7)


total <- icc_interaction_sig_dt[,.(count = .N), by = .(grp)]
total_less1pct <- icc_interaction_sig_dt[percent <= 1][,.(count_less_1pct = .N), by = .(grp)]
total_less5pct <- icc_interaction_sig_dt[percent <= 5][,.(count_less_5pct = .N), by = .(grp)]
total_less10pct <- icc_interaction_sig_dt[percent <= 10][,.(count_less_10pct = .N), by = .(grp)]
summary <- total[total_less1pct, on = "grp"]
summary <- summary[total_less5pct, on = "grp"]
summary <- summary[total_less10pct, on = "grp"]
summary$count_greater_1pct <- summary$count - summary$count_less_1pct
summary$count_greater_5pct <- summary$count - summary$count_less_5pct
summary$count_greater_10pct <- summary$count - summary$count_less_10pct
summary$percent_1pct <- (summary$count_less_1pct/summary$count)*100
summary$percent_5pct <- (summary$count_less_5pct/summary$count)*100
summary$percent_10pct <- (summary$count_less_10pct/summary$count)*100


pRaincloud_interaction_sig_1pct <- ggplot(icc_interaction_sig_dt[percent >= 1], aes(x = percent, y = grp_size, fill = factor(grp, levels = rev(vars)))) + 
                geom_density_ridges(stat = "binline", bins = 90, scale = 0.7, draw_baseline = FALSE, aes(height =..ndensity..), alpha = 0.75) +
                geom_boxplot(size = 0.5,width = .15, outlier.size = 0.25, position = position_nudge(y=-0.12), alpha = 0.75) +
                coord_cartesian(xlim = c(1.2, NA), clip = "off") +
                theme_classic() +
                theme(axis.title.y=element_blank()) +
                xlab("Percent Variance Explained") +
                scale_y_discrete(expand = c(0.03, 0)) +
                scale_fill_manual(values = var_colors) +
                geom_vline(xintercept = 1, linetype = "dashed", color = "firebrick3")

ggsave(pRaincloud_interaction_sig_1pct, filename = paste0(outdir, "variance_explained_raincloud_interaction_significant_1pct.png"), height = 8, width = 7)
ggsave(pRaincloud_interaction_sig_1pct, filename = paste0(outdir, "variance_explained_raincloud_interaction_significant_1pct.pdf"), height = 8, width = 7)



icc_interaction_sig_dt[gene == "ENSG00000106153"]
icc_dt[gene == "ENSG00000106153"]








##### Add gene IDs for easy identification downstream #####
GeneConversion1 <- read_delim("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/DRENEA_1/outs/filtered_feature_bc_matrix/features.tsv.gz", col_names = F, delim = "\t")
GeneConversion2 <- read_delim("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/Village_B_1_week/outs/filtered_feature_bc_matrix/features.tsv.gz", col_names = F, delim = "\t")

GeneConversion <- unique(rbind(GeneConversion1, GeneConversion2))
GeneConversion <- GeneConversion[!duplicated(GeneConversion$X1),]
GeneConversion$X3 <- NULL
colnames(GeneConversion) <- c("gene", "Gene_ID")

GeneConversion <- data.table(GeneConversion)


### Add the gene IDs to the icc_dt ###

icc_interaction_plus_dt <- GeneConversion[icc_interaction_plus_dt, on = "gene"]

icc_interaction_plus_dt[grp == "Cryopreserved"][rev(order(percent))]$Gene_ID
icc_interaction_plus_dt[grp == "Line"][rev(order(percent))]$Gene_ID
icc_interaction_plus_dt[grp == "Village"][rev(order(percent))]$Gene_ID

fwrite(icc_interaction_plus_dt, paste0(outdir, "sig_results.tsv.gz"), sep = "\t", compress = "gzip")

## Highlight 
##          X chromosome genes - wouldn't expect these to be Line-biased because expressed by both males and females
##          Y chromosome genes - should be line-biased because expressed by only males and have some male(s) and some female(s)
##          mt genes
##          ribosomal genes
##          look at gsea and kegg pathways for each

## Read in gtf used as reference and pull just X, Y or MT chromosome genes from it, use ribosomal file for rb genes
gtf <- fread("/directflow/GWCCGPipeline/projects/reference/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf", sep = "\t", autostart = 6, header = FALSE)

gtf_genes <- gtf[!(grep("transcript_id", V9))]

gtf_genes$V9 <- gsub("gene_id \"", "",gtf_genes$V9 ) %>%
            gsub("\"; gene_version \"", ";", .) %>%
            gsub("\"; gene_name \"", ";", .) %>%
            gsub("\"; gene_source \"", ";", .) %>%
            gsub("\"; gene_biotype \"", ";", .) %>%
            gsub("\"", "", .)

gtf_genes[, c("gene_id", "gene_version", "gene_name", "gene_source", "gene_biotype") := data.table(str_split_fixed(V9,";", 5))]


X_chromosome_genes <- gtf_genes[V1 == "X"]
X_genelist <- X_chromosome_genes$gene_id[X_chromosome_genes$gene_id %in% genes]
Y_chromosome_genes <- gtf_genes[V1 == "Y"]
Y_genelist <- Y_chromosome_genes$gene_id[Y_chromosome_genes$gene_id %in% genes]
MT_chromosome_genes <- gtf_genes[V1 == "MT"]
MT_genelist <- MT_chromosome_genes$gene_id[MT_chromosome_genes$gene_id %in% genes]
RbGeneList <- read.delim(file = "/directflow/SCCGGroupShare/projects/DrewNeavin/References/RibosomalGeneList_GeneID_ENSG.txt")
Rb_genelist <- RbGeneList$ENSG[RbGeneList$ENSG %in% genes]


### Make stacked bar plots of the variance explained by different factors for these gene groups




## Figure of x chromosome genes ##
icc_x <- icc_interaction_plus_dt[data.table(gene = icc_interaction_plus_dt[grp == "Residual"][gene %in% X_genelist][order(percent_round)]$gene), on = "gene"]
icc_x$grp <- factor(icc_x$grp, levels = rev(vars))
icc_x$gene <- factor(icc_x$gene, levels = unique(icc_x$gene))

bar_proportions_x <- ggplot(icc_x, aes(x = gene, y = percent, fill = grp)) +
    geom_bar(position="stack", stat="identity", alpha = 0.75) +
    theme_classic() +
    theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = var_colors) +
    ggtitle("Variance Explained of\nX Chromosome Genes")


ggsave(bar_proportions_x, filename = paste0(outdir, "variance_explained_bar_x_genes.png"), width = 20)


## Figure of y chromosome genes ##
icc_y <- icc_interaction_plus_dt[data.table(gene = icc_interaction_plus_dt[grp == "Residual"][gene %in% Y_genelist][order(percent_round)]$gene), on = "gene"]
icc_y$grp <- factor(icc_y$grp, levels = rev(vars))
icc_y$gene <- factor(icc_y$gene, levels = unique(icc_y$gene))
icc_y$Gene_ID <- factor(icc_y$Gene_ID, levels = unique(icc_y$Gene_ID))


bar_proportions_y <- ggplot(icc_y, aes(x = Gene_ID, y = percent, fill = grp)) +
    geom_bar(position="stack", stat="identity", alpha = 0.75) +
    theme_classic() +
    theme(axis.title.x=element_blank(),
        # axis.text.x=element_blank(),
        # axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_manual(values = var_colors) +
    scale_y_continuous(expand = c(0, 0)) +
    ggtitle("Variance Explained of\nY Chromosome Genes") +
    ylab("Percent")


ggsave(bar_proportions_y, filename = paste0(outdir, "variance_explained_bar_y_genes.png"), width = 4.5, height = 4.5)
ggsave(bar_proportions_y, filename = paste0(outdir, "variance_explained_bar_y_genes.pdf"), width = 4.5, height = 4.5)



## Figure of mt chromosome genes ##
icc_mt <- icc_interaction_plus_dt[data.table(gene = icc_interaction_plus_dt[grp == "Residual"][gene %in% MT_genelist][order(percent_round)]$gene), on = "gene"]
icc_mt$grp <- factor(icc_mt$grp, levels = rev(vars))
icc_mt$gene <- factor(icc_mt$gene, levels = unique(icc_mt$gene))
icc_mt$Gene_ID <- factor(icc_mt$Gene_ID, levels = unique(icc_mt$Gene_ID))


bar_proportions_mt <- ggplot(icc_mt, aes(x = Gene_ID, y = percent, fill = grp)) +
    geom_bar(position="stack", stat="identity", alpha = 0.75) +
    theme_classic() +
    theme(axis.title.x=element_blank(),
        # axis.text.x=element_blank(),
        # axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_manual(values = var_colors) +
    scale_y_continuous(expand = c(0, 0)) +
    ggtitle("Variance Explained of\nMitochondrial Genes") +
    ylab("Percent")

ggsave(bar_proportions_mt, filename = paste0(outdir, "variance_explained_bar_mt_genes.png"), width = 4.5, height = 4.5)
ggsave(bar_proportions_mt, filename = paste0(outdir, "variance_explained_bar_mt_genes.pdf"), width = 4.5, height = 4.5)



## Figure of mt chromosome genes ##
icc_rb <- icc_interaction_plus_dt[data.table(gene = icc_interaction_plus_dt[grp == "Residual"][gene %in% Rb_genelist][order(percent_round)]$gene), on = "gene"]
icc_rb$grp <- factor(icc_rb$grp, levels = rev(vars))
icc_rb$gene <- factor(icc_rb$gene, levels = unique(icc_rb$gene))


bar_proportions_rb <- ggplot(icc_rb, aes(x = gene, y = percent, fill = grp)) +
    geom_bar(position="stack", stat="identity", alpha = 0.75) +
    theme_classic() +
    theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = var_colors) +
    scale_y_continuous(expand = c(0, 0)) +
    ggtitle("Variance Explained\nof Ribosomal Genes")


ggsave(bar_proportions_rb, filename = paste0(outdir, "variance_explained_bar_rb_genes.png"), width = 10, height = 4)
ggsave(bar_proportions_rb, filename = paste0(outdir, "variance_explained_bar_rb_genes.pdf"), width = 10, height = 4)


### Plot Pluripotency Genes ###
pluri_genes <- fread(paste0(dir,"data/pluripotency_genes.tsv"), sep = "\t", col.names = "Gene_ID", header = FALSE)


pluri_genes <- GeneConversion[pluri_genes, on = "Gene_ID"]


icc_dt_pluri_genes <- icc_interaction_plus_dt[pluri_genes,on = c("gene", "Gene_ID")]
icc_dt_pluri_genes$grp <- factor(icc_dt_pluri_genes$grp, levels = rev(vars))
icc_dt_pluri_genes <- icc_dt_pluri_genes[data.table(gene = icc_dt_pluri_genes[grp == "Residual"][order(percent)]$gene), on = "gene"]
icc_dt_pluri_genes$Gene_ID <- factor(icc_dt_pluri_genes$Gene_ID, levels = unique(icc_dt_pluri_genes$Gene_ID))
 


pPluri_Genes_Cont <- ggplot() +
						geom_bar(data = icc_dt_pluri_genes, aes(Gene_ID, percent, fill = grp), position = "stack", stat = "identity", alpha = 0.75) +
						theme_classic() +
						# facet_wrap(Gene_ID ~ ., nrow = 3) +
						scale_fill_manual(values = var_colors) +
						theme(plot.title = element_text(hjust = 0.5),
                            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
						ylab("Percent Gene Expression Variance Explained") +
                        ggtitle("Variance Explained of\nPluripotency Genes") +
						theme(axis.title.x=element_blank())

ggsave(pPluri_Genes_Cont, filename = paste0(outdir, "Pluripotent_Gene_Variable_Contributions.png"), width = 6, height = 4.5)
ggsave(pPluri_Genes_Cont, filename = paste0(outdir, "Pluripotent_Gene_Variable_Contributions.pdf"), width = 6, height = 4.5)



##### check for variance explained for eQTL genes (from Kilpinen et al) that are #####
eqtls <- fread("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/KilpinenOverlap/gene_snp_list.tsv", sep = "\t")


eqtls_icc <- icc_interaction_plus_dt[unique(eqtls[,"gene"]), on = "gene"]
eqtls_icc$grp <- factor(eqtls_icc$grp, levels = rev(vars))
eqtls_icc <- eqtls_icc[data.table(gene = eqtls_icc[grp == "Residual"][order(percent)]$gene), on = "gene"]
eqtls_icc$Gene_ID <- factor(eqtls_icc$Gene_ID, levels = unique(eqtls_icc$Gene_ID))
 

group_size_eqtl  <- data.table(table(eqtls_icc$grp))
colnames(group_size_eqtl) <- c("grp", "size")
group_size_eqtl$grp_size <- paste0(group_size_eqtl$grp, "\nN = ", group_size_eqtl$size)

eqtls_icc <- group_size_eqtl[eqtls_icc, on = "grp"]
eqtls_icc$grp_size <- factor(eqtls_icc$grp_size, levels = c(unique(grep("^Line\nN = ", eqtls_icc$grp_size, value = TRUE)), 
                                                                                       unique(grep("^Village\nN = ",  eqtls_icc$grp_size, value = TRUE)), 
                                                                                       unique(grep("^Site\nN = ",  eqtls_icc$grp_size, value = TRUE)), 
                                                                                       unique(grep("^Cryopreserved\nN =",  eqtls_icc$grp_size, value = TRUE)), 
                                                                                       unique(grep("^Replicate\nN =",  eqtls_icc$grp_size, value = TRUE)), 
                                                                                       unique(grep("Line:Village\nN =",  eqtls_icc$grp_size, value = TRUE)), 
                                                                                       unique(grep("Line:Site\nN =",  eqtls_icc$grp_size, value = TRUE)), 
                                                                                       unique(grep("Line:Cryopreserved\nN =",  eqtls_icc$grp_size, value = TRUE)), 
                                                                                       unique(grep("Village:Site\nN =",  eqtls_icc$grp_size, value = TRUE)), 
                                                                                       unique(grep("Village:Cryopreserved\nN =",  eqtls_icc$grp_size, value = TRUE)), 
                                                                                       unique(grep("Replicate:Village\nN ",  eqtls_icc$grp_size, value = TRUE)), 
                                                                                       unique(grep("Replicate:Line\nN ",  eqtls_icc$grp_size, value = TRUE)), 
                                                                                       unique(grep("Replicate:Site\nN ",  eqtls_icc$grp_size, value = TRUE)), 
                                                                                       unique(grep("Replicate:Cryopreserved\nN ",  eqtls_icc$grp_size, value = TRUE)), 
                                                                                       unique(grep("Residual\nN = ", eqtls_icc$grp_size, value = TRUE))))


eqtls_icc_1pct <- eqtls_icc


for (ensg in unique(eqtls_icc$gene)){
    if (!any(eqtls_icc[gene == ensg & grp != "Residual"]$percent > 1)){
        eqtls_icc_1pct <- eqtls_icc_1pct[gene != ensg]
    }
}



eqtls_icc_1pct_grouped_list <- list()


for (ensg in unique(eqtls_icc_1pct$gene)){
    group <- eqtls_icc_1pct[gene == ensg & grp != "Residual"][which.max(percent)]$grp
    eqtls_icc_1pct_grouped_list[[group]][[ensg]] <- eqtls_icc_1pct[gene == ensg]
}

eqtls_icc_1pct_grouped <- lapply(eqtls_icc_1pct_grouped_list, function(x) do.call(rbind, x))
eqtls_icc_1pct_grouped <- lapply(names(eqtls_icc_1pct_grouped), function(x){
    eqtls_icc_1pct_grouped[[x]]$largest_contributor <- x
    return(eqtls_icc_1pct_grouped[[x]])
})

eqtls_icc_1pct_grouped_dt <- do.call(rbind, eqtls_icc_1pct_grouped)

eqtls_icc_1pct_grouped_dt$largest_contributor <- factor(eqtls_icc_1pct_grouped_dt$largest_contributor, levels = vars)



pPluri_Genes_largest_Cont_eqtl <- ggplot() +
						geom_bar(data = eqtls_icc_1pct_grouped_dt, aes(Gene_ID, percent, fill = factor(grp, levels = rev(vars))), position = "stack", stat = "identity", alpha = 0.75) +
						theme_classic() +
						facet_grid(. ~ largest_contributor, scales = "free_x", space = "free_x") +
						scale_fill_manual(values = var_colors) +
						ylab("Percent Gene Expression Variance Explained") +
						theme(axis.title.x=element_blank(),
                            axis.text.x = element_blank(),
                            panel.spacing.x=unit(0, "lines"),
                            axis.ticks.x = element_blank()) +
                        geom_hline(yintercept = 1, linetype = "dashed") 
                        # scale_y_discrete(expand = c(0.03, 0))
                        # scale_x_discrete(expand = c(0.03, 0))


ggsave(pPluri_Genes_largest_Cont_eqtl, filename = paste0(outdir, "eQTL_Genes_Variance_Contributions_1pct_largest_cont.png"), width = 10, height = 4)
ggsave(pPluri_Genes_largest_Cont_eqtl, filename = paste0(outdir, "eQTL_Genes_Variance_Contributions_1pct_largest_cont.pdf"), width = 10, height = 4)







# ### Count number of each chromosome/gene category type 
# ### numbers are the numbers of that category that where a significant percent of variance is explained by this variable
# ### percent of that category that where a significant percent of variance is explained by this variable

# x_number <- lapply(genes_list, function(x){
#     length(which(x %in% X_chromosome_genes$gene_id))
# })

# x_percent <- lapply(genes_list, function(x){
#     length(which(x %in% X_chromosome_genes$gene_id))/length(X_chromosome_genes$gene_id)
# })

# y_number <- lapply(genes_list, function(x){
#     length(which(x %in% Y_chromosome_genes$gene_id))
# })

# y_percent <- lapply(genes_list, function(x){
#     length(which(x %in% Y_chromosome_genes$gene_id))/length(Y_chromosome_genes$gene_id)
# })

# mt_number <- lapply(genes_list, function(x){
#     length(which(x %in% MT_chromosome_genes$gene_id))
# })

# mt_percent <- lapply(genes_list, function(x){
#     length(which(x %in% MT_chromosome_genes$gene_id))/length(MT_chromosome_genes$gene_id)
# })

# rb_number <- lapply(genes_list, function(x){
#     length(which(x %in% RbGeneList$ENSG))
# })

# rb_percent <- lapply(genes_list, function(x){
#     length(which(x %in% RbGeneList$ENSG))/length(RbGeneList$ENSG)
# })



# ### for > 1% var explained
# ### Count number of each chromosome/gene category type 
# x_number_1 <- lapply(vars, function(group){
#     genes <- icc_dt[grp == group][rev(order(percent_round))][percent > 1]$gene
#     length(which(genes %in% X_chromosome_genes$gene_id))
# })
# names(x_number_1) <- vars

# x_percent_1 <- lapply(vars, function(group){
#     genes <- icc_dt[grp == group][rev(order(percent_round))][percent > 1]$gene
#     length(which(genes %in% X_chromosome_genes$gene_id))/length(X_chromosome_genes$gene_id)
# })
# names(x_percent_1) <- vars

# y_number_1 <- lapply(vars, function(group){
#     genes <- icc_dt[grp == group][rev(order(percent_round))][percent > 1]$gene
#     length(which(genes %in% Y_chromosome_genes$gene_id))
# })
# names(y_number_1) <- vars

# y_percent_1 <- lapply(vars, function(group){
#     genes <- icc_dt[grp == group][rev(order(percent_round))][percent > 1]$gene
#     length(which(genes %in% Y_chromosome_genes$gene_id))/length(Y_chromosome_genes$gene_id)
# })
# names(y_percent_1) <- vars

# mt_number_1 <- lapply(vars, function(group){
#     genes <- icc_dt[grp == group][rev(order(percent_round))][percent > 1]$gene
#     length(which(genes %in% MT_chromosome_genes$gene_id))
# })
# names(mt_number_1) <- vars

# mt_percent_1 <- lapply(vars, function(group){
#     genes <- icc_dt[grp == group][rev(order(percent_round))][percent > 1]$gene
#     length(which(genes %in% MT_chromosome_genes$gene_id))/length( MT_chromosome_genes$gene_id)
# })
# names(mt_percent_1) <- vars

# rb_number_1 <- lapply(vars, function(group){
#     genes <- icc_dt[grp == group][rev(order(percent_round))][percent > 1]$gene
#     length(which(genes %in% RbGeneList$ENSG))
# })
# names(rb_number_1) <- vars

# rb_percent_1 <- lapply(vars, function(group){
#     genes <- icc_dt[grp == group][rev(order(percent_round))][percent > 1]$gene
#     length(which(genes %in% RbGeneList$ENSG))/length(RbGeneList$ENSG)
# })
# names(rb_percent_1) <- vars


### Pathway analysis ###
geneList <- lapply(genes_list, function(x){
    tmp <- bitr(x, fromType = "ENSEMBL",
        toType = c("ENTREZID"),
        OrgDb = org.Hs.eg.db)$ENTREZID
    tmp[!is.na(tmp)]
})


gg <- list()
kk <- list()



df = as.data.frame(org.Hs.egGO)
go_gene_list = unique(sort(df$gene_id))

dfk = as.data.frame(org.Hs.egPATH)
kegg_gene_list = unique(sort(dfk$gene_id))



for (group in vars){
    kk[[group]] <- enrichKEGG(gene = geneList[[group]],
                    universe = geneList[[group]],
                    organism     = 'hsa',
                    pvalueCutoff = 0.05,
                    keyType = 'ncbi-geneid')

    gg[[group]] <- groupGO(gene = geneList[[group]],
                OrgDb = org.Hs.eg.db,
                readable = TRUE)

}


hsGO <- godata('org.Hs.eg.db', ont="MF")


sim_results <- list()

vars <- c("Line", "Village", "Cryopreserved",  "Replicate","Village:Line", "Line:Cryopreserved", "Village:Cryopreserved", "Replicate:Village", "Replicate:Line", "Replicate:Cryopreserved")

for (group1 in vars){
    print(group1)
    for (group2 in vars[(grep(paste0("^",group1, "$"), vars) + 1): length(vars)]){
        print(group2)
        sim_results[[group1]][[group2]] <- clusterSim(geneList[[group1]], geneList[[group2]], semData=hsGO, measure="Wang", combine="BMA")
    }
}





# genes2rerun <- c(character())

# for (g in unique(icc_dt$gene)){
#     print(g)
#     genes2rerun <- c(genes2rerun, as.character(unique(icc_dt[gene == g][P > 0.05/(nrow(icc_dt[gene == g])-1)]$gene)))
# }


# for (gene in genes2rerun){
#     print(gene)
#     # unlink(paste0(icc_dir,gene,"_icc.rds"))
#     print(file.exists(paste0(icc_dir,gene,"_icc.rds")))
#     unlink(paste0("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/variance_partition_post_review/gene_separated/fit_models/",gene,"_fitted_models.rds"))
#     unlink(paste0("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/variance_partition_post_review/gene_separated/residuals4qtl/",gene,"_residuals4qtl.rds"))
# }


