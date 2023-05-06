### Reason: test eQTLs identified in Kilpinen et al in our dataset
### Author: Drew Neavin
### Date: 1 April, 2022


library(data.table)
library(Seurat)
library(vcfR)
library(glmmTMB)
library(dplyr)
library(ggplot2)



##### Bring in variables #####
### Bring in arguments
args <- commandArgs(trailingOnly = TRUE)
ensg <- args[1]
outdir <- args[2]
bed <- args[3]


cell_line_colors <- c("FSA0006" = "#F79E29", "MBE1006" = "#9B2C99", "TOB0421"= "#35369C")


##### Read in dataframe of pairs to test #####
snp_gene_pairs <- fread(bed)
snp_gene_pairs$gene_id <- gsub("\\..+", "", snp_gene_pairs$gene_id)




##### Subset the correct  snp #####
snp_gene_pairs_subset <- snp_gene_pairs[gene_id == ensg]



##### Read in data #####
residuals <- readRDS(paste0("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/variance_partition_post_review_cryo2/village/gene_separated/residuals4qtl/", ensg, "_residuals4qtl.rds"))
residual_dt <- data.table(residuals)
colnames(residual_dt) <- c("residual")
residual_dt$Barcode <- names(residuals)

vcf <- read.vcfR("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/uni_village/cryopreserved/village/deboever_finalized_snps.recode.vcf")
meta <- fread("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/uni_village/cryopreserved/cell_meta.tsv", sep = "\t")
data <- meta[residual_dt, on = "Barcode"]
data$Barcode <- NULL




### Summarize data by site and line and replicatae
data_sum <- data[, .(residual=mean(residual)), by = .(Village, Line, Cryopreserved, Replicate)]



##### Get the snp from the vcf #####
## GT ##
snp_dt <- data.table(extract.gt(element = "GT",vcf, IDtoRowNames = F))


if (!all(colSums(is.na(snp_dt)) == nrow(snp_dt))){
    message("Found GT genotype format in cluster vcf. Will use that metric for cluster correlation.")
    format_clust = "GT"

    if (any(grepl("\\|",snp_dt[,1]))){
        separator = "|"
        message("Detected | separator for GT genotype format in cluster vcf")
    } else if (any(grepl("/",snp_dt[,1]))) {
        separator = "/"
        message("Detected / separator for GT genotype format in cluster vcf")
    } else {
        format_clust = NA
        message("Can't identify a separator for the GT field in cluster vcf, moving on to using GP.")
    }
    if (!is.na(format_clust)){
        snp_dt <- data.table(as_tibble(lapply(snp_dt, function(x) {gsub(paste0("0\\",separator,"0"),0, x)}) %>%
                                lapply(., function(x) {gsub(paste0("0\\",separator,"1"),1, x)}) %>%
                                lapply(., function(x) {gsub(paste0("1\\",separator,"0"),1, x)}) %>%
                                lapply(., function(x) {gsub(paste0("1\\",separator,"1"),2, x)})))
    }
}


colnames(snp_dt) <- gsub("36_TOB00421_i_E8", "TOB0421", colnames(snp_dt)) %>%
                        gsub("22_FSA", "FSA0006", .) %>%
                        gsub("29_MBE", "MBE1006", .)

    
## Get specific SNP ##
snp_dt$ID_ref_alt <- paste0(vcf@fix[,'CHROM'],":", vcf@fix[,'POS'],"_", vcf@fix[,'REF'],"_", vcf@fix[,'ALT'])

snp_dt_subset <- snp_dt[ID_ref_alt %in% snp_gene_pairs_subset$ID_ref_alt]



### Add SNP data to data data.table ###
snp_dt_subset_long <- melt(snp_dt_subset, id.vars = c("ID_ref_alt"), measure.vars = c("TOB0421", "FSA0006", "MBE1006"))



##### Make dataframe with deboever results, snp genotypes and our beta + pvalue #####
results_dt <- snp_gene_pairs_subset[,c("chrom", "start", "end", "rsid", "maf", "stat", "pvalue", "beta", "sebeta", "gene_id", "gene_name", "ref","alt", "ID", "REF", "ALT", "ID_ref_alt")]
colnames(results_dt) <- c(colnames(results_dt)[1:4],paste0(colnames(results_dt)[5:13], "_deboever"), colnames(results_dt)[14:17])
results_dt <- results_dt[snp_dt_subset, on = "ID_ref_alt"]

results_dt$dataset_beta <- as.numeric(NA)
results_dt$dataset_beta_se <- as.numeric(NA)
results_dt$dataset_z <- as.numeric(NA)
results_dt$dataset_p <- as.numeric(NA)
results_dt$direction <- as.character(NA)


for (snp in unique(snp_dt_subset_long$ID_ref_alt)){

    data_sum_snp <- snp_dt_subset_long[ID_ref_alt == snp][data_sum, on = c("variable" = "Line")]

    gt_long <- melt(results_dt[ID_ref_alt == snp][,c("ID_ref_alt", "TOB0421", "FSA0006", "MBE1006", "REF", "ALT")], measure.vars = c("TOB0421", "FSA0006", "MBE1006"))
    gt_long$Genotype <- ifelse(gt_long$value == 0, paste0(gt_long$REF, "/", gt_long$REF), 
                            ifelse(gt_long$value == 1, paste0(gt_long$REF, "/", gt_long$ALT),
                                ifelse(gt_long$value == 2, paste0(gt_long$ALT, "/", gt_long$ALT), NA)))

    gt_long <- unique(gt_long)

    data_sum_snp <- gt_long[data_sum_snp, on = c("value", "variable", "ID_ref_alt")]

    data_sum_snp$Genotype <- factor(data_sum_snp$Genotype, levels = unique(data_sum_snp[order(value)]$Genotype))


    ##### Check for beta #####
    base_model <- glmmTMB(residual ~ 1, data = data_sum_snp, REML = TRUE)
    snp_model <- glmmTMB(as.numeric(residual) ~ as.numeric(value), data = data_sum_snp, REML = TRUE)

    anova_results <- anova(base_model, snp_model)

    results_dt[ID_ref_alt == snp]$dataset_beta <- summary(snp_model)$coefficients$cond[2,1]
    results_dt[ID_ref_alt == snp]$dataset_beta_se <- summary(snp_model)$coefficients$cond[2,2]
    results_dt[ID_ref_alt == snp]$dataset_z <- summary(snp_model)$coefficients$cond[2,3]
    results_dt[ID_ref_alt == snp]$dataset_p <- anova_results$`Pr(>Chisq)`[2]

    if (results_dt[ID_ref_alt == snp]$ref_deboever == results_dt[ID_ref_alt == snp]$REF & results_dt[ID_ref_alt == snp]$alt_deboever == results_dt[ID_ref_alt == snp]$ALT) {
        if ((results_dt[ID_ref_alt == snp]$beta_deboever * results_dt[ID_ref_alt == snp]$dataset_beta) > 0){
            results_dt[ID_ref_alt == snp]$direction <- "match"
        } else {
            results_dt[ID_ref_alt == snp]$direction <- "opposite"
        }
    } else if (results_dt[ID_ref_alt == snp]$ref_deboever == results_dt[ID_ref_alt == snp]$ALT & results_dt[ID_ref_alt == snp]$REF == results_dt[ID_ref_alt == snp]$alt_deboever) {
        if ((results_dt[ID_ref_alt == snp]$beta_deboever * results_dt[ID_ref_alt == snp]$dataset_beta) < 0){
            results_dt[ID_ref_alt == snp]$direction <- "match"
        } else {
            results_dt[ID_ref_alt == snp]$direction <- "opposite"
        }
    } else {
        results_dt[ID_ref_alt == snp]$direction <- "different_snp"
    }

    ###### Make a figure of the results
    ### Different shapes for location
    ### Different fill for village and uniculture
    ### Color by line
    if (!results_dt[ID_ref_alt == snp]$direction ==  "different_snp"){
        shapes <- c(0,15,1,16,2,17)
        names(shapes) <- paste0(c("Sydney", "Sydney", "Melbourne", "Melbourne", "Brisbane", "Brisbane"), " ", rep(c("Uni-culture","Village")))

        data_sum_snp$shapes <- paste0(data_sum_snp$Site, " ", gsub(0, "Uni-culture",gsub(1, "Village", data_sum_snp$Village)))
        data_sum_snp$value <- factor(data_sum_snp$value, levels = sort(unique(data_sum_snp$value)))

        labels <- levels(data_sum_snp$Genotype)
        names(labels) <- as.numeric(levels(data_sum_snp$value))

        plot <- ggplot(data_sum_snp, aes(value, residual)) +
            geom_point(aes(shape = shapes, color = variable)) +
            theme_classic() +
            scale_color_manual(values = cell_line_colors) +
            scale_fill_manual(values = cell_line_colors) +
            scale_shape_manual(values = shapes) +
            geom_smooth(aes(as.numeric(value), residual), position = "identity",method = "lm", color = "black", se=FALSE) +
            scale_x_discrete(labels=labels) +
            ylab("Normalized Expression") +
            ggtitle(paste0(results_dt[ID_ref_alt == snp]$rsid, "-", results_dt[ID_ref_alt == snp]$gene_name_deboever, " eQTL"),
                    subtitle = paste0("beta = ", round(results_dt[ID_ref_alt == snp]$dataset_beta, 3))) +
            xlab(paste0(results_dt[ID_ref_alt == snp]$rsid, '\nGenotype')) +
            theme(plot.title = element_text(hjust = 0.5),
                    plot.subtitle = element_text(hjust = 0.5)) +
            labs(color = "iPSC Line", shape = "Site & Village", fill = NULL)

        ggsave(plot, filename = paste0(outdir,"plots/",ensg, "_", snp,"_deboever_eQTL_results.png"), width = 5, height = 3.5)
        ggsave(plot, filename = paste0(outdir,"plots/",ensg, "_", snp,"_deboever_eQTL_results.pdf"), width = 5, height = 3.5)
    }
}

fwrite(results_dt, paste0(outdir,"beds/",ensg, "_deboever_eQTL_results.bed"), sep = "\t")




