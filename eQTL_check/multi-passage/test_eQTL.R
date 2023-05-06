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


cell_line_colors <- c("TOB0199" = "#f44336", "TOB0220" = "#e81f63", "MBE2900" = "#9c27b0", "MBE0953" = "#673ab7", "TOB0421" = "#3f51b5", "MBE2817" = "#2096f3","FSA0004" = "#2096f3", "TOB0435" = "#009688", "WAB0004" = "#4caf50", "IST1877" = "#8bc34a", "WAB0103" = "#cddc39", "TOB0198" = "#ffeb3b", "TOB0205" = "#ffc108", "WAB0038" = "#ff9801", "IST3323" = "#ff5723" , "FSA0001" = "#795548", "180N" = "#9e9e9e", "166" = "#607d8b")

##### Read in dataframe of pairs to test #####
snp_gene_pairs <- fread(bed)




##### Subset the correct  snp #####
snp_gene_pairs_subset <- snp_gene_pairs[gene_id == ensg]



##### Read in data #####
residuals <- readRDS(paste0("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/Variance/variance_integratedSCT/gene_separated/residuals4qtl/", ensg, "_residuals4qtl.rds"))
residual_dt <- data.table(residuals)
colnames(residual_dt) <- c("residual")
residual_dt$Barcode <- names(residuals)

vcf <- read.vcfR("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/multi-passage/deboever_finalized_snps.vcf")
meta <- fread("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/multi-passage/cell_meta.tsv", sep = "\t")
data <- meta[residual_dt, on = "Barcode"]
data$Barcode <- NULL




### Summarize data by site and line and replicatae
data_sum <- data[, .(residual=mean(residual)), by = .(Passage, Line)]
data_sum$Line <- gsub("^0_", "", data_sum$Line) %>%
                        gsub("D-", "", .) %>%
                            gsub("\\.\\d\\.", "", .) %>%
                                gsub("N-", "", .) %>%
                                    gsub("-P36", "", .)  %>%
                                        gsub("-", "", .)


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


colnames(snp_dt) <- gsub("^0_", "", colnames(snp_dt)) %>%
                        gsub("D-", "", .) %>%
                            gsub("\\.\\d\\.", "", .) %>%
                                gsub("N-", "", .) %>%
                                    gsub("-P36", "", .)  %>%
                                        gsub("-", "", .)

    
## Get specific SNP ##
# snp_dt$ID <- paste0(vcf@fix[,'CHROM'],":", vcf@fix[,'POS'],":", vcf@fix[,'REF'],":", vcf@fix[,'ALT'])
snp_dt$ID <- vcf@fix[,'ID']

snp_dt_subset <- snp_dt[ID %in% snp_gene_pairs_subset$ID]


### Add SNP data to data data.table ###
snp_dt_subset_long <- melt(snp_dt_subset, id.vars = c("ID"), measure.vars = colnames(snp_dt_subset)[1:(ncol(snp_dt_subset) - 1)])
snp_dt_subset_long$variable <- gsub("^0_", "", snp_dt_subset_long$variable) %>%
                        gsub("D-", "", .) %>%
                            gsub("\\.\\d\\.", "", .) %>%
                                gsub("N-", "", .) %>%
                                    gsub("-P36", "", .)  %>%
                                        gsub("-", "", .)



##### Make dataframe with deboever results, snp genotypes and our beta + pvalue #####
results_dt <- snp_gene_pairs_subset[,c("chrom", "start", "end", "rsid", "maf", "stat", "pvalue", "beta", "sebeta", "gene_id", "gene_name", "ref","alt", "ID", "REF", "ALT")]
colnames(results_dt)[5:13] <- paste0(colnames(results_dt)[5:13], "_deboever")
results_dt <- results_dt[snp_dt_subset, on = "ID"]

results_dt$dataset_beta <- as.numeric(NA)
results_dt$dataset_beta_se <- as.numeric(NA)
results_dt$dataset_z <- as.numeric(NA)
results_dt$dataset_p <- as.numeric(NA)
results_dt$direction <- as.character(NA)


for (snp in unique(snp_dt_subset_long$ID)){

    data_sum_snp <- snp_dt_subset_long[ID == snp][data_sum, on = c("variable" = "Line")]

    columns <- c("ID",  unique(snp_dt_subset_long$variable), "REF", "ALT")
    gt_long <- melt(results_dt[ID == snp][,..columns], measure.vars =  unique(snp_dt_subset_long$variable))
    gt_long$Genotype <- ifelse(gt_long$value == 0, paste0(gt_long$REF, "/", gt_long$REF), 
                            ifelse(gt_long$value == 1, paste0(gt_long$REF, "/", gt_long$ALT),
                                ifelse(gt_long$value == 2, paste0(gt_long$ALT, "/", gt_long$ALT), NA)))

    gt_long <- unique(gt_long)

    data_sum_snp <- gt_long[data_sum_snp, on = c("value", "variable", "ID")]

    data_sum_snp$Genotype <- factor(data_sum_snp$Genotype, levels = unique(data_sum_snp[order(value)]$Genotype))


    ##### Check for beta #####
    base_model <- glmmTMB(residual ~ 1, data = data_sum_snp, REML = TRUE)
    snp_model <- glmmTMB(as.numeric(residual) ~ as.numeric(value), data = data_sum_snp, REML = TRUE)

    anova_results <- anova(base_model, snp_model)

    results_dt[ID == snp]$dataset_beta <- summary(snp_model)$coefficients$cond[2,1]
    results_dt[ID == snp]$dataset_beta_se <- summary(snp_model)$coefficients$cond[2,2]
    results_dt[ID == snp]$dataset_z <- summary(snp_model)$coefficients$cond[2,3]
    results_dt[ID == snp]$dataset_p <- anova_results$`Pr(>Chisq)`[2]

    if (results_dt[ID == snp]$ref_deboever == results_dt[ID == snp]$REF & results_dt[ID == snp]$alt_deboever == results_dt[ID == snp]$ALT) {
        if ((results_dt[ID == snp]$beta_deboever * results_dt[ID == snp]$dataset_beta) > 0){
            results_dt[ID == snp]$direction <- "match"
        } else {
            results_dt[ID == snp]$direction <- "opposite"
        }
    } else if (results_dt[ID == snp]$ref_deboever == results_dt[ID == snp]$ALT & results_dt[ID == snp]$REF == results_dt[ID == snp]$alt_deboever) {
        if ((results_dt[ID == snp]$beta_deboever * results_dt[ID == snp]$dataset_beta) < 0){
            results_dt[ID == snp]$direction <- "match"
        } else {
            results_dt[ID == snp]$direction <- "opposite"
        }
    } else {
        results_dt[ID == snp]$direction <- "different_snp"
    }

    ###### Make a figure of the results
    ### Different shapes for location
    ### Different fill for village and uniculture
    ### Color by line
    if (!results_dt[ID == snp]$direction ==  "different_snp"){

        data_sum_snp$value <- factor(data_sum_snp$value, levels = sort(unique(data_sum_snp$value)))

        labels <- levels(data_sum_snp$Genotype)
        names(labels) <- as.numeric(levels(data_sum_snp$value))

        plot <- ggplot(data_sum_snp, aes(value, residual)) +
            geom_point(aes(color = variable)) +
            theme_classic() +
            scale_color_manual(values = cell_line_colors) +
            scale_fill_manual(values = cell_line_colors) +
            geom_smooth(aes(as.numeric(value), residual), position = "identity",method = "lm", color = "black", se=FALSE) +
            scale_x_discrete(labels=labels) +
            ylab("Normalized Expression") +
            ggtitle(paste0(results_dt[ID == snp]$rsid, "-", results_dt[ID == snp]$gene_name_deboever, " eQTL"),
                    subtitle = paste0("beta = ", round(results_dt[ID == snp]$dataset_beta, 3))) +
            xlab(paste0(results_dt[ID == snp]$rsid, '\nGenotype')) +
            theme(plot.title = element_text(hjust = 0.5),
                    plot.subtitle = element_text(hjust = 0.5)) +
            labs(color = "iPSC Line", shape = "Site & Village", fill = NULL)

        ggsave(plot, filename = paste0(outdir,"plots/",ensg, "_", snp,"_deboever_eQTL_results.png"), width = 5, height = 3.5)
        ggsave(plot, filename = paste0(outdir,"plots/",ensg, "_", snp,"_deboever_eQTL_results.pdf"), width = 5, height = 3.5)
    }
}

fwrite(results_dt, paste0(outdir,"beds/",ensg, "_deboever_eQTL_results.bed"), sep = "\t")




