### Author: Drew neavin
### Date: 25 June, 2022
### Rational: Reviewer asked if eQTL effects are consistently higher in a village so need to look at village*snp effect for eQTLs


##### Read in Libraries #####
library(data.table)
library(tidyverse)
library(Seurat)
library(vcfR)
library(pkgcond)
library(glmmTMB)



##### Set up directories #####
args <- commandArgs(trailingOnly = TRUE)
ensg <- as.character(args[1])
outdir <- args[2]
bed <- args[3]
datadir <- args[4]



##### Read in data #####
### significant effects for all genes dataframe ###
effects_dt <- fread("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/variance_partition_post_review2/combined/effect_results.tsv")

### Seurat object ###
seurat <- readRDS("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/variance_partition_post_review/seurat_integrated_noncryo_1pct_expressing.rds")


### SNP vcf ###
vcf <- read.vcfR("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/KilpinenOverlap/finalized_snps.recode.vcf")


##### Read in dataframe of pairs to test #####
snp_gene_pairs <- fread(bed)
snp_gene_pairs$gene_id <- gsub("\\..+", "", snp_gene_pairs$gene_id)


##### Subset the correct snp #####
snp_gene_pairs_subset <- snp_gene_pairs[gene_id == ensg]


##### Make model dataframe #####
df_hier_unscale <- data.frame("Expression" = seurat[["SCT"]]@scale.data[ensg,], "Village" = as.factor(ifelse(seurat@meta.data$Time == "Baseline", 0, 1)), "Line" = seurat@meta.data$Final_Assignment, "Replicate" = as.factor(gsub("[A-Z][a-z]+", "", seurat@meta.data$MULTI_ID)), "Site" = seurat$Location)
colnames(df_hier_unscale)[1] <- "Expression"




##### Get the Correct Variables for Modeling #####
variables <- unique(effects_dt[gene == ensg]$grp)
variables <- variables[!(variables %in% c("Line", "Line:Village", "Residual"))]



##### Fit effects except line and village*line effects ######
model_all <- as.formula(paste0("Expression ~ (1|", paste0(variables, collapse = ") + (1|"), ")"))
model_glmmtmb <- suppress_warnings(glmmTMB(formula = noquote(model_all), data = df_hier_unscale, REML = TRUE), "giveCsparse")



###### Summarize data to site x replicate x village #####
data <- data.table(df_hier_unscale)
data$Residuals <- resid(model_glmmtmb)
data_sum <- data[, .(residual=mean(Residuals)), by = .(Village, Line, Site, Replicate)]



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
snp_dt_subset_long <- melt(snp_dt_subset, id.vars = c("ID_ref_alt"), measure.vars = c("TOB0421", "FSA0006", "MBE1006"))



##### Make dataframe with results, snp genotypes and our beta + pvalue #####
results_dt <- snp_gene_pairs_subset[,c("ID", "REF", "ALT", "ID_ref_alt")]
results_dt <- results_dt[snp_dt_subset, on = "ID_ref_alt"]

results_dt$snp_beta <- as.numeric(NA)
results_dt$snp_beta_se <- as.numeric(NA)
results_dt$snp_z <- as.numeric(NA)
results_dt$snp_p <- as.numeric(NA)
results_dt$snpXvillage_beta <- as.numeric(NA)
results_dt$snpXvillage_beta_se <- as.numeric(NA)
results_dt$snpXvillage_z <- as.numeric(NA)
results_dt$snpXvillage_p <- as.numeric(NA)




for (snp in unique(snp_dt_subset_long$ID_ref_alt)){

    ### Add SNP data to data data.table ###
    data_sum_snp <- snp_dt_subset_long[ID_ref_alt == snp][data_sum, on = c("variable" = "Line")]

    gt_long <- melt(results_dt[ID_ref_alt == snp][,c("ID_ref_alt", "TOB0421", "FSA0006", "MBE1006", "REF", "ALT")], measure.vars = c("TOB0421", "FSA0006", "MBE1006"))
    gt_long$Genotype <- ifelse(gt_long$value == 0, paste0(gt_long$REF, "/", gt_long$REF), 
                            ifelse(gt_long$value == 1, paste0(gt_long$REF, "/", gt_long$ALT),
                                ifelse(gt_long$value == 2, paste0(gt_long$ALT, "/", gt_long$ALT), NA)))

    gt_long <- unique(gt_long)

    data_sum_snp <- gt_long[data_sum_snp, on = c("value", "variable", "ID_ref_alt")]

    data_sum_snp$Genotype <- factor(data_sum_snp$Genotype, levels = unique(data_sum_snp[order(value)]$Genotype))


    ##### Fit SNP => store resids #####
    ##### Check for beta #####
    base_model <- glmmTMB(residual ~ 1, data = data_sum_snp, REML = TRUE)
    snp_model <- glmmTMB(as.numeric(residual) ~ as.numeric(value), data = data_sum_snp, REML = TRUE)

    anova_results <- anova(base_model, snp_model)

    results_dt[ID_ref_alt == snp]$snp_beta <- summary(snp_model)$coefficients$cond[2,1]
    results_dt[ID_ref_alt == snp]$snp_beta_se <- summary(snp_model)$coefficients$cond[2,2]
    results_dt[ID_ref_alt == snp]$snp_z <- summary(snp_model)$coefficients$cond[2,3]
    results_dt[ID_ref_alt == snp]$snp_p <- anova_results$`Pr(>Chisq)`[2]


    ##### resids fit SNP*village #####
    data_sum_snp$residual2 <- resid(snp_model)

    snpXvillage_model <- glmmTMB(as.numeric(residual) ~ as.numeric(value)*Village, data = data_sum_snp, REML = TRUE)

    anova_results2 <- anova(snp_model, snpXvillage_model)


    results_dt[ID_ref_alt == snp]$snpXvillage_beta <- summary(snpXvillage_model)$coefficients$cond[2,1]
    results_dt[ID_ref_alt == snp]$snpXvillage_beta_se <- summary(snpXvillage_model)$coefficients$cond[2,2]
    results_dt[ID_ref_alt == snp]$snpXvillage_z <- summary(snpXvillage_model)$coefficients$cond[2,3]
    results_dt[ID_ref_alt == snp]$snpXvillage_p <- anova_results2$`Pr(>Chisq)`[2]

    fwrite(data_sum_snp, paste0(datadir, ensg, "_",snp,"_dataframe"), sep = "\t")
}

fwrite(results_dt, paste0(outdir, ensg, "_snpXvillage_interactions.tsv"), sep = "\t")


