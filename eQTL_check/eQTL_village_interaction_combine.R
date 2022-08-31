library(data.table)
library(dplyr)
library(ggplot2)
library(colorspace)
library(DEGreport)
library(ggpp)


deboever_indir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/eQTL_village_interaction/deboever/results/"
kilpinen_indir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/eQTL_village_interaction/kilpinen/results/"
outdir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/eQTL_village_interaction/combined/"

dir.create(outdir, recursive = TRUE)


##### Read in results #####
### Deboever
deboever_files <- list.files(deboever_indir)


deboever_results_list <- lapply(deboever_files, function(x){
    tmp <- fread(paste0(deboever_indir, x), sep = "\t")
    tmp$ensg <- gsub("_snpXvillage_interactions.tsv", "", x)
    return(tmp)
})



deboever_results_dt <- do.call(rbind, deboever_results_list)
deboever_results_dt$snp_fdr <- p.adjust(deboever_results_dt$snp_p, method = "fdr")
deboever_results_dt$snpXvillage_fdr <- p.adjust(deboever_results_dt$snpXvillage_p, method = "fdr")
deboever_results_dt$study <- "DeBoever"


### Kilpinen
kilpinen_files <- list.files(kilpinen_indir)

kilpinen_results_list <- lapply(kilpinen_files, function(x){
    tmp <- fread(paste0(kilpinen_indir, x), sep = "\t")
    tmp$ensg <- gsub("_snpXvillage_interactions.tsv", "", x)
    return(tmp)
})


kilpinen_results_dt <- do.call(rbind, kilpinen_results_list)
kilpinen_results_dt$snp_fdr <- p.adjust(kilpinen_results_dt$snp_p, method = "fdr")
kilpinen_results_dt$snpXvillage_fdr <- p.adjust(kilpinen_results_dt$snpXvillage_p, method = "fdr")
kilpinen_results_dt$study <- "Kilpinen"


results_dt <- rbind(deboever_results_dt,kilpinen_results_dt)


results_dt[order(snp_fdr)]
results_dt[order(snpXvillage_fdr)]

### See how many snps have all different genotypes for each cell line
results_dt[TOB0421 != FSA0006 & TOB0421 != MBE1006 & FSA0006 != MBE1006] ## 1108 total

fwrite(results_dt, paste0(outdir, "snpXvillage_results.tsv"), sep = "\t")


results_dt$significant <- ifelse(results_dt$snpXvillage_fdr < 0.05, TRUE, FALSE)

results_dt[ensg == "ENSG00000106153"]

results_dt$snpXvillage_sig <- ifelse(results_dt$snpXvillage_fdr < 0.05, "significant", "not significant")


correlation <- list()
correlation_dt <- data.table(study = unique(results_dt$study), rho = as.numeric(NA), p = as.numeric(NA))

for (Study in results_dt$study){
    correlation[[Study]] <- cor.test(results_dt[study == Study]$snp_beta, results_dt[study == Study]$snpXvillage_beta, method = "spearman", exact = TRUE)
    correlation_dt[study == Study, "rho"] <- correlation[[Study]]$estimate[1]
    correlation_dt[study == Study, "p"] <- correlation[[Study]]$p.value
}


correlation_plot <- ggplot(results_dt[study == "DeBoever"]) +
        geom_hline(yintercept=0) +
        geom_vline(xintercept=0) +
        geom_point(aes(snp_beta, snpXvillage_beta)) +
        theme_classic() +
        # facet_wrap(vars(study)) +
        xlab("SNP Effect Size") +
        ylab("SNP x Village Interaction Effect Size") +
        geom_text(data = correlation_dt, aes(x = -3.85, y = 5.5, label = paste0("Rho = ",round(rho,2),"***")))
        #  scale_color_discrete_sequential(palette = "ag_GrnYl", rev = FALSE) +

ggsave(correlation_plot, filename = paste0(outdir, "correlation_plot.png"), height = 3, width = 5.5)



### add number of genes in each quadrant
dim(results_dt[snp_beta < 0 & snpXvillage_beta < 0]) # 2216
dim(results_dt[snp_beta > 0 & snpXvillage_beta > 0]) # 2103

dim(results_dt[snp_beta < 0 & snpXvillage_beta > 0]) # 1155
dim(results_dt[snp_beta > 0 & snpXvillage_beta < 0]) # 1160


dim(results_dt[snp_beta < 0 & snpXvillage_beta < 0 & study == "DeBoever"]) # 1670
dim(results_dt[snp_beta > 0 & snpXvillage_beta > 0 & study == "DeBoever"]) # 1582

dim(results_dt[snp_beta < 0 & snpXvillage_beta > 0 & study == "DeBoever"]) # 546
dim(results_dt[snp_beta > 0 & snpXvillage_beta < 0 & study == "DeBoever"]) # 586


dim(results_dt[snp_beta < 0 & snpXvillage_beta < 0 & study == "Kilpinen"]) # 546
dim(results_dt[snp_beta > 0 & snpXvillage_beta > 0 & study == "Kilpinen"]) # 521

dim(results_dt[snp_beta < 0 & snpXvillage_beta > 0 & study == "Kilpinen"]) # 609
dim(results_dt[snp_beta > 0 & snpXvillage_beta < 0 & study == "Kilpinen"]) # 574


### Sig
dim(results_dt[snp_beta < 0 & snpXvillage_beta < 0 & snpXvillage_sig == "significant"]) # 15
dim(results_dt[snp_beta > 0 & snpXvillage_beta > 0 & snpXvillage_sig == "significant"]) # 18

dim(results_dt[snp_beta < 0 & snpXvillage_beta > 0 & snpXvillage_sig == "significant"]) # 67
dim(results_dt[snp_beta > 0 & snpXvillage_beta < 0 & snpXvillage_sig == "significant"]) # 72


dim(results_dt[snp_beta < 0 & snpXvillage_beta < 0 & study == "DeBoever" & snpXvillage_sig == "significant"]) # 13
dim(results_dt[snp_beta > 0 & snpXvillage_beta > 0 & study == "DeBoever" & snpXvillage_sig == "significant"]) # 14

dim(results_dt[snp_beta < 0 & snpXvillage_beta > 0 & study == "DeBoever" & snpXvillage_sig == "significant"]) # 47
dim(results_dt[snp_beta > 0 & snpXvillage_beta < 0 & study == "DeBoever" & snpXvillage_sig == "significant"]) # 46


dim(results_dt[snp_beta < 0 & snpXvillage_beta < 0 & study == "Kilpinen" & snpXvillage_sig == "significant"]) # 2
dim(results_dt[snp_beta > 0 & snpXvillage_beta > 0 & study == "Kilpinen" & snpXvillage_sig == "significant"]) # 4

dim(results_dt[snp_beta < 0 & snpXvillage_beta > 0 & study == "Kilpinen" & snpXvillage_sig == "significant"]) # 20
dim(results_dt[snp_beta > 0 & snpXvillage_beta < 0 & study == "Kilpinen" & snpXvillage_sig == "significant"]) # 26


### Filter by genes that matched direction in original dataset ###
effect_match_dt <- fread( "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/eQTL_overlap/combined/overlap_results.tsv")

effect_match_dt_sig_match <- effect_match_dt[direction == "match" & fdr < 0.05]


results_dt_sig_match <- results_dt[effect_match_dt_sig_match[,c("study", "ensg", "ID_ref_alt")], on = c("study", "ensg", "ID_ref_alt")]



correlation_sig_match <- list()
correlation_sig_match_dt <- data.table(study = unique(results_dt_sig_match$study), rho = as.numeric(NA), p = as.numeric(NA))

for (Study in results_dt_sig_match$study){
    correlation_sig_match[[Study]] <- cor.test(results_dt_sig_match[study == Study]$snp_beta, results_dt_sig_match[study == Study]$snpXvillage_beta, method = "spearman", exact = TRUE)
    correlation_sig_match_dt[study == Study, "rho"] <- correlation_sig_match[[Study]]$estimate[1]
    correlation_sig_match_dt[study == Study, "p"] <- correlation_sig_match[[Study]]$p.value
}


correlation_sig_match_plot <- ggplot(results_dt_sig_match[study == "DeBoever"], aes(snp_beta, snpXvillage_beta)) +
        geom_hline(yintercept=0) +
        geom_vline(xintercept=0) +
        geom_point() +
        theme_classic() +
        # facet_wrap(vars(study)) +
        xlab("SNP Effect Size") +
        ylab("SNP*Village Interaction\nEffect Size") +
        geom_text(data = correlation_sig_match_dt, aes(x = -4, y = 0.5, label = paste0("Rho = ",round(rho,2),"***"))) +
        stat_quadrant_counts() +
        xlim(-7,3) +
        ylim(-7,3)
        #  scale_color_discrete_sequential(palette = "ag_GrnYl", rev = FALSE) +

ggsave(correlation_sig_match_plot, filename = paste0(outdir, "correlation_sig_match_plot.png"), height = 2.5, width = 2.75)
ggsave(correlation_sig_match_plot, filename = paste0(outdir, "correlation_sig_match_plot.pdf"), height = 2.5, width = 2.75)



### add number of genes in each quadrant
dim(results_dt_sig_match[snp_beta < 0 & snpXvillage_beta < 0]) # 2216
dim(results_dt_sig_match[snp_beta > 0 & snpXvillage_beta > 0]) # 2103

dim(results_dt_sig_match[snp_beta < 0 & snpXvillage_beta > 0]) # 1155
dim(results_dt_sig_match[snp_beta > 0 & snpXvillage_beta < 0]) # 1160


dim(results_dt_sig_match[snp_beta < 0 & snpXvillage_beta < 0 & study == "DeBoever"]) # 1670
dim(results_dt_sig_match[snp_beta > 0 & snpXvillage_beta > 0 & study == "DeBoever"]) # 1582

dim(results_dt_sig_match[snp_beta < 0 & snpXvillage_beta > 0 & study == "DeBoever"]) # 546
dim(results_dt_sig_match[snp_beta > 0 & snpXvillage_beta < 0 & study == "DeBoever"]) # 586


dim(results_dt_sig_match[snp_beta < 0 & snpXvillage_beta < 0 & study == "Kilpinen"]) # 546
dim(results_dt_sig_match[snp_beta > 0 & snpXvillage_beta > 0 & study == "Kilpinen"]) # 521

dim(results_dt_sig_match[snp_beta < 0 & snpXvillage_beta > 0 & study == "Kilpinen"]) # 609
dim(results_dt_sig_match[snp_beta > 0 & snpXvillage_beta < 0 & study == "Kilpinen"]) # 574


### Sig
dim(results_dt_sig_match[snp_beta < 0 & snpXvillage_beta < 0 & snpXvillage_sig == "significant"]) # 15
dim(results_dt_sig_match[snp_beta > 0 & snpXvillage_beta > 0 & snpXvillage_sig == "significant"]) # 18

dim(results_dt_sig_match[snp_beta < 0 & snpXvillage_beta > 0 & snpXvillage_sig == "significant"]) # 67
dim(results_dt_sig_match[snp_beta > 0 & snpXvillage_beta < 0 & snpXvillage_sig == "significant"]) # 72


dim(results_dt_sig_match[snp_beta < 0 & snpXvillage_beta < 0 & study == "DeBoever" & snpXvillage_sig == "significant"]) # 13
dim(results_dt_sig_match[snp_beta > 0 & snpXvillage_beta > 0 & study == "DeBoever" & snpXvillage_sig == "significant"]) # 14

dim(results_dt_sig_match[snp_beta < 0 & snpXvillage_beta > 0 & study == "DeBoever" & snpXvillage_sig == "significant"]) # 47
dim(results_dt_sig_match[snp_beta > 0 & snpXvillage_beta < 0 & study == "DeBoever" & snpXvillage_sig == "significant"]) # 46


dim(results_dt_sig_match[snp_beta < 0 & snpXvillage_beta < 0 & study == "Kilpinen" & snpXvillage_sig == "significant"]) # 2
dim(results_dt_sig_match[snp_beta > 0 & snpXvillage_beta > 0 & study == "Kilpinen" & snpXvillage_sig == "significant"]) # 4

dim(results_dt_sig_match[snp_beta < 0 & snpXvillage_beta > 0 & study == "Kilpinen" & snpXvillage_sig == "significant"]) # 20
dim(results_dt_sig_match[snp_beta > 0 & snpXvillage_beta < 0 & study == "Kilpinen" & snpXvillage_sig == "significant"]) # 26















### Filter by original interaction significant genes ###
effects_dt <- fread("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/variance_partition_post_review2/combined/effect_results.tsv")

interaction_sig_genes <- unique(effects_dt[grp == "Line:Village", "gene"])

results_dt_subset <- results_dt[ensg %in% interaction_sig_genes$gene]

correlation_plot_sub <- ggplot(results_dt_subset, aes(snp_beta, snpXvillage_beta, color = snpXvillage_sig)) +
        geom_point() +
        theme_classic() +
         scale_color_discrete_sequential(palette = "ag_GrnYl", rev = FALSE) +
         facet_wrap(vars(study)) +
         geom_hline(yintercept=0) +
         geom_vline(xintercept=0)
        # scale_color_manual(values = c("grey","black"))

ggsave(correlation_plot_sub, filename = paste0(outdir, "correlation_plot_var_sig.png"), height = 2.25, width = 5.5)



### add number of genes in each quadrant
dim(results_dt_subset[snp_beta < 0 & snpXvillage_beta < 0]) # 202
dim(results_dt_subset[snp_beta > 0 & snpXvillage_beta > 0]) # 232

dim(results_dt_subset[snp_beta < 0 & snpXvillage_beta > 0]) # 304
dim(results_dt_subset[snp_beta > 0 & snpXvillage_beta < 0]) # 323


dim(results_dt_subset[snp_beta < 0 & snpXvillage_beta < 0 & study == "DeBoever"]) # 144
dim(results_dt_subset[snp_beta > 0 & snpXvillage_beta > 0 & study == "DeBoever"]) # 153

dim(results_dt_subset[snp_beta < 0 & snpXvillage_beta > 0 & study == "DeBoever"]) # 199
dim(results_dt_subset[snp_beta > 0 & snpXvillage_beta < 0 & study == "DeBoever"]) # 219


dim(results_dt_subset[snp_beta < 0 & snpXvillage_beta < 0 & study == "Kilpinen"]) # 58
dim(results_dt_subset[snp_beta > 0 & snpXvillage_beta > 0 & study == "Kilpinen"]) # 79

dim(results_dt_subset[snp_beta < 0 & snpXvillage_beta > 0 & study == "Kilpinen"]) # 105
dim(results_dt_subset[snp_beta > 0 & snpXvillage_beta < 0 & study == "Kilpinen"]) # 104



deboever_eqtls <- fread("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/deboever/combined/deboever_overlap_results.tsv", sep = "\t")
deboever_eqtls_dif_snp <- deboever_eqtls[direction == "different_snp"]

results_dt_subset_correct_direction <- results_dt_subset[!(ID_ref_alt %in% deboever_eqtls_dif_snp$gene_id_deboever)]


deboever_results_dt[deboever_eqtls[, .(gene_id_deboever, ID_ref_alt)], on = c("ensg" = "gene_id_deboever", "ID_ref_alt")][ensg == "ENSG00000106153"]



### add number of genes in each quadrant
dim(results_dt_subset_correct_direction[snp_beta < 0 & snpXvillage_beta < 0]) # 144
dim(results_dt_subset_correct_direction[snp_beta > 0 & snpXvillage_beta > 0]) # 153

dim(results_dt_subset_correct_direction[snp_beta < 0 & snpXvillage_beta > 0]) # 199
dim(results_dt_subset_correct_direction[snp_beta > 0 & snpXvillage_beta < 0]) # 219


dim(results_dt_subset_correct_direction[snp_beta < 0 & snpXvillage_beta < 0 & study == "DeBoever"]) # 144
dim(results_dt_subset_correct_direction[snp_beta > 0 & snpXvillage_beta > 0 & study == "DeBoever"]) # 153

dim(results_dt_subset_correct_direction[snp_beta < 0 & snpXvillage_beta > 0 & study == "DeBoever"]) # 199
dim(results_dt_subset_correct_direction[snp_beta > 0 & snpXvillage_beta < 0 & study == "DeBoever"]) # 219


dim(results_dt_subset_correct_direction[snp_beta < 0 & snpXvillage_beta < 0 & study == "Kilpinen"]) # 58
dim(results_dt_subset_correct_direction[snp_beta > 0 & snpXvillage_beta > 0 & study == "Kilpinen"]) # 79

dim(results_dt_subset_correct_direction[snp_beta < 0 & snpXvillage_beta > 0 & study == "Kilpinen"]) # 105
dim(results_dt_subset_correct_direction[snp_beta > 0 & snpXvillage_beta < 0 & study == "Kilpinen"]) # 104





histogram_plot <- ggplot(results_dt, aes(fdr, fill = direction)) +
        geom_histogram(position = "stack", binwidth = 0.05) +
        theme_classic() +
        geom_vline(aes(xintercept=0.05),
            color="red", linetype="dashed", size=1) +
        scale_fill_manual(values = c("grey","black"))
        
ggsave(histogram_plot, filename = paste0(outdir, "histogram_plot.png"))

histogram_plot_log <- ggplot(results_dt, aes(log10(fdr), fill = direction)) +
        geom_histogram(position = "stack", binwidth = 1) +
        theme_classic() +
        geom_vline(aes(xintercept=log10(0.05)),
            color="red", linetype="dashed", size=0.25) +
        scale_fill_manual(values = c("grey","black")) +
        scale_y_continuous(expand = c(0, 0)) +
        ylab("Number of Associations") +
        labs(fill = "eQTL\nEffect\nDirection\nAgreement")

        
ggsave(histogram_plot_log, filename = paste0(outdir, "histogram_plot_log.png"))


popout_histogram_plot_log <- ggplot(results_dt[log10(fdr) < -10], aes(log10(fdr), fill = direction)) +
        geom_histogram(position = "stack", binwidth = 1) +
        theme_classic() +
        scale_fill_manual(values = c("grey","black")) +
        theme(axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            legend.position="none",
            axis.text.x=element_text(size=7),
            axis.text.y=element_text(size=7)) +
        scale_y_continuous(expand = c(0, 0))

ggsave(popout_histogram_plot_log, filename = paste0(outdir, "popout_histogram_plot_log.png"))

combined_plot_log <- histogram_plot_log + annotation_custom(ggplotGrob(popout_histogram_plot_log), xmin = -35, xmax = -15, 
                       ymin = 500, ymax = 1225)

ggsave(combined_plot_log, filename = paste0(outdir, "histogram_plot_log_w_popout.png"), width = 5, height = 3)
ggsave(combined_plot_log, filename = paste0(outdir, "histogram_plot_log_w_popout.pdf"), width = 5, height = 3)



### Calculate number of eQTLs and write in notion for writing paper time ###
