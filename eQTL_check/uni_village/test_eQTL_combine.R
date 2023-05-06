library(data.table)
library(dplyr)
library(ggplot2)
library(colorspace)
library(ggpp)
library(ggbreak) 


indir_uni <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/uni_village/uniculture/gene_separate/beds/"
indir_village <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/uni_village/village/gene_separate/beds/"
outdir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/uni_village/combined/"


dir.create(outdir, recursive = TRUE)


### Read in results
uni_files <- list.files(indir_uni)

uni_results_list <- lapply(uni_files, function(x){
    fread(paste0(indir_uni, x), sep = "\t")
})


uni_results_dt <- do.call(rbind, uni_results_list)
uni_results_dt$fdr <- p.adjust(uni_results_dt$dataset_p, method = "fdr")
uni_results_dt$group <- "Uni-culture"


### Kilpinen
village_files <- list.files(indir_village)

village_results_list <- lapply(village_files, function(x){
    fread(paste0(indir_village, x), sep = "\t")

})


village_results_dt <- do.call(rbind, village_results_list)
village_results_dt$fdr <- p.adjust(village_results_dt$dataset_p, method = "fdr")
village_results_dt$group <- "Village"


uni_results_dt_prep <- uni_results_dt[,c("chrom", "start", "end", "maf_deboever", "ref_deboever", "alt_deboever", "gene_id_deboever", "beta_deboever", "ID_ref_alt", "REF", "ALT", "TOB0421", "FSA0006", "MBE1006", "dataset_beta", "dataset_p", "direction", "fdr", "group")]
colnames(uni_results_dt_prep) <- c("chrom", "start", "end", "maf", "REF", "ALT", "ensg", "beta", "ID_ref_alt", "dataset_REF", "dataset_ALT", "TOB0421", "FSA0006", "MBE1006", "dataset_beta", "dataset_p", "direction", "fdr", "group")
village_results_dt_prep <- village_results_dt[,c("chrom", "start", "end", "maf_deboever", "ref_deboever", "alt_deboever", "gene_id_deboever", "beta_deboever", "ID_ref_alt", "REF", "ALT", "TOB0421", "FSA0006", "MBE1006", "dataset_beta", "dataset_p", "direction", "fdr", "group")]
colnames(village_results_dt_prep) <- c("chrom", "start", "end", "maf", "REF", "ALT", "ensg", "beta", "ID_ref_alt", "dataset_REF", "dataset_ALT", "TOB0421", "FSA0006", "MBE1006", "dataset_beta", "dataset_p", "direction", "fdr", "group")


##### find common genes and snps #####
uni_results_dt_prep_gene_snp <- uni_results_dt_prep[!direction == "different_snp",c("ensg", "ID_ref_alt")]

village_results_dt_prep_gene_snp <- village_results_dt_prep[!direction == "different_snp",c("ensg", "ID_ref_alt")]

common_gene_snp <- unique(uni_results_dt_prep_gene_snp[village_results_dt_prep_gene_snp, on = c("ensg", "ID_ref_alt"), nomatch = NULL])

### Check inner join worked ###
uni_has <- c()
village_has <- c()

for (row in 1:nrow(common_gene_snp)){
    uni_has <- c(uni_has, ifelse(nrow(uni_results_dt_prep_gene_snp[ensg == common_gene_snp$ensg[row] & ID_ref_alt == common_gene_snp$ID_ref_alt[row]] > 0), TRUE, FALSE))
    village_has <- c(village_has, ifelse(nrow(uni_results_dt_prep_gene_snp[ensg == common_gene_snp$ensg[row] & ID_ref_alt == common_gene_snp$ID_ref_alt[row]] > 0), TRUE, FALSE))
}

all(uni_has)
all(village_has)

uni_results_dt_sub <- unique(uni_results_dt_prep[common_gene_snp, on = c("ensg", "ID_ref_alt")])
village_results_dt_sub <- unique(village_results_dt_prep[common_gene_snp, on = c("ensg", "ID_ref_alt")])


results_dt <- rbind(uni_results_dt_sub, village_results_dt_sub)
results_dt$direction <- factor(results_dt$direction, levels = c("opposite","match"))


results_dt$updated_beta <- ifelse(results_dt$REF != results_dt$dataset_REF, results_dt$beta * -1, results_dt$beta)

fwrite(results_dt, paste0(outdir, "eQTL_results.tsv"), sep = "\t")


# correlation <- list()
# correlation_dt <- data.table(study = unique(results_dt$study), rho = as.numeric(NA), p = as.numeric(NA))

# for (Study in results_dt$study){
#     correlation[[Study]] <- cor.test(results_dt[study == Study]$updated_beta, results_dt[study == Study]$dataset_beta, method = "spearman", exact = TRUE)
#     correlation_dt[study == Study, "rho"] <- correlation[[Study]]$estimate[1]
#     correlation_dt[study == Study, "p"] <- correlation[[Study]]$p.value
# }



correlation_plot <- ggplot(results_dt, aes(updated_beta, dataset_beta, color = direction)) +
        geom_hline(yintercept=0) +
        geom_vline(xintercept=0) +
        geom_point(alpha = 0.25) +
        theme_classic() +
        ylab("This Dataset Effect Size") +
        xlab("DeBoever Effect Size") +
        scale_color_manual(values = c("grey60", "black"), name = "Effect\nDirections") +
        stat_quadrant_counts() +
        xlim(-2.25,2.25) +
        # ylim(-8.75,8.75) +
        # ylim(-2.25,2.25) +
        theme(legend.position = "none") +
        facet_wrap(vars(group), nrow = 1) +
        scale_y_break(c(-2.9, -2), scales = 4) 
                # scale_y_break(c(2, 2.9), scale=1)

ggsave(correlation_plot, filename = paste0(outdir, "quadrant_plot_facet.png"), height = 2.5, width = 4.2)
ggsave(correlation_plot, filename = paste0(outdir, "quadrant_plot_facet.pdf"), height = 2.5, width = 4.2)


correlation_plot_noN <- ggplot(results_dt, aes(updated_beta, dataset_beta, color = direction)) +
        geom_hline(yintercept=0) +
        geom_vline(xintercept=0) +
        geom_point(size = 0.5) +
        theme_classic() +
        ylab("This Dataset Effect Size") +
        xlab("DeBoever Effect Size") +
        scale_color_manual(values = c("grey60", "black"), name = "Effect\nDirections") +
        xlim(-2.25,2.25) +
        ylim(-8,6) +
        theme(legend.position = "none") +
        facet_wrap(vars(group), nrow = 1) +
        scale_y_break(c(-5.5, -2.1), scales = 4) +
                scale_y_break(c(2, 2.85), scale=1)

ggsave(correlation_plot_noN, filename = paste0(outdir, "quadrant_plot_facet_noN.png"), height = 2.5, width = 4.2)
ggsave(correlation_plot_noN, filename = paste0(outdir, "quadrant_plot_facet_noN.pdf"), height = 2.5, width = 4.2)



### X^2 test of significance ###
xisq_df_uni <- data.frame("DeBoever" = ifelse((results_dt[group == "Uni-culture"]$dataset_beta < 0), "negative", "positive"), "Dataset" = ifelse((results_dt[group == "Uni-culture"]$updated_beta < 0), "negative", "positive"))

xisq_results_uni <- chisq.test(xisq_df_uni$DeBoever, xisq_df_uni$Dataset, correct=FALSE)
xisq_results_uni$p.value



xisq_df_vil <- data.frame("DeBoever" = ifelse((results_dt[group == "Village"]$dataset_beta < 0), "negative", "positive"), "Dataset" = ifelse((results_dt[group == "Village"]$updated_beta < 0), "negative", "positive"))

xisq_results_vil <- chisq.test(xisq_df_vil$DeBoever, xisq_df_vil$Dataset, correct=FALSE)
xisq_results_vil$p.value





##### Correlate Uni-culture to village #####
uni_results_dt_sub4wide <- uni_results_dt_sub[,c("chrom", "start", "end", "ID_ref_alt", "REF", "ALT", "ensg", "TOB0421", "FSA0006", "MBE1006", "dataset_beta", "dataset_p", "direction", "fdr")]
colnames(uni_results_dt_sub4wide) <- c("chrom", "start", "end", "ID_ref_alt", "REF", "ALT", "ensg", "TOB0421", "FSA0006", "MBE1006", "Uni-culture_beta", "Uni-culture_p", "Uni-culture_direction", "Uni-culture_fdr")

village_results_dt_sub4wide <- village_results_dt_sub[,c("chrom", "start", "end", "ID_ref_alt", "REF", "ALT", "ensg", "TOB0421", "FSA0006", "MBE1006", "dataset_beta", "dataset_p", "direction", "fdr")]
colnames(village_results_dt_sub4wide) <- c("chrom", "start", "end", "ID_ref_alt", "REF", "ALT", "ensg", "TOB0421", "FSA0006", "MBE1006", "village_beta", "village_p", "village_direction", "village_fdr")



results_wide_dt <- uni_results_dt_sub4wide[village_results_dt_sub4wide, on = c("chrom", "start", "end", "ID_ref_alt", "REF", "ALT", "ensg", "TOB0421", "FSA0006", "MBE1006")]


cor_spear <- cor.test(results_wide_dt$`Uni-culture_beta`, results_wide_dt$village_beta, method = "spearman") ### rho = 0.8828232 and p = 0


quadrant_uni_village <- ggplot(results_wide_dt, aes(`Uni-culture_beta`, village_beta)) +
        geom_hline(yintercept=0) +
        geom_vline(xintercept=0) +
        geom_point(size = 0.5) +
        theme_classic() +
        xlab("Uni-culture Effect Size") +
        ylab("Village Effect Size") +
        stat_quadrant_counts() +
        xlim(-8.5,8.5) +
        ylim(-8.5,8.5) +
        theme(legend.position = "none")

ggsave(quadrant_uni_village, filename = paste0(outdir, "quadrant_uni_village.png"), height = 2.5, width = 2.5)
ggsave(quadrant_uni_village, filename = paste0(outdir, "quadrant_uni_village.pdf"), height = 2.5, width = 2.5)




##### Make a plot of the difference in effect size between village and uni-culture #####
results_wide_dt$abs_difference <- abs(results_wide_dt$village_beta) - abs(results_wide_dt$Uni-culture_beta)


effect_diff_hist <- ggplot(results_wide_dt, aes(abs_difference)) +
    geom_histogram(bins = 100) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme_classic()

ggsave(effect_diff_hist, filename = paste0(outdir, "effect_size_diff_histogram.png"))





