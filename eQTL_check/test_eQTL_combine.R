library(data.table)
library(dplyr)
library(ggplot2)
library(colorspace)
library(ggpp)


indir_deboever <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/deboever/gene_separate/beds/"
indir_kilpinen <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/kilpinen/gene_separate/beds/"
outdir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/eQTL_overlap/combined/"


dir.create(outdir, recursive = TRUE)


### Read in results
deboever_files <- list.files(indir_deboever)

bed_results_list <- lapply(deboever_files, function(x){
    fread(paste0(indir_deboever, x), sep = "\t")
})


deboever_results_dt <- do.call(rbind, bed_results_list)
deboever_results_dt$fdr <- p.adjust(deboever_results_dt$dataset_p, method = "fdr")
deboever_results_dt$study <- "DeBoever"


### Kilpinen
kilpinen_files <- list.files(indir_kilpinen)

kilpinen_results_list <- lapply(kilpinen_files, function(x){
    fread(paste0(indir_kilpinen, x), sep = "\t")

})


kilpinen_results_dt <- do.call(rbind, kilpinen_results_list)
kilpinen_results_dt$fdr <- p.adjust(kilpinen_results_dt$dataset_p, method = "fdr")
kilpinen_results_dt$study <- "Kilpinen"


deboever_results_dt_prep <- deboever_results_dt[,c("chrom", "start", "end", "maf_deboever", "ref_deboever", "alt_deboever", "gene_id_deboever", "beta_deboever", "ID_ref_alt", "REF", "ALT", "TOB0421", "FSA0006", "MBE1006", "dataset_beta", "dataset_p", "direction", "fdr", "study")]
colnames(deboever_results_dt_prep) <- c("chrom", "start", "end", "maf", "REF", "ALT", "ensg", "beta", "ID_ref_alt", "dataset_REF", "dataset_ALT", "TOB0421", "FSA0006", "MBE1006", "dataset_beta", "dataset_p", "direction", "fdr", "study")
kilpinen_results_dt_prep <- kilpinen_results_dt[,c("#chrom", "start", "end", "maf_kilpinen", "REF", "ALT", "gene_id_kilpinen", "lmm_peer_beta_kilpinen", "ID_ref_alt", "REF", "ALT", "TOB0421", "FSA0006", "MBE1006", "dataset_beta",  "dataset_p", "direction", "fdr", "study")]
colnames(kilpinen_results_dt_prep) <- c("chrom", "start", "end", "maf", "REF", "ALT", "ensg", "beta", "ID_ref_alt", "dataset_REF", "dataset_ALT", "TOB0421", "FSA0006", "MBE1006", "dataset_beta",  "dataset_p", "direction", "fdr", "study")
kilpinen_results_dt_prep$REF <- NA
kilpinen_results_dt_prep$ALT <- NA


results_dt <- rbind(deboever_results_dt_prep, kilpinen_results_dt_prep)


results_dt <- results_dt[order(fdr)]

head(results_dt, n =100)

results_dt[TOB0421 != FSA0006 & TOB0421 != MBE1006 & FSA0006 != MBE1006]

fwrite(results_dt, paste0(outdir, "overlap_results.tsv"), sep = "\t")

results_dt <- results_dt[direction != "different_snp"]
results_dt$significant <- ifelse(results_dt$fdr < 0.05, TRUE, FALSE)

results_dt[ensg == "ENSG00000106153"]

results_dt$direction <- factor(results_dt$direction, levels = c("opposite","match"))


results_dt$updated_beta <- ifelse(results_dt$study == "Kilpinen" & results_dt$direction == "match" & results_dt$beta * results_dt$dataset_beta < 0, results_dt$beta * -1, 
    ifelse(results_dt$study == "Kilpinen" & results_dt$direction == "opposite" & results_dt$beta * results_dt$dataset_beta > 0, results_dt$beta * -1, 
        ifelse(results_dt$study == "DeBoever" & results_dt$REF != results_dt$dataset_REF, results_dt$beta * -1, results_dt$beta)))


correlation <- list()
correlation_dt <- data.table(study = unique(results_dt$study), rho = as.numeric(NA), p = as.numeric(NA))

for (Study in results_dt$study){
    correlation[[Study]] <- cor.test(results_dt[study == Study]$updated_beta, results_dt[study == Study]$dataset_beta, method = "spearman", exact = TRUE)
    correlation_dt[study == Study, "rho"] <- correlation[[Study]]$estimate[1]
    correlation_dt[study == Study, "p"] <- correlation[[Study]]$p.value
}

correlation_plot <- ggplot(results_dt, aes(updated_beta, dataset_beta, color = significant)) +
        geom_hline(yintercept=0) +
        geom_vline(xintercept=0) +
        geom_point(alpha = 0.25) +
        theme_classic() +
        # scale_color_continuous_sequential(palette = "Purp", rev = FALSE) +
        scale_color_manual(values = c("lightgrey", "black")) +
        facet_wrap(vars(study))
        # scale_color_manual(values = c("grey","black"))

ggsave(correlation_plot, filename = paste0(outdir, "correlation_plot.png"))



correlation_plot <- ggplot(results_dt[study == "DeBoever" & significant == TRUE], aes(updated_beta, dataset_beta, color = direction)) +
        geom_hline(yintercept=0) +
        geom_vline(xintercept=0) +
        geom_point(alpha = 0.25) +
        theme_classic() +
        ylab("This Dataset Effect Size") +
        xlab("DeBoever Effect Size") +
        scale_color_manual(values = c("grey60", "black"), name = "Effect\nDirections") +
        stat_quadrant_counts() +
        xlim(-2.25,2.25) +
        ylim(-8.5,8.5) +
        theme(legend.position = "none")

ggsave(correlation_plot, filename = paste0(outdir, "correlation_plot_DeBoever.png"), height = 2.5, width = 2.5)
ggsave(correlation_plot, filename = paste0(outdir, "correlation_plot_DeBoever.pdf"), height = 2.5, width = 2.5)




### add number of genes in each quadrant
dim(results_dt[updated_beta < 0 & dataset_beta < 0]) # 2255
dim(results_dt[updated_beta > 0 & dataset_beta > 0]) # 2116

dim(results_dt[updated_beta < 0 & dataset_beta > 0]) # 1071
dim(results_dt[updated_beta > 0 & dataset_beta < 0]) # 1008


dim(results_dt[updated_beta < 0 & dataset_beta < 0 & study == "DeBoever"]) # 1534
dim(results_dt[updated_beta > 0 & dataset_beta > 0 & study == "DeBoever"]) # 1378

dim(results_dt[updated_beta < 0 & dataset_beta > 0 & study == "DeBoever"]) # 688
dim(results_dt[updated_beta > 0 & dataset_beta < 0 & study == "DeBoever"]) # 600


dim(results_dt[updated_beta < 0 & dataset_beta < 0 & study == "Kilpinen"]) # 721
dim(results_dt[updated_beta > 0 & dataset_beta > 0 & study == "Kilpinen"]) # 738

dim(results_dt[updated_beta < 0 & dataset_beta > 0 & study == "Kilpinen"]) # 383
dim(results_dt[updated_beta > 0 & dataset_beta < 0 & study == "Kilpinen"]) # 408







histogram_plot <- ggplot(results_dt, aes(fdr, fill = direction)) +
        geom_histogram(position = "stack", binwidth = 0.05) +
        theme_classic() +
        geom_vline(aes(xintercept=0.05),
            color="red", linetype="dashed", size=1) +
        scale_fill_manual(values = c("grey","black")) +
        facet_wrap(vars(study))
        
ggsave(histogram_plot, filename = paste0(outdir, "histogram_plot.png"))

histogram_plot_log <- ggplot(results_dt, aes(log10(fdr), fill = direction)) +
        geom_histogram(position = "stack", binwidth = 1) +
        theme_classic() +
        geom_vline(aes(xintercept=log10(0.05)),
            color="red", linetype="dashed", size=0.25) +
        scale_fill_manual(values = c("grey","black")) +
        scale_y_continuous(expand = c(0, 0)) +
        ylab("Number of Associations") +
        labs(fill = "eQTL\nEffect\nDirection\nAgreement") +
        facet_wrap(vars(study), scales = "free")

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
        scale_y_continuous(expand = c(0, 0)) +
        facet_wrap(vars(study), scales = "free")


ggsave(popout_histogram_plot_log, filename = paste0(outdir, "popout_histogram_plot_log.png"))

annotation_custom2 <- function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data) 
{
  layer(data = data, stat = StatIdentity, position = PositionIdentity, 
        geom = ggplot2:::GeomCustomAnn,
        inherit.aes = TRUE, params = list(grob = grob, 
                                          xmin = xmin, xmax = xmax, 
                                          ymin = ymin, ymax = ymax))
}

get_inset <- function(df){
  p <- ggplot(data=df %>% 
                group_by(study) %>% 
                slice(1),
               aes(x = fdr, fill = direction)) +
    geom_histogram(position = "stack", binwidth = 1) +
        theme_classic() +
        scale_fill_manual(values = c("grey","black")) +
        theme(axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            legend.position="none",
            axis.text.x=element_text(size=7),
            axis.text.y=element_text(size=7)) +
        scale_y_continuous(expand = c(0, 0)) 
  return(p)
}


insets <- results_dt %>% 
  split(f = .$study) %>%
  purrr::map(~annotation_custom2(
    grob = ggplotGrob(get_inset(.)), 
    data = data.frame(category=unique(.$study)),
    ymin = -500, ymax=1225, xmin=-32, xmax=-15)
  )


combined_plot_log <- histogram_plot_log + annotation_custom2(ggplotGrob(popout_histogram_plot_log), xmin = -35, xmax = -15, 
                       ymin = 500, ymax = 1225)

combined_plot_log <- histogram_plot_log + annotation_custom(grob=ggplotGrob(insets), xmin = -35, xmax = -15, 
                       ymin = 500, ymax = 1225)


ggsave(combined_plot_log, filename = paste0(outdir, "histogram_plot_log_w_popout.png"), width = 5, height = 2.5)
ggsave(combined_plot_log, filename = paste0(outdir, "histogram_plot_log_w_popout.pdf"), width = 5, height = 2.5)



### Calculate number of eQTLs and write in notion for writing paper time ###
