library(data.table)
library(dplyr)
library(ggplot2)
library(colorspace)


indir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/kilpinen/gene_separate/beds/"
outdir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/kilpinen/combined/"

dir.create(outdir, recursive = TRUE)


### Read in results
files <- list.files(indir)

bed_results_list <- lapply(files, function(x){
    fread(paste0(indir, x), sep = "\t")
})


results_dt <- do.call(rbind, bed_results_list)


results_dt$fdr <- p.adjust(results_dt$dataset_p, method = "fdr")


results_dt <- results_dt[order(fdr)]

head(results_dt, n =100)

results_dt[TOB0421 != FSA0006 & TOB0421 != MBE1006 & FSA0006 != MBE1006]

fwrite(results_dt, paste0(outdir, "kilpinen_overlap_results.tsv"), sep = "\t")

results_dt <- results_dt[direction != "different_snp"]
results_dt$significant <- ifelse(results_dt$fdr < 0.05, TRUE, FALSE)

results_dt[Gene_ID == "CHCHD2"]

results_dt$direction <- factor(results_dt$direction, levels = c("opposite","match"))

results_dt$lmm_peer_beta_kilpinen4plotting <- ifelse(results_dt$maf_kilpinen > 0.5, (-1 * results_dt$lmm_peer_beta_kilpinen), results_dt$lmm_peer_beta_kilpinen)
results_dt$dataset_beta4plotting <- ifelse(results_dt$dataset_maf > 0.5, (-1 * results_dt$dataset_beta), results_dt$dataset_beta)



correlation_plot <- ggplot(results_dt, aes(lmm_peer_beta_kilpinen4plotting, dataset_beta4plotting, color = log(fdr))) +
        geom_point() +
        theme_classic() +
         scale_color_continuous_sequential(palette = "Purp", rev = FALSE)
        # scale_color_manual(values = c("grey","black"))

ggsave(correlation_plot, filename = paste0(outdir, "kilpinen_correlation_plot.png"))


histogram_plot <- ggplot(results_dt, aes(fdr, fill = direction)) +
        geom_histogram(position = "stack", binwidth = 0.05) +
        theme_classic() +
        geom_vline(aes(xintercept=0.05),
            color="red", linetype="dashed", size=1) +
        scale_fill_manual(values = c("grey","black"))
        
ggsave(histogram_plot, filename = paste0(outdir, "kilpinen_histogram_plot.png"))


histogram_plot_log <- ggplot(results_dt, aes(log10(fdr), fill = direction)) +
        geom_histogram(position = "stack", binwidth = 1) +
        theme_classic() +
        geom_vline(aes(xintercept=log10(0.05)),
            color="red", linetype="dashed", size=0.25) +
        scale_fill_manual(values = c("grey","black")) +
        scale_y_continuous(expand = c(0, 0)) +
        ylab("Number of Associations") +
        labs(fill = "eQTL\nEffect\nDirection\nAgreement")

        
ggsave(histogram_plot_log, filename = paste0(outdir, "kilpinen_histogram_plot_log.png"))


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

ggsave(popout_histogram_plot_log, filename = paste0(outdir, "kilpinen_popout_histogram_plot_log.png"))

combined_plot_log <- histogram_plot_log + annotation_custom(ggplotGrob(popout_histogram_plot_log), xmin = -29, xmax = -12, 
                       ymin = 275, ymax = 750)

ggsave(combined_plot_log, filename = paste0(outdir, "kilpinen_histogram_plot_log_w_popout.png"), width = 5, height = 3)
ggsave(combined_plot_log, filename = paste0(outdir, "kilpinen_histogram_plot_log_w_popout.pdf"), width = 5, height = 3)



### Calculate number of eQTLs and write in notion for writing paper time ###
nrow(results_dt)

table(results_dt$direction)

nrow(results_dt[fdr < 0.05])

table(results_dt[fdr < 0.05]$direction)

