library(data.table)
library(dplyr)
library(ggplot2)
library(colorspace)
library(ggpp)


indir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/deboever/gene_separate/beds/"
outdir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/deboever/combined/"

dir.create(outdir, recursive = TRUE)


### Read in results
files <- list.files(indir)

bed_results_list <- lapply(files, function(x){
    fread(paste0(indir, x), sep = "\t", header = TRUE)
})


results_dt <- do.call(rbind, bed_results_list)


results_dt$fdr <- p.adjust(results_dt$dataset_p, method = "fdr")


results_dt <- results_dt[order(fdr)]

head(results_dt, n =100)


fwrite(results_dt, paste0(outdir, "deboever_overlap_results.tsv"), sep = "\t")

results_dt <- results_dt[direction != "different_snp"]
results_dt$significant <- ifelse(results_dt$fdr < 0.05, TRUE, FALSE)


results_dt$direction <- factor(results_dt$direction, levels = c("opposite","match"))

fwrite(results_dt, paste0(outdir,"deboever_comparison_multipassage.tsv"), "\t")


correlation_plot <- ggplot(results_dt, aes(beta_deboever, dataset_beta, color = log(fdr))) +
        geom_point() +
        theme_classic() +
         scale_color_continuous_sequential(palette = "Purp", rev = FALSE)
        # scale_color_manual(values = c("grey","black"))

ggsave(correlation_plot, filename = paste0(outdir, "correlation_plot.png"))


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


results_dt$updated_beta <- ifelse(results_dt$ref_deboever != results_dt$REF, results_dt$beta_deboever * -1, results_dt$beta_deboever)




correlation_plot <- ggplot(results_dt, aes(updated_beta, dataset_beta, color = direction)) +
        geom_hline(yintercept=0) +
        geom_vline(xintercept=0) +
        geom_point(alpha = 0.5) +
        theme_classic() +
        ylab("This Dataset Effect Size") +
        xlab("DeBoever Effect Size") +
        scale_color_manual(values = c("grey60", "black"), name = "Effect\nDirections") +
        stat_quadrant_counts() +
        xlim(-2,2) +
        ylim(-2.5,2.5) +
        theme(legend.position = "none")

ggsave(correlation_plot, filename = paste0(outdir, "correlation_plot_DeBoever.png"), height = 2.5, width = 2.5)
ggsave(correlation_plot, filename = paste0(outdir, "correlation_plot_DeBoever.pdf"), height = 2.5, width = 2.5)

