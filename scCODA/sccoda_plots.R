library(data.table)
library(tidyverse)



dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/scCODA/"

line_colors <- c(FSA0006 = "#F79E29", MBE1006 = "#9B2C99", TOB0421 = "#35369C")
groups <- c("brisbane", "melbourne", "sydney", "sydney_cryo")

##### Read in Results #####
results <- lapply(groups, function(group){
    tmp <- fread(paste0(dir, group, "_logfc_results.tsv"))
    tmp$group <- group
    return(tmp)
})

results_dt <- do.call(rbind, results)
results_dt$V1 <- NULL
results_dt$group <- gsub("brisbane", "Brisbane", results_dt$group) %>%
                        gsub("melbourne", "Melbourne", .) %>%
                        gsub("^sydney$", "Sydney", .) %>%
                        gsub("sydney_cryo", "Sydney\nCryopreserved", .)


results_long <- melt(results_dt, id.vars = c("Cell Type", "is_credible", "group"),
                measure.vars = c("FSA0006_ref", "MBE1006_ref", "TOB0421_ref"))



plot <- ggplot(results_long, aes(`Cell Type`, value, color = `Cell Type`)) +
            geom_hline(yintercept = 0, linetype = "dashed") +
            geom_jitter(width = 0.2, alpha = 0.7) +
            theme_classic() +
            scale_color_manual(values = line_colors) +
            facet_wrap(vars(group), nrow = 1) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                axis.title.x = element_blank()) +
            ylab("log2(Fold Change)") +
            geom_text(
                data = unique(results_long[,c("Cell Type", "is_credible", "group")]),
                aes(y = 1.6,label = gsub(TRUE, "\\^", is_credible) %>% gsub(FALSE, " ", .)), 
                color = "black",
                # position = position_dodge(width = 1), 
                # vjust = -0.5, 
                size = 3, 
                # stat = "unique", 
                # parse = TRUE
            )

ggsave(plot, filename = paste0(dir, "fold_change_plot.png"), height = 3, width = 6)
ggsave(plot, filename = paste0(dir, "fold_change_plot.pdf"), height = 3, width = 6)




plot_fresh <- ggplot(results_long[group != "Sydney\nCryopreserved"], aes(`Cell Type`, value, color = `Cell Type`)) +
            geom_hline(yintercept = 0, linetype = "dashed") +
            geom_jitter(width = 0.2, alpha = 0.7) +
            theme_classic() +
            scale_color_manual(values = line_colors) +
            facet_wrap(vars(group), nrow = 1) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                axis.title.x = element_blank()) +
            ylab("log2(Fold Change)") +
            geom_text(
                data = unique(results_long[group != "Sydney\nCryopreserved",c("Cell Type", "is_credible", "group")]),
                aes(y = 1.6,label = gsub(TRUE, "\\^", is_credible) %>% gsub(FALSE, " ", .)), 
                color = "black",
                # position = position_dodge(width = 1), 
                # vjust = -0.5, 
                size = 3, 
                # stat = "unique", 
                # parse = TRUE
            )

ggsave(plot_fresh, filename = paste0(dir, "fold_change_plot_fresh.png"), height = 3, width = 4.5)
ggsave(plot_fresh, filename = paste0(dir, "fold_change_plot_fresh.pdf"), height = 3, width = 4.5)




plot_cryo <- ggplot(results_long[grep("Sydney", group)], aes(`Cell Type`, value, color = `Cell Type`)) +
            geom_hline(yintercept = 0, linetype = "dashed") +
            geom_jitter(width = 0.2, alpha = 0.7) +
            theme_classic() +
            scale_color_manual(values = line_colors) +
            facet_wrap(vars(group), nrow = 1) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                axis.title.x = element_blank()) +
            ylab("log2(Fold Change)") +
            geom_text(
                data = unique(results_long[grep("Sydney", group),c("Cell Type", "is_credible", "group")]),
                aes(y = 1.6,label = gsub(TRUE, "\\^", is_credible) %>% gsub(FALSE, " ", .)), 
                color = "black",
                # position = position_dodge(width = 1), 
                # vjust = -0.5, 
                size = 3, 
                # stat = "unique", 
                # parse = TRUE
            )

ggsave(plot_cryo, filename = paste0(dir, "fold_change_plot_cryo.png"), height = 3, width = 3.75)
ggsave(plot_cryo, filename = paste0(dir, "fold_change_plot_cryo.pdf"), height = 3, width = 3.75)

