library(data.table)
library(tidyverse)
library(ggplot2)
library(Seurat)


##### Set up directories #####
indir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/preQC/"
outdir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/growth_rate/single_cell_comparison/"
rate_dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/growth_rate/18_line_village/ratrack/more_growth_rate_points/results/"

dir.create(outdir, recursive = TRUE)



##### Read in proportions #####
diff_props <- fread(paste0(indir, "cardiac_diff_prop_lines.tsv"), sep = "\t")
multi_passage_props <- fread(paste0(indir, "multi_passage_prop_lines.tsv"), sep = "\t")

colors <- readRDS("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Nona_multiome/line_colors")



##### Calculate the ratio of proportions to baseline "proportion #####
diff_ratios_dt <- unique(diff_props[,c("Assignment", "Day")])
diff_ratios_dt$ratio <- as.numeric(NA)

multi_passage_ratios_dt <-  unique(multi_passage_props[,c("Assignment", "Day")])
multi_passage_ratios_dt$ratio <- as.numeric(NA)

for (sample in unique(diff_props$Assignment)){
    for (day in unique(diff_props$Day)){
        diff_ratios_dt[Assignment == sample & Day == day]$ratio <- diff_props[Assignment == sample & Day == day]$N/diff_props[Assignment == sample & Day == 0]$N
    }
    for (day in unique(multi_passage_ratios_dt$Day)){
        multi_passage_ratios_dt[Assignment == sample & Day == day]$ratio <- multi_passage_props[Assignment == sample & Day == day]$N/diff_props[Assignment == sample & Day == 0]$N
    }
}

multi_passage_ratios_dt <- rbind(diff_ratios_dt[Day ==0], multi_passage_ratios_dt)


pRatio_diff <- ggplot(diff_ratios_dt, aes(Day, ratio, color = Assignment)) +
    geom_point() +
    geom_line() +
    theme_classic() +
    scale_color_manual(values = colors) +
    xlab("Days") +
    ylab("Ratio of Proportions\nto Day 0")

ggsave(pRatio_diff, filename = paste0(outdir, "cardiac_diff_proportion_ratios.png"))



pRatio_multi_passage <- ggplot(multi_passage_ratios_dt, aes(Day, ratio, color = Assignment)) +
    geom_point() +
    geom_line() +
    theme_classic() +
    scale_color_manual(values = colors) +
    xlab("Passage") +
    ylab("Ratio of Proportions\nto Day 0")

ggsave(pRatio_multi_passage, filename = paste0(outdir, "multi_passage_proportion_ratios.png"))





##### Read in results #####
rate_results <- lapply(unique(diff_props$Assignment), function(sample){
    fread(paste0(rate_dir, sample, "_ratrack.fit.csv"))
})
names(rate_results) <- unique(diff_props$Assignment)


rate_results_2 <- lapply(rate_results, function(x){
    x[model_index == 1]
})


rate_results_3 <- lapply(rate_results, function(x){
    x[model_index == 2]
})





##### Read in results #####
rate_dir_2 <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/growth_rate/18_line_village/ratrack/results/"

rate_results.2 <- lapply(unique(diff_props$Assignment), function(sample){
    fread(paste0(rate_dir_2, sample, "_ratrack.fit.csv"))
})
names(rate_results.2) <- unique(diff_props$Assignment)


rate_results_2.2 <- lapply(rate_results.2, function(x){
    x[model_index == 1]
})


rate_results_3.2 <- lapply(rate_results.2, function(x){
    x[model_index == 2]
})



for (x in names(rate_results_2)){
    print(x)
    # print(dim(rate_results_2[[x]]))
    print(dim(rate_results_2.2[[x]]))
    # print(dim(rate_results_3[[x]]))
    # print(dim(rate_results_3.2[[x]])) ### This one worked for all of them so use this
}


##### Calculate the rate in the middle of the growth curve #####
rate_results_3.2 <- lapply(rate_results_3.2, function(x){
    x$name <- gsub(" ", "", x$name)
    return(x)
})

rate_results_3.2_middle <- lapply(rate_results_3.2, function(x){
    x[rate_position == 1]
})


rate_results_3.2_middle_dt <- do.call(rbind, rate_results_3.2_middle)

rate_results_3.2_middle_dt[order(rate_mean)]


rate_results_3.2_middle_dt$Assignment <- gsub("_rep[1-3]", "", rate_results_3.2_middle_dt$name)


combined_dt <- multi_passage_ratios_dt[rate_results_3.2_middle_dt, on = "Assignment", allow.cartesian=TRUE]

for (day in unique(combined_dt$Day)){
    print(day)
    combined_dt[Day == day, c("ratio_rank")] <- rank(combined_dt[Day == day]$ratio)
    combined_dt[Day == day, c("rate_rank")] <- rank(combined_dt[Day == day]$rate_mean)
    print(cor.test(combined_dt[Day == day]$ratio, combined_dt[Day == day]$rate_mean, method = "spearman"))
}




p_ratio_mean <- ggplot(combined_dt[Day != 0], aes(ratio, rate_mean, color = Assignment)) +
        geom_point() +
        facet_wrap(vars(Day), nrow = 1, scales = "free_x") +
        scale_color_manual(values = colors) +
        theme_classic()

ggsave(p_ratio_mean, filename = paste0(outdir, "ratio_rate_comparison.png"), width = 10, height = 4)



p_ratio_mean_rank <- ggplot(combined_dt[Day != 0], aes(ratio_rank, rate_rank, color = Assignment)) +
        geom_point() +
        facet_wrap(vars(Day), nrow = 1, scales = "free_x") +
        scale_color_manual(values = colors) +
        theme_classic()

ggsave(p_ratio_mean_rank, filename = paste0(outdir, "ratio_rate_comparison_rank.png"), width = 10, height = 4)








rate_results_2.2 <- lapply(rate_results_2.2, function(x){
    x$name <- gsub(" ", "", x$name)
    return(x)
})

rate_results_2.2_middle <- lapply(rate_results_2.2, function(x){
    tmp <- data.table(name = unique(x$name), Assignment = gsub("_rep[1-3]","",unique(x$name)), rate = as.numeric(NA))
    for (line in tmp$name){
        tmp[name == line, c("rate")] <- (x[rate_position == 1 & name == line]$rate_mean + x[rate_position == 0 & name == line]$rate_mean)/2
    }
    return(tmp)
})


multi_passage_ratios_dt_2 <- do.call(rbind, rate_results_2.2_middle)
combined_dt_2 <- multi_passage_ratios_dt[multi_passage_ratios_dt_2, on = "Assignment", allow.cartesian=TRUE]

starting_prop <- diff_props[Day == 0]
starting_prop$Day <- NULL

combined_dt_2 <- starting_prop[combined_dt_2, on = "Assignment"]


for (day in unique(combined_dt_2$Day)){
    print(day)
    combined_dt_2[Day == day, c("ratio_rank")] <- rank(combined_dt_2[Day == day]$ratio)
    combined_dt_2[Day == day, c("rate_rank")] <- rank(combined_dt_2[Day == day]$rate)
    print(cor.test(combined_dt_2[Day == day]$ratio, combined_dt_2[Day == day]$rate, method = "spearman"))
}


p_ratio_mean2 <- ggplot(combined_dt_2[Day != 0], aes(ratio, rate, color = Assignment, size = N)) +
        geom_point(alpha = 0.8) +
        facet_wrap(vars(Day), nrow = 1, scales = "free_x") +
        scale_color_manual(values = colors) +
        theme_classic()+
        geom_smooth() +
        geom_smooth(aes(ratio, rate), color = "black", span = 5)

ggsave(p_ratio_mean2, filename = paste0(outdir, "ratio_rate_comparison2.png"), width = 10, height = 4)
ggsave(p_ratio_mean2, filename = paste0(outdir, "ratio_rate_comparison2.pdf"), width = 10, height = 4)



p_ratio_mean_rank2 <- ggplot(combined_dt_2[Day != 0], aes(ratio_rank, rate_rank, color = Assignment, size = N)) +
        geom_point() +
        facet_wrap(vars(Day), nrow = 1, scales = "free_x") +
        scale_color_manual(values = colors) +
        theme_classic()+
        geom_smooth(method = "lm",aes(ratio_rank, rate_rank), color = "black")

ggsave(p_ratio_mean_rank2, filename = paste0(outdir, "ratio_rate_comparison_rank2.png"), width = 10, height = 4)







diff_combined_dt_2 <- diff_ratios_dt[multi_passage_ratios_dt_2, on = "Assignment", allow.cartesian=TRUE]
diff_combined_dt_2 <- starting_prop[diff_combined_dt_2, on = "Assignment"]


for (day in unique(diff_combined_dt_2$Day)){
    print(day)
    diff_combined_dt_2[Day == day, c("ratio_rank")] <- rank(diff_combined_dt_2[Day == day]$ratio)
    diff_combined_dt_2[Day == day, c("rate_rank")] <- rank(diff_combined_dt_2[Day == day]$rate)
    print(cor.test(diff_combined_dt_2[Day == day]$ratio, diff_combined_dt_2[Day == day]$rate, method = "spearman"))
}


p_diff_ratio_mean2 <- ggplot(diff_combined_dt_2[Day != 0], aes(ratio, rate, color = Assignment, size = N)) +
        geom_point(alpha = 0.8) +
        facet_wrap(vars(paste0("Day ", Day)), nrow = 1, scales = "free_x") +
        scale_color_manual(values = colors) +
        theme_classic() +
        geom_smooth() +
        geom_smooth(aes(ratio, rate), color = "black", span = 5) +
        xlab("Proportion of Village Relative to Day 0 Proportion") +
        ylab("Growth Rate")

ggsave(p_diff_ratio_mean2, filename = paste0(outdir, "dif_ratio_rate_comparison2.png"), width = 15, height = 3)
ggsave(p_diff_ratio_mean2, filename = paste0(outdir, "dif_ratio_rate_comparison2.pdf"), width = 15, height = 3)



p_diff_ratio_mean_rank2 <- ggplot(diff_combined_dt_2[Day != 0], aes(ratio_rank, rate_rank, color = Assignment, size = N)) +
        geom_point() +
        facet_wrap(vars(Day), nrow = 1, scales = "free_x") +
        scale_color_manual(values = colors) +
        theme_classic() +
        geom_smooth(method = "lm",aes(ratio_rank, rate_rank), color = "black")

ggsave(p_diff_ratio_mean_rank2, filename = paste0(outdir, "dif_ratio_rate_comparison_rank2.png"), width = 10, height = 3)






##### Just proportion of the village at that time #####
multi_passage_proportion_combined_dt <- multi_passage_props[multi_passage_ratios_dt_2, on = "Assignment", allow.cartesian=TRUE]
colnames(multi_passage_proportion_combined_dt) <- gsub("N", "Proportion", colnames(multi_passage_proportion_combined_dt))
multi_passage_proportion_combined_dt <- starting_prop[multi_passage_proportion_combined_dt, on = "Assignment"]


for (day in unique(multi_passage_proportion_combined_dt$Day)){
    print(day)
    multi_passage_proportion_combined_dt[Day == day, c("ratio_rank")] <- rank(multi_passage_proportion_combined_dt[Day == day]$Proportion)
    multi_passage_proportion_combined_dt[Day == day, c("rate_rank")] <- rank(multi_passage_proportion_combined_dt[Day == day]$rate)
    print(cor.test(multi_passage_proportion_combined_dt[Day == day]$Proportion, multi_passage_proportion_combined_dt[Day == day]$rate, method = "spearman"))
}


multi_passage_proportion_cor <- data.table(Day = unique(multi_passage_proportion_combined_dt$Day), rho = as.numeric(NA), P = as.numeric(NA))

for (day in unique(multi_passage_proportion_combined_dt$Day)){
    print(day)
    multi_passage_proportion_combined_dt[Day == day, c("ratio_rank")] <- rank(multi_passage_proportion_combined_dt[Day == day]$Proportion)
    multi_passage_proportion_combined_dt[Day == day, c("rate_rank")] <- rank(multi_passage_proportion_combined_dt[Day == day]$rate)
    test <- cor.test(multi_passage_proportion_combined_dt[Day == day]$Proportion, multi_passage_proportion_combined_dt[Day == day]$rate, method = "spearman")
    multi_passage_proportion_cor[Day == day, c("rho")] <- test$estimate
    multi_passage_proportion_cor[Day == day, c("P")] <- test$p.value
}

multi_passage_proportion_cor$FDR <- p.adjust (multi_passage_proportion_cor$P, method="BH")



p_multi_passage_prop_mean2 <- ggplot(multi_passage_proportion_combined_dt[Day != 0], aes(Proportion, rate, color = Assignment, size = N)) +
        geom_point(alpha = 0.8) +
        facet_wrap(vars(paste0("Passage ", Day)), nrow = 1, scales = "free_x") +
        scale_color_manual(values = colors) +
        theme_classic() +
        geom_smooth() +
        geom_smooth(aes(Proportion, rate), color = "black", span = 5) +
        xlab("Proportion of Village") +
        ylab("Growth Rate") +
        theme(axis.text.x = element_text(size = 12, angle = 90, hjust = .5, vjust = .5),
                axis.text.y = element_text(size = 12),
                axis.title.x = element_text(size = 16),
                axis.title.y = element_text(size = 16),
                strip.text.x = element_text(size = 20))

ggsave(p_multi_passage_prop_mean2, filename = paste0(outdir, "multi_passage_prop_rate_comparison2.png"), width =8, height = 3)
ggsave(p_multi_passage_prop_mean2, filename = paste0(outdir, "multi_passage_prop_rate_comparison2.pdf"), width =8, height = 3)
ggsave(p_multi_passage_prop_mean2, filename = paste0(outdir, "multi_passage_prop_rate_comparison2_tall.pdf"), width =8, height = 10)



p_multi_passage_prop_mean_rank2 <- ggplot(multi_passage_proportion_combined_dt[Day != 0], aes(ratio_rank, rate_rank, color = Assignment, size = N)) +
        geom_point() +
        facet_wrap(vars(Day), nrow = 1, scales = "free_x") +
        scale_color_manual(values = colors) +
        theme_classic() +
        geom_smooth(method = "lm",aes(ratio_rank, rate_rank), color = "black")

ggsave(p_multi_passage_prop_mean_rank2, filename = paste0(outdir, "multi_passage_prop_rate_comparison_rank2.png"), width = 10, height = 3)


fwrite(multi_passage_proportion_combined_dt, paste0(outdir,"multi-passage_growth_rate_prop.tsv"), sep = "\t")



diff_proportion_combined_dt <- diff_props[multi_passage_ratios_dt_2, on = "Assignment", allow.cartesian=TRUE]
colnames(diff_proportion_combined_dt) <- gsub("N", "Proportion", colnames(diff_proportion_combined_dt))
diff_proportion_combined_dt <- starting_prop[diff_proportion_combined_dt, on = "Assignment"]


diff_proportion_cor_list <- list()
diff_proportion_cor <- data.table(Day = unique(diff_proportion_combined_dt$Day), rho = as.numeric(NA), P = as.numeric(NA))

for (day in unique(diff_proportion_combined_dt$Day)){
    print(day)
    diff_proportion_combined_dt[Day == day, c("ratio_rank")] <- rank(diff_proportion_combined_dt[Day == day]$Proportion)
    diff_proportion_combined_dt[Day == day, c("rate_rank")] <- rank(diff_proportion_combined_dt[Day == day]$rate)
    test <- cor.test(diff_proportion_combined_dt[Day == day]$Proportion, diff_proportion_combined_dt[Day == day]$rate, method = "spearman")
    diff_proportion_cor[Day == day, c("rho")] <- test$estimate
    diff_proportion_cor[Day == day, c("P")] <- test$p.value
}

diff_proportion_cor$FDR <- p.adjust (diff_proportion_cor$P, method="BH")



p_diff_prop_mean2 <- ggplot(diff_proportion_combined_dt[Day != 0], aes(Proportion, rate, color = Assignment, size = N)) +
        geom_point(alpha = 0.8) +
        facet_wrap(vars(paste0("Day ",Day)), nrow = 1, scales = "free_x") +
        scale_color_manual(values = colors) +
        theme_classic() +
        geom_smooth() +
        geom_smooth(aes(Proportion, rate), color = "black", span = 5) +
        xlab("Proportion of Village") +
        ylab("Growth Rate") +
        theme(axis.text.x = element_text(size = 12, angle = 90, hjust = .5, vjust = .5),
                axis.text.y = element_text(size = 12),
                axis.title.x = element_text(size = 16),
                axis.title.y = element_text(size = 16),
                strip.text.x = element_text(size = 20))

ggsave(p_diff_prop_mean2, filename = paste0(outdir, "dif_prop_rate_comparison2.png"), width = 15, height = 3)
ggsave(p_diff_prop_mean2, filename = paste0(outdir, "dif_prop_rate_comparison2.pdf"), width = 15, height = 3)
ggsave(p_diff_prop_mean2, filename = paste0(outdir, "dif_prop_rate_comparison2_tall.pdf"), width = 15, height = 10)




p_diff_prop_mean2_d0 <- ggplot(diff_proportion_combined_dt, aes(Proportion, rate, color = Assignment, size = N)) +
        geom_point(alpha = 0.8) +
        facet_wrap(vars(factor(paste0("Day ",Day), levels = paste0("Day ",c(0:5,7,15)))), nrow = 1, scales = "free_x") +
        scale_color_manual(values = colors) +
        theme_classic() +
        geom_smooth() +
        geom_smooth(aes(Proportion, rate), color = "black", span = 5) +
        xlab("Proportion of Village") +
        ylab("Growth Rate") +
        theme(axis.text.x = element_text(size = 12, angle = 90, hjust = .5, vjust = .5),
                axis.text.y = element_text(size = 12),
                axis.title.x = element_text(size = 16),
                axis.title.y = element_text(size = 16),
                strip.text.x = element_text(size = 20))

ggsave(p_diff_prop_mean2_d0, filename = paste0(outdir, "dif_prop_rate_comparison2_w_d0.png"), width = 15, height = 3)
ggsave(p_diff_prop_mean2_d0, filename = paste0(outdir, "dif_prop_rate_comparison2_w_d0.pdf"), width = 15, height = 3)
ggsave(p_diff_prop_mean2_d0, filename = paste0(outdir, "dif_prop_rate_comparison2_w_d0_tall.pdf"), width = 15, height = 10)




p_diff_prop_mean_rank2 <- ggplot(diff_proportion_combined_dt[Day != 0], aes(ratio_rank, rate_rank, color = Assignment, size = N)) +
        geom_point() +
        facet_wrap(vars(Day), nrow = 1, scales = "free_x") +
        scale_color_manual(values = colors) +
        theme_classic() +
        geom_smooth(method = "lm",aes(ratio_rank, rate_rank), color = "black")

ggsave(p_diff_prop_mean_rank2, filename = paste0(outdir, "dif_prop_rate_comparison_rank2.png"), width = 10, height = 3)



fwrite(diff_proportion_combined_dt, paste0(outdir, "18-line_growth_prop.tsv"), sep = "\t")