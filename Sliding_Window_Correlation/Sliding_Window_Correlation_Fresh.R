library(data.table)
library(Seurat)
library(tidyverse) 
library(colorspace)


dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
datadir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/nb_partitioning_village_separate/combined/"
faultdir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/nb_partitioning_village_separate/gene_separated/"
outdir <- paste0(dir,"output/Sliding_Window_Correlation/Fresh/")
dir.create(outdir, recursive = TRUE)


##### Set up colors #####
variable_colors <- c(Village = "#A2B0D0", Replicate = "#64A66B", Line = "#68319B") 
line_colors <- c(FSA0006 = "#F79E29", MBE1006 = "#9B2C99", TOB0421 = "#35369C")

site_updates <- c("Brisbane" = "Site 1", "Melbourne" = "Site 2" ,"Sydney" = "Site 3")


##### Read in Data #####
fits_df <- fread(paste0(datadir,"nb_variance_partitioning_df.tsv"), sep = "\t")



##### Remove faulty model genes and locations #####
### Read in faulty log files ###
faulty <- fread(paste0(faultdir,"faulty_runs.txt"), header = FALSE)


### Read in the seurat data objects
seurat_list <- list()
for (location in unique(fits_df$Location)){
	for (time in unique(fits_df$Village)){
		seurat_list[[paste0(location, "_", time)]] <- readRDS(paste0("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/nb_partitioning_village_separate/data/",time, "_", location,"_seurat.rds"))
	}
}

### Make dataframe of location x number for faulty models ###
faulty_df <- data.table(Location = gsub("_nb\\.o.+", "", faulty$V1), Number = gsub(".+nb\\.o\\d+\\.", "", faulty$V1))

## Add gene id to the faulty dataframe
faulty_df$ENSG <- NA
for (row in 1:nrow(faulty_df)){
	print(row)
	faulty_df$ENSG[row] <- rownames(seurat_list[[faulty_df$Location[row]]])[as.numeric(faulty_df$Number[row])] 
}


### Remove the faulty genes ###
fits_df <- setDT(fits_df)[!faulty_df, on = c("Location", "Gene" = "ENSG")]


for (location in names(site_updates)){
	fits_df$Location <- gsub(location, site_updates[location], fits_df$Location) %>% gsub("_Cryo", " Cryo",.) %>% gsub("_Fresh", "", .)
}



##### Separate the fits_df into list of dfs - separate by Site
groups <- unique(fits_df[,c("Location", "Village")])
fits_df_list <- list()

for (location in unique(groups$Location)){
	for (time in unique(groups$Village)){
		fits_df_list[[location]][[time]] <- fits_df[Location == location & Village == time]
	}
}



### Order the genes by the line effect ###
fits_df_line_list <- lapply(fits_df_list, function(x){
	lapply(x, function(y){
		tmp <- y[!is.na(y$Line)]
		tmp <- tmp[order(tmp$Line)]
		rownames(tmp) <- tmp$Gene
		return(tmp)
	})
})

# fits_df_village_list <- lapply(fits_df_list, function(x){
# 	lapply(x, function(y){
# 		tmp <- y[order(y$Village)]
# 		rownames(tmp) <- tmp$Gene
# 		return(tmp)
# 	})
# })



combinations <- data.table(Location = rep(c("Site 1", "Site 2", "Site 3", "Site 3 Cryopreserved"), 2),
							Time1 = c(rep("Village",4),rep("Baseline", 4)),
							Time2 =c(rep("Baseline",4),rep("Village", 4)))



### Make df ###
cor_list <- list()
df_list <- list()
plots_line <- list()
spearman_plot_list <- list()


for (row in 1:nrow(combinations)){
	location <- combinations[row]$Location
	time1 <- combinations[row]$Time1
	time2 <- combinations[row]$Time2

	all_genes_line <- fits_df_line_list[[location]][[time1]]$Gene[fits_df_line_list[[location]][[time1]]$Gene %in% fits_df_line_list[[location]][[time2]]$Gene]
	df_list[[location]][[time1]][[time2]] <- data.table(Window = seq(1:(length(all_genes_line) - 100)), R = 0,  pvalue = NA, Min_Line_Variance_Explained = 0, Max_Line_Variance_Explained = 0)

	for (window in df_list[[location]][[time1]][[time2]]$Window){
		print(paste(row,window))
		genes_line <- all_genes_line[window:(window+100)]
		cor_list[[location]][[time1]][[time2]][[window]] <- cor.test(fits_df_line_list[[location]][[time1]][Gene %in% genes_line]$Line, fits_df_line_list[[location]][[time2]][,Line[match(genes_line, Gene)]], method = "spearman")
		df_list[[location]][[time1]][[time2]][Window == window, "R"] <- cor_list[[location]][[time1]][[time2]][[window]]$estimate
		if (!is.null(cor_list[[location]][[time1]][[time2]][[window]]$p.value)){
			df_list[[location]][[time1]][[time2]][Window == window, "pvalue"] <- cor_list[[location]][[time1]][[time2]][[window]]$p.value
		}
		df_list[[location]][[time1]][[time2]][Window == window, "Min_Line_Variance_Explained"] <- min(fits_df_line_list[[location]][[time1]][Gene %in% genes_line]$Line)
		df_list[[location]][[time1]][[time2]][Window == window, "Max_Line_Variance_Explained"] <- max(fits_df_line_list[[location]][[time1]][Gene %in% genes_line]$Line)

		if (window %in% c(1,14340,14419) & location == "Site 1" & time1 == "Baseline" & time2 == "Village"){
			df_tmp <- data.table(Baseline = fits_df_line_list[[location]][[time1]][Gene %in% genes_line]$Line, Village = fits_df_line_list[[location]][[time2]][,Line[match(genes_line, Gene)]])
			spearman_plot_list[[window]] <- ggplot(df_tmp, aes(Baseline, Village)) +
												geom_point() +
												theme_classic() +
												theme(axis.text.x=element_blank(),
													axis.ticks.x=element_blank(),
													axis.text.y=element_blank(),
													axis.ticks.y=element_blank())
			ggsave(spearman_plot_list[[window]], filename = gsub(" ", "_",paste0(outdir,location, "_", time1, "_",time2,"_window", window, "_correlation.png")), width = 2, height = 2)
			ggsave(spearman_plot_list[[window]], filename = gsub(" ", "_",paste0(outdir,location, "_", time1, "_",time2,"_window", window, "_correlation.pdf")), width = 2, height = 2)						
		}
	}

	plots_line[[location]][[time1]][[time2]] <- ggplot(df_list[[location]][[time1]][[time2]], aes(Min_Line_Variance_Explained*100, R, color = Max_Line_Variance_Explained*100)) +
		geom_point() +
		theme_classic() +
		scale_color_continuous_sequential(begin = 0, end = 1, palette = "SunsetDark", rev = FALSE, name = "Maximum\nPercent\nGene\nExpression\nExplained\nin\nWindow") +
		geom_segment(x=min(df_list[[location]][[time1]][[time2]][R > 0.5]$Min_Line_Variance_Explained)*100, xend=min(df_list[[location]][[time1]][[time2]][R > 0.5]$Min_Line_Variance_Explained)*100, y = -Inf, yend = 0.5,linetype = "longdash", color = "black")+
		geom_segment(x=-Inf, xend=min(df_list[[location]][[time1]][[time2]][R > 0.5]$Min_Line_Variance_Explained)*100, y = 0.5, yend = 0.5,linetype = "longdash", color = "black") +
		xlab("Minimum Percent Gene Expression\nExplained in Window") +
		ylab("Spearman rho")
	ggsave(plots_line[[location]][[time1]][[time2]], filename = gsub(" ", "_",paste0(outdir,location, "_", time1, "_",time2,"_window_spearman_line.png")), width = 4, height = 3)
	ggsave(plots_line[[location]][[time1]][[time2]], filename = gsub(" ", "_",paste0(outdir,location, "_", time1, "_",time2,"_window_spearman_line.pdf")), width = 4, height = 3)

	combinations$Percent_Pass_R0.5[row] <- min(df_list[[location]][[time1]][[time2]][R > 0.5]$Min_Line_Variance_Explained)*100
	combinations$Percent_sig[row] <- min(df_list[[location]][[time1]][[time2]][pvalue < 0.05/nrow(df_list[[location]][[time1]][[time2]])]$Min_Line_Variance_Explained)*100

}

saveRDS(df_list, paste0(outdir, "window_results_village_baseline.rds"))
fwrite(combinations, paste0(outdir, "window_results_village_baseline_summary.tsv"), sep = "\t")


df_sub_combined <- rbind(df_list[["Site 1"]][["Baseline"]][["Village"]][Window == 1], df_list[["Site 1"]][["Baseline"]][["Village"]][Window == 14340], df_list[["Site 1"]][["Baseline"]][["Village"]][Window == 14419])


combinations2 <- data.table(Location1 = rep(c(rep("Site 1",2), rep("Site 2",2), rep("Site 3",3), "Site 3 Cryopreserved"), 2),
							Location2 = rep(c("Site 2", "Site 3", "Site 1", "Site 3", "Site 1", "Site 2", "Site 3 Cryopreserved","Site 3"),2),
							Time = c(rep("Village",8),rep("Baseline", 8)))


### Make df ###
cor_list2 <- list()
df_list2 <- list()
plots_line2 <- list()


for (row in 1:nrow(combinations2)){
	location1 <- combinations2[row]$Location1
	location2 <- combinations2[row]$Location2
	time <- combinations2[row]$Time

	all_genes_line <- fits_df_line_list[[location1]][[time]]$Gene[fits_df_line_list[[location1]][[time]]$Gene %in% fits_df_line_list[[location2]][[time]]$Gene]
	df_list2[[location1]][[location2]][[time]] <- data.table(Window = seq(1:(length(all_genes_line) - 100)), R = 0, pvalue = NA, Min_Line_Variance_Explained = 0, Max_Line_Variance_Explained = 0)
	df_list2[[location1]][[location2]][[time]]$pvalue <- as.numeric(df_list2[[location1]][[location2]][[time]]$pvalue)

	for (window in df_list2[[location1]][[location2]][[time]]$Window){
		print(paste(row,window))
		genes_line <- all_genes_line[window:(window+100)]
		cor_list2[[location1]][[location2]][[time]][[window]] <- cor.test(fits_df_line_list[[location1]][[time]][Gene %in% genes_line]$Line, fits_df_line_list[[location2]][[time]][,Line[match(genes_line, Gene)]], method = "spearman")
		df_list2[[location1]][[location2]][[time]][Window == window, "R"] <- cor_list2[[location1]][[location2]][[time]][[window]]$estimate
		if (!is.null(cor_list2[[location1]][[location2]][[time]][[window]]$p.value)){
			df_list2[[location1]][[location2]][[time]][Window == window, "pvalue"] <- cor_list2[[location1]][[location2]][[time]][[window]]$p.value
		}
		df_list2[[location1]][[location2]][[time]][Window == window, "Min_Line_Variance_Explained"] <- min(fits_df_line_list[[location1]][[time]][Gene %in% genes_line]$Line)
		df_list2[[location1]][[location2]][[time]][Window == window, "Max_Line_Variance_Explained"] <- max(fits_df_line_list[[location1]][[time]][Gene %in% genes_line]$Line)
	}

	plots_line2[[location]][[time]][[time2]] <- ggplot(df_list2[[location1]][[location2]][[time]], aes(Min_Line_Variance_Explained*100, R, color = Max_Line_Variance_Explained*100)) +
		geom_point() +
		theme_classic() +
		scale_color_continuous_sequential(begin = 0, end = 1, palette = "SunsetDark", rev = FALSE, name = "Maximum\nPercent\nGene\nExpression\nExplained\nin\nWindow") +
		geom_segment(x=min(df_list2[[location1]][[location2]][[time]][R > 0.5]$Min_Line_Variance_Explained)*100, xend=min(df_list2[[location1]][[location2]][[time]][R > 0.5]$Min_Line_Variance_Explained)*100, y = -Inf, yend = 0.5,linetype = "longdash", color = "black")+
		geom_segment(x=-Inf, xend=min(df_list2[[location1]][[location2]][[time]][R > 0.5]$Min_Line_Variance_Explained)*100, y = 0.5, yend = 0.5,linetype = "longdash", color = "black") +
		xlab("Minimum Percent Gene Expression\nExplained in Window") +
		ylab("Spearman rho")
	ggsave(plots_line2[[location]][[time]][[time2]], filename = gsub(" ", "_",paste0(outdir,location1, "_", location2, "_",time,"_window_spearman_line.png")), width = 4, height = 3)
	ggsave(plots_line2[[location]][[time]][[time2]], filename = gsub(" ", "_",paste0(outdir,location1, "_", location2, "_",time,"_window_spearman_line.pdf")), width = 4, height = 3)

	
	combinations2$Percent_Pass_R0.5[row] <- min(df_list2[[location1]][[location2]][[time]][R > 0.5]$Min_Line_Variance_Explained)*100
	combinations2$Percent_sig[row] <- min(df_list2[[location1]][[location2]][[time]][pvalue < 0.05/nrow(df_list2[[location1]][[location2]][[time]])]$Min_Line_Variance_Explained)*100
}

saveRDS(df_list2, paste0(outdir, "window_results_between_site.rds"))
fwrite(combinations2, paste0(outdir, "window_results_between_site_summary.tsv"), sep = "\t")



combinations_merged <- rbind(data.table(combinations[c(Location %in% c("Site 1", "Site 2", "Site 3")),.(Location, Percent_Pass_R0.5)], Time = "Between Time"), data.table(combinations2[c(Location1 %in% c("Site 1", "Site 2", "Site 3") & c(Location2 %in% c("Site 1", "Site 2", "Site 3"))),.(Time, Percent_Pass_R0.5)], Location = NA))
combinations_merged$Time <- factor(gsub("Baseline", "Baseline\nBetween\nSite", combinations_merged$Time) %>% gsub("Village", "Village\nBetween\nSite", .) %>% gsub("Between Time", "Same Site\nBetween\nTime",.), levels = c("Baseline\nBetween\nSite", "Village\nBetween\nSite", "Same Site\nBetween\nTime"))

Rcomparison <- ggplot(combinations_merged, aes(Time, Percent_Pass_R0.5)) +
	geom_boxplot(aes = 0.6, outlier.shape = NA) +
	geom_jitter(width = 0.25) +
	theme_classic() +
	ylim(0,25) +
	theme(axis.title.x=element_blank()) +
	ylab("Expression Variance Explained for rho > 0.5")
ggsave(Rcomparison, filename = paste0(outdir, "minimum_variance_R0.05.png"), width = 2.5, height = 4.5)
ggsave(Rcomparison, filename = paste0(outdir, "minimum_variance_R0.05.pdf"), width = 3, height = 4)

pairwise.wilcox.test(combinations_merged$Percent_Pass_R0.5, combinations_merged$Time,p.adjust.method = "BH")
summary(aov(Percent_Pass_R0.5 ~ Time, data = combinations_merged))



combinations_merged_cryo <- rbind(data.table(combinations[c(Location %in% c("Site 3", "Site 3 Cryopreserved")),.(Location, Percent_Pass_R0.5)], Time = "Between Time"), data.table(combinations2[c(Location1 %in% c("Site 3", "Site 3 Cryopreserved") & c(Location2 %in% c("Site 3", "Site 3 Cryopreserved"))),.(Time, Percent_Pass_R0.5)], Location = NA))
combinations_merged_cryo$Time <- factor(gsub("Baseline", "Between\nSite", combinations_merged_cryo$Time) %>% gsub("Village", "Between\nSite", .) %>% gsub("Between Time", "Between\nTime",.), levels = c("Between\nSite", "Between\nTime"))



Rcomparison_cryo <- ggplot(combinations_merged_cryo, aes(Time, Percent_Pass_R0.5)) +
	geom_boxplot(aes = 0.6, outlier.shape = NA) +
	geom_jitter(width = 0.25) +
	theme_classic() +
	ylim(0,35) +
	theme(axis.title.x=element_blank()) +
	ylab("Expression Variance Explained for rho > 0.5")
ggsave(Rcomparison_cryo, filename = paste0(outdir, "minimum_variance_R0.05_cryo.png"), width = 2.5, height = 3.5)
ggsave(Rcomparison_cryo, filename = paste0(outdir, "minimum_variance_R0.05_cryo.pdf"), width = 2.5, height = 3.5)

pairwise.wilcox.test(combinations_merged_cryo$Percent_Pass_R0.5, combinations_merged_cryo$Time,p.adjust.method = "BH")
t.test(combinations_merged_cryo[Time == "Between\nTime"]$Percent_Pass_R0.5, combinations_merged_cryo[Time == "Between\nSite"]$Percent_Pass_R0.5)



### Order the genes by the village effect ###
village_df_list <- list()

for (row in 1:nrow(combinations)){
	print(row)
	location_x <- combinations[row]$Location.x
	location_y <- combinations[row]$Location.y

	all_genes <- fits_df_line_list[[location_x]]$Gene[fits_df_line_list[[location_x]]$Gene %in% fits_df_line_list[[location_y]]$Gene]
	village_df_list[[location_x]][[location_y]] <- data.table(Window = seq(1:(length(all_genes) - 1000)), R = 0, Min_Variance_Explained = 0)

	for (window in village_df_list[[location_x]][[location_y]]$Window){
		genes <- all_genes[window:(window+1000)]
		village_df_list[[location_x]][[location_y]][Window == window, "R"] <- cor(fits_df_line_list[[location_x]][Gene %in% genes]$Village, fits_df_line_list[[location_y]][,Village[match(genes, Gene)]], method = "spearman")
		village_df_list[[location_x]][[location_y]][Window == window, "Min_Variance_Explained"] <- min(fits_df_line_list[[location_x]][Gene %in% genes]$Village)
	}

}
