library(data.table)
library(Seurat)
library(tidyverse) 
library(colorspace)
library(powerEQTL)



dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
datadir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/nb_variance_partitioning_w_cryo/combined/"
faultdir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/nb_variance_partitioning_w_cryo/gene_separated/"
outdir <- paste0(dir,"output/Power_Calculations/")
dir.create(outdir, recursive = TRUE)


##### Set up colors #####
variable_colors <- c(Village = "#A2B0D0", Replicate = "#64A66B", Line = "#68319B") 
line_colors <- c(FSA0006 = "#F79E29", MBE1006 = "#9B2C99", TOB0421 = "#35369C")

site_updates <- c("Brisbane" = "Site 1", "Melbourne" = "Site 2" ,"Sydney" = "Site 3")

##### Test power for changes in std - test 3 stdev and plot #####
power_df_list <- list()

for (maf in seq(0.002, 0.5, 0.001)){
	for (n in seq(100,700,100)){
		for (sigma in c(0.13, 0.14,0.15)){
			power_df_list[[paste0(maf,n,sigma)]] <- data.table(MAF = maf, N = n, Sigma = sigma, Power = powerEQTL.SLR(
				MAF = maf,
				slope = 0.13,
				n = n,
				power = NULL,
				sigma.y = sigma,
				FWER = 0.05,
       			nTests = 2e+05,
			))
		}
	}
}

power_df <- do.call(rbind, power_df_list)

power_df_wide <- dcast(power_df, MAF + N ~ Sigma, value.var = "Power")
colnames(power_df_wide) <- c(colnames(power_df_wide)[1:2], paste0("sigma_", colnames(power_df_wide)[3:5]))

power_plot <- ggplot(power_df_wide, aes(color = factor(N, levels = seq(100,700,100)))) +
	# geom_point(aes(log10(MAF*100), sigma_0.13)) +
	# geom_point(aes(log10(MAF*100), sigma_0.16)) +
	# geom_point(aes(log10(MAF*100), sigma_0.19)) +
	geom_line(aes(log10(MAF*100), sigma_0.13), linetype = "solid") +
	geom_line(aes(log10(MAF*100), sigma_0.16), linetype = "dashed") +
	geom_line(aes(log10(MAF*100), sigma_0.19), linetype = "dotted") +
	ylab("Power") +
	geom_hline(yintercept = 0.8, linetype = "longdash") +
	theme_classic()
ggsave(power_plot, filename = paste0(outdir, "power_plot.png"), height = 3)


power_df$Sigma <- gsub("^", "Standard Deviation = ", power_df$Sigma)

power_plot_facet <- ggplot(power_df, aes(color = factor(N, levels = seq(100,700,10)))) +
	geom_line(aes(log10(MAF*100), Power), linetype = "solid") +
	facet_wrap(vars(Sigma), ncol = 1) +
	ylab("Power") +
	xlab("MAF (%)") +
	labs(color = "Sample\nSize (N)") +
	scale_color_discrete_sequential(palette = "ag_Sunset", ) +
	geom_hline(yintercept = 0.8, linetype = "longdash") +
	theme_classic() + 
	scale_y_continuous(breaks=seq(0, 1, 0.1)) +
	scale_x_continuous(breaks = c(log10(0.5),log10(1), log10(2), log10(5), log10(10), log10(20), log10(50)), labels = c(0.5, 1,2,5,10,20,50))
ggsave(power_plot_facet, filename = paste0(outdir, "power_plot_facet.png"), height = 7, width = 5)
ggsave(power_plot_facet, filename = paste0(outdir, "power_plot_facet.pdf"), height = 7, width = 5)




# ### Read in the seurat data objects
# seurat_file_list <- list.files("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/nb_partitioning_village_separate/data/", pattern = "_SCT_seurat_1pct.rds")
# seurat_list <- list()
# for (location in seurat_file_list){
# 	seurat_list[[location]] <- readRDS(paste0("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/nb_partitioning_village_separate/data/",location))
# }


# ### Separate by individual to look at sd within each separately ###
# seurat_list_subset <- lapply(seurat_list, function(x){
# 	tmp <- list()
# 	for (line in unique(x$Final_Assignment)){
# 		x$Time <- gsub("Village Day 4", "Village",x$Time)
# 		for (time in unique(x$Time)){
# 			tmp[[line]][[time]] <- subset(x, subset = (Final_Assignment == line & Time == time))
# 		}
# 	}
# 	return(tmp)
# })

# seurat_list_subset <- lapply(seurat_list_subset, function(x){
# 	lapply(x, function(y){
# 		lapply(y, function(z){
# 			subset(z, features = rownames(z)[which(rowSums(z[["SCT"]]@counts > 0)/ncol(z[["SCT"]]@counts) >= 0.01)])
# 		})
# 	})
# })

# seurat_list_subset <- lapply(seurat_list_subset, function(x){
# 	lapply(x, function(y){
# 		lapply(y, function(z){
# 			SCTransform(z, verbose = TRUE, vars.to.regress = c("scores.G1", "scores.S", "scores.G2M", "percent.mt", "percent.rb"), return.only.var.genes = FALSE)
# 		})
# 	})
# })

# ### Need to have the same number in each group for sd calculations
# ## Identify the minimum number of genes for each baseline-village pair
# number_list <- lapply(seurat_list_subset, function(x){
# 	lapply(x, function(y){
# 		min(ncol(y[[1]]), ncol(y[[2]]))
# 	})
# })


# ### Calculate SD for each gene ###
# sd_list <- lapply(names(seurat_list_subset), function(x){
# 	tmp <- lapply(names(seurat_list_subset[[x]]), function(y){
# 		lapply(seurat_list_subset[[x]][[y]], function(z){
# 			print(z)
# 			if (ncol(z[["SCT"]]@counts) > number_list[[x]][[y]]){ ### subsample for the longer dataframe to make comparison fair
# 				apply(z[["SCT"]]@counts[,sample(1:ncol(z[["SCT"]]@counts), number_list[[x]][[y]])],1, sd, na.rm = TRUE)
# 			} else {
# 				apply(z[["SCT"]]@counts,1, sd, na.rm = TRUE)
# 			}
# 		})
# 	})
# 	names(tmp) <- names(seurat_list_subset[[x]])
# 	return(tmp)
# })
# names(sd_list) <- names(seurat_list_subset)

# mean_list <- lapply(names(seurat_list_subset), function(x){
# 	tmp <- lapply(names(seurat_list_subset[[x]]), function(y){
# 		lapply(seurat_list_subset[[x]][[y]], function(z){
# 			print(z)
# 			if (ncol(z[["SCT"]]@counts) > number_list[[x]][[y]]){ ### subsample for the longer dataframe to make comparison fair
# 				apply(z[["SCT"]]@counts[,sample(1:ncol(z[["SCT"]]@counts), number_list[[x]][[y]])],1, mean, na.rm = TRUE)
# 			} else {
# 				apply(z[["SCT"]]@counts,1, sd, na.rm = TRUE)
# 			}
# 		})
# 	})
# 	names(tmp) <- names(seurat_list_subset[[x]])
# 	return(tmp)
# })
# names(mean_list) <- names(seurat_list_subset)

# lapply(sd_list, function(x) lapply(x, function(y) lapply(y, min)))
# lapply(sd_list, function(x) lapply(x, function(y) lapply(y, max)))
# lapply(sd_list, function(x) lapply(x, function(y) lapply(y, mean)))

# ### for variable features
# lapply(names(sd_list), function(x) lapply(names(sd_list[[x]]), function(y) lapply(names(sd_list[[x]][[y]]), function(z) min(sd_list[[x]][[y]][[z]][seurat_list_subset[[x]][[y]][[z]]$SCT@var.features]))))
# lapply(names(sd_list), function(x) lapply(names(sd_list[[x]]), function(y) lapply(names(sd_list[[x]][[y]]), function(z) max(sd_list[[x]][[y]][[z]][seurat_list_subset[[x]][[y]][[z]]$SCT@var.features]))))
# lapply(names(sd_list), function(x) lapply(names(sd_list[[x]]), function(y) lapply(names(sd_list[[x]][[y]]), function(z) mean(sd_list[[x]][[y]][[z]][seurat_list_subset[[x]][[y]][[z]]$SCT@var.features]))))


# ### Combine the data into a data.table to plot and test ###
# sd_df_list <- lapply(names(sd_list), function(x){
# 	tmp2 <- lapply(names(sd_list[[x]]), function(y){
# 		tmp <- lapply(names(sd_list[[x]][[y]]), function(z){
# 			data.table(Gene = names(sd_list[[x]][[y]][[z]]), SD = sd_list[[x]][[y]][[z]], Site = gsub("_SCT_seurat_1pct.rds", "", x), Cell_Line = y, Time = z)
# 		})
# 		do.call(rbind, tmp)
# 	})
# 	do.call(rbind,tmp2)
# })

# sd_df <- do.call(rbind, sd_df_list)

# mean_df_list <- lapply(names(mean_list), function(x){
# 	tmp2 <- lapply(names(mean_list[[x]]), function(y){
# 		tmp <- lapply(names(mean_list[[x]][[y]]), function(z){
# 			data.table(Gene = names(mean_list[[x]][[y]][[z]]), mean = mean_list[[x]][[y]][[z]], Site = gsub("_SCT_seurat_1pct.rds", "", x), Cell_Line = y, Time = z)
# 		})
# 		do.call(rbind, tmp)
# 	})
# 	do.call(rbind,tmp2)
# })

# mean_df <- do.call(rbind, mean_df_list)


# sd_df_wide <- dcast(sd_df, Gene + Site + Cell_Line ~ Time, value.var = "SD")
# na.omit(sd_df_wide)

# mean_df_wide <- dcast(mean_df, Gene + Site + Cell_Line ~ Time, value.var = "mean")
# na.omit(mean_df_wide)

# p_boxplot_facet <- ggplot(sd_df, aes(Time, SD)) +
# 						geom_boxplot() +
# 						theme_classic() +
# 						facet_grid(Cell_Line ~ Site)
# ggsave(p_boxplot_facet, filename = paste0(outdir,"sd_boxplot_facet.png"))

# p_density_facet <- ggplot(sd_df, aes(log(SD), color = Time)) +
# 						geom_density() +
# 						theme_classic() +
# 						facet_grid(Cell_Line ~ Site)
# ggsave(p_density_facet, filename = paste0(outdir,"sd_density_facet.png"))

# p_mean_density_facet <- ggplot(mean_df, aes(log(mean), color = Time)) +
# 						geom_density() +
# 						theme_classic() +
# 						facet_grid(Cell_Line ~ Site)
# ggsave(p_mean_density_facet, filename = paste0(outdir,"mean_density_facet.png"))

# ### Subset for genes where < 10% of variance is described by cell line ###
# fits_df <- read_delim("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/nb_variance_partitioning_w_cryo/combined/nb_variance_partitionin_df.tsv", delim = "\t")


# ##### Remove faulty model genes and locations #####
# ### Read in faulty log files ###
# faulty <- fread("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/nb_variance_partitioning_w_cryo/gene_separated/faulty_runs.txt", header = FALSE)


# ### Make dataframe of location x number for faulty models ###
# faulty_df <- data.table(Location = gsub("_nb\\.o.+", "", faulty$V1), Number = gsub(".+nb\\.o\\d+\\.", "", faulty$V1))

# ## Add gene id to the faulty dataframe
# faulty_df$ENSG <- NA
# for (row in 1:nrow(faulty_df)){
# 	print(row)
# 	faulty_df$ENSG[row] <- rownames(seurat_list[[faulty_df$Location[row]]])[as.numeric(faulty_df$Number[row])] 
# }


# ### Remove the faulty genes ###
# fits_df <- setDT(fits_df)[!faulty_df, on = c("Location", "Gene" = "ENSG")]

# ### Remove mt genes
# MtGeneList <- read.delim(file = "/directflow/SCCGGroupShare/projects/DrewNeavin/References/MtGeneList_GeneID_ENSG.txt")
# fits_df_noMT <- fits_df[-c(which(fits_df$Gene %in% MtGeneList$ENSG)),]

# sd_df_wide <- fits_df_noMT[sd_df_wide, on = c("Gene", "Location" = "Site")]
# sd_df_wide <- na.omit(sd_df_wide)

# mean_df_wide <- fits_df_noMT[mean_df_wide, on = c("Gene", "Location" = "Site")]
# mean_df_wide <- na.omit(mean_df_wide)


# # sd_df_low_list <- list()

# # for (perc in seq(0.01, 0.1, 0.01)){
# # 	sd_df_low_list <- 
# # }
# sd_df_5_10 <- sd_df_wide[0.05 < Line & Line < 0.1]
# sd_df_5_10$Difference <- sd_df_5_10$Baseline - sd_df_5_10$i.Village
# sd_df_5_10 <- sd_df_5_10[order(sd_df_5_10$Difference),]

# mean_df_5_10 <- mean_df_wide[0.05 < Line & Line < 0.1]
# mean_df_5_10$Difference <- mean_df_5_10$Baseline - sd_df_5_10$i.Village
# mean_df_5_10 <- mean_df_5_10[order(mean_df_5_10$Difference),]

# sd_df_5_10_long <- melt(sd_df_5_10, id.vars = colnames(sd_df_5_10)[1:6], measure.vars = c("Baseline", "i.Village"), variable.name = "Time", value.name = "SD")
# mean_df_5_10_long <- melt(mean_df_5_10, id.vars = colnames(mean_df_5_10)[1:6], measure.vars = c("Baseline", "i.Village"), variable.name = "Time", value.name = "SD")


# p_density_facet_5_10 <- ggplot(sd_df_5_10_long, aes(log(SD), color = Time)) +
# 						geom_density() +
# 						theme_classic() +
# 						facet_grid(Cell_Line ~ Location)
# ggsave(p_density_facet_5_10, filename = paste0(outdir,"sd_boxplot_facet_5_10.png"))


# p_density_facet_5_10_dif <- ggplot(sd_df_5_10, aes(Baseline-i.Village)) +
# 						geom_density() +
# 						theme_classic() +
# 						facet_grid(Cell_Line ~ Location)
# ggsave(p_density_facet_5_10_dif, filename = paste0(outdir,"sd_density_facet_5_10_diff.png"))


# wilcox_result_5_10 <- list()
# spearman_result_5_10 <- list()


# for (site in unique(sd_df_5_10$Location)){
# 	for (line in unique(sd_df_5_10$Cell_Line)){
# 		wilcox_result_5_10[[site]][[line]] <- wilcox.test(sd_df_5_10[Location == site & Cell_Line == line]$Baseline, sd_df_5_10[Location == site & Cell_Line == line]$i.Village, paired = TRUE, alternative = "less")
# 		spearman_result_5_10[[site]][[line]] <- cor.test(sd_df_5_10[Location == site & Cell_Line == line]$Baseline,sd_df_5_10[Location == site & Cell_Line == line]$i.Village, method = "spearman")
# 	}
# }

# wilcox_result <- list()

# for (site in unique(sd_df_wide$Location)){
# 	for (line in unique(sd_df_wide$Cell_Line)){
# 		wilcox_result[[site]][[line]] <- wilcox.test(sd_df_wide[Location == site & Cell_Line == line]$Baseline, sd_df_wide[Location == site & Cell_Line == line]$i.Village, paired = TRUE, alternative = "less")
# 	}
# }

# wilcox_mean_result_5_10 <- list()
# spearman_mean_result_5_10 <- list()

# for (site in unique(mean_df_5_10$Location)){
# 	for (line in unique(mean_df_5_10$Cell_Line)){
# 		wilcox_mean_result_5_10[[site]][[line]] <- wilcox.test(mean_df_5_10[Location == site & Cell_Line == line]$Baseline, mean_df_5_10[Location == site & Cell_Line == line]$i.Village, paired = TRUE, alternative = "two.sided")
# 		`[[site]][[line]] <- cor(mean_df_5_10[Location == site & Cell_Line == line]$Baseline,mean_df_5_10[Location == site & Cell_Line == line]$i.Village, method = "spearman")
# 	}
# }


# wilcox_mean_result <- list()

# for (site in unique(sd_df_wide$Location)){
# 	for (line in unique(sd_df_wide$Cell_Line)){
# 		wilcox_mean_result[[site]][[line]] <- wilcox.test(sd_df_wide[Location == site & Cell_Line == line]$Baseline, sd_df_wide[Location == site & Cell_Line == line]$i.Village, paired = TRUE, alternative = "less")
# 	}
# }


