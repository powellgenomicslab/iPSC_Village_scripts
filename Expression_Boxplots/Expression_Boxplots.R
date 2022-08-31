library(tidyverse)
library(Seurat)
library(ggplot2)
library(ggforce)
library(ggsignif)
library(data.table)
library(glmmTMB)
library(ggnewscale)
library(MASS)
library(ggpubr)



dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
datadir <- paste0(dir,"output/nb_partitioning_village_separate/data/")
outdir <- paste0(dir,"output/Expression_Boxplots/")

dir.create(outdir)


##### Make funciton to make dataframees for each gene
expression_df <- function(seurat_list, gene, meta_data){
	Expression_list <- lapply(seurat_list, function(x){
		data.frame(Counts = x[["SCT"]]@counts[gene,], Normalized = x[["SCT"]]@scale.data[gene,], x@meta.data[,c(meta_data)])
	})
	do.call(rbind, Expression_list)
}



save_figs <- function(plot, basename, width = 17, height = 17, units = "cm"){
    ggsave(plot, filename = paste0(basename,".png"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".pdf"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".eps"), height = height, width = width, units = units)
}




##### Set up colors #####
variable_colors <- c(Village = "#A2B0D0", Replicate = "#64A66B", Line = "#68319B") 
line_colors <- c(FSA0006 = "#F79E29", MBE1006 = "#9B2C99", TOB0421 = "#35369C")
village_colors <- c("Uni-Culture" = "#613246", "Village" = "#A286AA")

site_updates <- c("Brisbane" = "Site 1" ,"Melbourne" = "Site 2", "Sydney" = "Site 3")



##### Read in Data #####
files <- list.files(paste0(datadir))

seurat_list <- lapply(files, function(x){
	readRDS(paste0(datadir,x))
})
names(seurat_list) <- gsub("_seurat.rds", "", files)


seurat_list <- lapply(seurat_list, function(x){
	x$Location <- ifelse(x$Time %in% c("Thawed Village Day 0", "Thawed Village Day 7"), "Sydney_Cryopreserved", x$Location)
	x$Time <- ifelse(x$Time %in% c("Village Day 4", "Thawed Village Day 7"), "Village", "Uni-Culture")
	for (location in names(site_updates)){
		x$Location <- gsub(location, site_updates[location], x$Location)
	}
	return(x)
})





##### CHCHD2 #####
CHCHD2 <- expression_df(seurat_list, "ENSG00000106153", c("Location", "Time", "Final_Assignment", "Site_rep"))
CHCHD2$Site_rep <- gsub("[a-z]+", "", CHCHD2$Site_rep) %>% gsub("[A-Z]", "", .)
dir.create(paste0(outdir,"CHCHD2"))

pCHCHD2_counts <- ggplot(CHCHD2[which(CHCHD2$Location != "Site 3_Cryopreserved"),], aes(x = Time, y = Counts, color = Final_Assignment)) +
					geom_boxplot(outlier.size = 0.5) +
					theme_classic() +
					facet_wrap(vars(Location), scales = "free_y", ncol = 1, strip.position = "right") +
					scale_color_manual(values = line_colors) +
					ylab("CHCHD2 Counts")
save_figs(pCHCHD2_counts, paste0(outdir,"CHCHD2/counts"), width = 9, height = 6)

pCHCHD2_scaled <- ggplot(CHCHD2[which(CHCHD2$Location != "Site 3_Cryopreserved"),], aes(x = Time, y = Normalized, color = Final_Assignment)) +
					geom_boxplot(outlier.size = 0.5) +
					theme_classic() +
					facet_wrap(vars(Location), scales = "free_y", ncol = 1, strip.position = "right") +
					scale_color_manual(values = line_colors)
save_figs(pCHCHD2_scaled, paste0(outdir,"CHCHD2/normalized"), width = 9, height = 6)


CHCHD2_cryo <- CHCHD2[grepl("Site 3", CHCHD2$Location),]
CHCHD2_cryo$Location <- gsub("Site 3_Cryopreserved", "Cryopreserved", CHCHD2_cryo$Location) %>% gsub("Site 3", "Fresh", .)


pCHCHD2_counts_cryo <- ggplot(CHCHD2_cryo, aes(x = Time, y = Counts, color = Final_Assignment)) +
					geom_boxplot(outlier.size = 0.5) +
					theme_classic() +
					facet_wrap(vars(Location), scales = "free_y", ncol = 1, strip.position = "right") +
					scale_color_manual(values = line_colors) +
					ylab("CHCHD2 Counts")
save_figs(pCHCHD2_counts_cryo, paste0(outdir,"CHCHD2/counts_cryo"), width = 9, height = 6)

pCHCHD2_scaled_cryo <- ggplot(CHCHD2_cryo, aes(x = Time, y = Normalized, color = Final_Assignment)) +
					geom_boxplot(outlier.size = 0.5) +
					theme_classic() +
					facet_wrap(vars(Location), scales = "free_y", ncol = 1, strip.position = "right") +
					scale_color_manual(values = line_colors)
save_figs(pCHCHD2_scaled_cryo, paste0(outdir,"CHCHD2/normalized_cryo"), width = 9, height = 6)



### Line Effect ###
pCHCHD2_counts_line <- ggplot(CHCHD2[which(CHCHD2$Location == "Site 1"),], aes(x = Final_Assignment, y = Counts, color = Final_Assignment)) +
					geom_boxplot(outlier.size = 0.5) +
					theme_classic() +
					scale_color_manual(values = line_colors) +
					ylab("Gene Counts") +
					xlab("hiPSC Line") +
					theme(legend.position = "none",
						axis.text.x=element_blank(),
						axis.ticks.x=element_blank()) +
					# geom_smooth(method = "lm", inherit.aes = FALSE, aes(x = Final_Assignment, y = Counts), color = "black", size = 0.5) +
					# geom_text(aes(label="Variance Explained = 95.5%"), x = -Inf, y = Inf,  hjust = -0.1, vjust = 1.4,inherit.aes = FALSE, size = 2)
save_figs(pCHCHD2_counts_line, paste0(outdir,"CHCHD2/counts_line"), width = 4.5, height = 6)

### Village Effect ###
pCHCHD2_counts_village <- ggplot(CHCHD2[which(CHCHD2$Location != "Site 1"),], aes(x = Time, y = Counts)) +
					geom_boxplot(outlier.size = 0.5) +
					theme_classic() +
					ylab("CHCHD2 Counts") +
					theme(legend.position = "none",
						axis.title.x = element_blank()) +
					# geom_text(aes(label = "Variance Explained = 2.3%"), x = -Inf, y = Inf,  hjust = -0.1, vjust = 1.4,inherit.aes = FALSE, size = 2)
save_figs(pCHCHD2_counts_village, paste0(outdir,"CHCHD2/counts_time"), width = 4.5, height = 6)

### Replicate Effect ###
pCHCHD2_counts_village <- ggplot(CHCHD2[which(CHCHD2$Location != "Site 1"),], aes(x = Site_rep, y = Counts)) +
					geom_boxplot(outlier.size = 0.5) +
					theme_classic() +
					ylab("CHCHD2 Counts") +
					xlab("Replicate") +
					theme(legend.position = "none",
						axis.text.x=element_blank(),
						axis.ticks.x=element_blank()) +
					# geom_text(aes(label = "Variance Explained = 2.1%"), x = -Inf, y = Inf,  hjust = -0.1, vjust = 1.4,inherit.aes = FALSE, size = 2)
save_figs(pCHCHD2_counts_village, paste0(outdir,"CHCHD2/counts_effect"), width = 4.5, height = 6)


##### CHCHD2 eQTL Fig #####
CHCHD2 <- expression_df(seurat_list, "ENSG00000106153", c("Location", "Time", "Final_Assignment", "Site_rep"))
### 7:56049019 SNP genotype (genotyped SNP in LD block)
CHCHD2_geno <- data.frame("Final_Assignment" = c("FSA0006", "MBE1006", "TOB0421"), "Genotype" = c("A/A", "A/A", "A/G"))
CHCHD2 <- left_join(CHCHD2, CHCHD2_geno)
CHCHD2 <- data.table(CHCHD2)

CHCHD2_means <- CHCHD2[,.(Mean=mean(Counts)),.(Final_Assignment, Genotype, Time, Site_rep, Location)]
CHCHD2_means$Genotype_Numeric <- ifelse(CHCHD2_means$Genotype == "A/A", 0, 1)

CHCHD2_means_norm <- CHCHD2[,.(Mean=mean(Normalized)),.(Final_Assignment, Genotype, Time, Site_rep, Location)]
CHCHD2_means_norm$Genotype_Numeric <- ifelse(CHCHD2_means$Genotype == "A/A", 0, 1)

dir.create(paste0(outdir,"CHCHD2"))

### Calculate t test significance ###
t_test_df <- unique(data.table(Location = CHCHD2$Location, Time = CHCHD2$Time))
t_test_df$Beta_nb <- 0
t_test_df$Beta_lm <- 0

t_test_norm_df <- unique(data.table(Location = CHCHD2$Location, Time = CHCHD2$Time))
t_test_norm_df$Beta_nb <- 0
t_test_norm_df$Beta_lm <- 0

model_nb <- list()
model_lm <- list()
model_nb_norm <- list()
model_lm_norm <- list()

for (row in 1:nrow(t_test_df)){
	model_nb[[paste0(t_test_df$Location[row], "_", t_test_df$Time[row])]] <- summary(glmmTMB(Mean ~ Genotype_Numeric + (1|Site_rep), data = CHCHD2_means[Location == t_test_df$Location[row] & Time == t_test_df$Time[row]], family = nbinom2))
	model_lm[[paste0(t_test_df$Location[row], "_", t_test_df$Time[row])]] <- summary(lm(Mean ~ Genotype_Numeric, data = CHCHD2_means[Location == t_test_df$Location[row] & Time == t_test_df$Time[row]]))
	model_lm_norm[[paste0(t_test_df$Location[row], "_", t_test_df$Time[row])]] <- summary(lm(Mean ~ Genotype_Numeric, data = CHCHD2_means_norm[Location == t_test_df$Location[row] & Time == t_test_df$Time[row]]))
	t_test_df[row]$Beta_nb <- round(model_nb[[paste0(t_test_df$Location[row], "_", t_test_df$Time[row])]]$coefficients$cond[2,"Estimate"], 1)
	t_test_df[row]$Beta_lm <- round(model_lm[[paste0(t_test_df$Location[row], "_", t_test_df$Time[row])]]$coefficients[2,"Estimate"], 1)
	t_test_norm_df[row]$Beta_lm <- round(model_lm_norm[[paste0(t_test_norm_df$Location[row], "_", t_test_norm_df$Time[row])]]$coefficients[2,"Estimate"], 1)
}



pCHCHD2_counts_eQTL <- ggplot(CHCHD2[which(CHCHD2$Location != "Site 3_Cryopreserved"),], aes(x = Genotype, y = Counts)) +
					geom_boxplot(outlier.size = 0.3) +
					theme_classic() +
					facet_grid(Location ~ Time) +
					ylab("CHCHD2 Counts")+
					xlab("rs2304376 Genotype") +
					geom_signif(comparisons = list(c("A/A", "A/G")), 
								map_signif_level=TRUE, y = max(CHCHD2[which(CHCHD2$Location != "Site 3_Cryopreserved"),]$Counts) + 5,
								test = "t.test")+
					ylim(0,max(CHCHD2[which(CHCHD2$Location != "Site 3_Cryopreserved"),]$Counts) + 25)
save_figs(pCHCHD2_counts_eQTL, paste0(outdir,"CHCHD2/counts_eQTL"), width = 5, height = 6)



pCHCHD2_counts_eQTL_mean <- ggplot(CHCHD2_means[which(CHCHD2_means$Location != "Site 3_Cryopreserved"),], aes(x = Genotype_Numeric, y = Mean, color = Final_Assignment)) +
					geom_point(size = 0.2, position = position_jitter(width = 0.05)) +
					theme_classic() +
					facet_grid(Location ~ Time) +
					ylab("CHCHD2 Counts")+
					xlab("rs2304376 Genotype") +
					scale_color_manual(values = line_colors) +
					new_scale("color") +
					geom_smooth(method = "lm", inherit.aes = FALSE, aes(x = Genotype_Numeric, y = Mean), color = "black", size = 0.5) +
					scale_x_continuous(breaks  = c(0,0.5), labels=c("A/A", "A/G"))  +
					geom_text(data = t_test_df[which(t_test_df$Location != "Site 3_Cryopreserved"),], aes(label=paste0("Beta = ",Beta_lm)), x = -Inf, y = Inf,  hjust = -0.1, vjust = 1.4,
            inherit.aes = FALSE, size = 2)
save_figs(pCHCHD2_counts_eQTL_mean, paste0(outdir,"CHCHD2/counts_eQTL_mean"), width = 9, height = 6)


pCHCHD2_counts_eQTL_mean_norm <- ggplot(CHCHD2_means_norm[which(CHCHD2_means_norm$Location != "Site 3_Cryopreserved"),], aes(x = Genotype_Numeric, y = Mean, color = Final_Assignment)) +
					geom_point(size = 1, position = position_jitter(width = 0.05)) +
					theme_classic() +
					facet_grid(Location ~ Time, scales = "free_y") +
					ylab("CHCHD2 Normalized Counts")+
					xlab("rs2304376 Genotype") +
					scale_color_manual(values = line_colors) +
					new_scale("color") +
					geom_smooth(method = "lm", inherit.aes = FALSE, aes(x = Genotype_Numeric, y = Mean), color = "black", size = 0.5) +
					scale_x_continuous(breaks  = c(0,0.5), labels=c("A/A", "A/G"))  +
					geom_text(data = t_test_norm_df[which(t_test_norm_df$Location != "Site 3_Cryopreserved"),], aes(label=paste0("beta = ", Beta_lm)), x = -Inf, y = Inf,  hjust = -0.1, vjust = 1.4,inherit.aes = FALSE, size = 2.5)
save_figs(pCHCHD2_counts_eQTL_mean_norm, paste0(outdir,"CHCHD2/counts_eQTL_mean_norm"), width = 9, height = 6)



CHCHD2_cryo <- CHCHD2[grepl("Site 3", CHCHD2$Location),]
CHCHD2_cryo$Location <- gsub("Site 3_Cryopreserved", "Cryopreserved", CHCHD2_cryo$Location) %>% gsub("Site 3", "Fresh", .)


pCHCHD2_counts_cryo_eQTL <- ggplot(CHCHD2_cryo, aes(x = Genotype, y = Counts)) +
					geom_boxplot(outlier.size = 0.3) +
					theme_classic() +
					facet_grid(Location ~ Time) +
					ylab("CHCHD2 Counts") +
					xlab("rs2304376 Genotype") +
					geom_signif(comparisons = list(c("A/A", "A/G")), 
								map_signif_level=TRUE, y = max(CHCHD2[which(CHCHD2$Location != "Site 3_Cryopreserved"),]$Counts) + 5,
								test = "t.test")+
					ylim(0,max(CHCHD2[which(CHCHD2$Location != "Site 3_Cryopreserved"),]$Counts) + 25)
save_figs(pCHCHD2_counts_cryo_eQTL, paste0(outdir,"CHCHD2/counts_cryo_eQTL"), width = 5, height = 7)


CHCHD2_means_norm_cryo <- CHCHD2_means_norm[grepl("Site 3", CHCHD2_means_norm$Location),]
CHCHD2_means_norm_cryo$Location <- gsub("Site 3_Cryopreserved", "Cryopreserved", CHCHD2_means_norm_cryo$Location) %>% gsub("Site 3", "Fresh", .)
CHCHD2_means_norm_cryo$Location <- factor(CHCHD2_means_norm_cryo$Location, levels = c("Fresh", "Cryopreserved"))

t_test_norm_df_cryo <- t_test_norm_df[grepl("Site 3", t_test_norm_df$Location),]
t_test_norm_df_cryo$Location <- gsub("Site 3_Cryopreserved", "Cryopreserved", t_test_norm_df_cryo$Location) %>% gsub("Site 3", "Fresh", .)



pCHCHD2_counts_eQTL_mean_norm_cryo <- ggplot(CHCHD2_means_norm_cryo, aes(x = Genotype_Numeric, y = Mean, color = Final_Assignment)) +
					geom_point(size = 1, position = position_jitter(width = 0.05)) +
					theme_classic() +
					facet_grid(factor(Location, levels = c("Fresh", "Cryopreserved")) ~ Time) +
					ylab("CHCHD2 Normalized Counts")+
					xlab("rs2304376 Genotype") +
					scale_color_manual(values = line_colors) +
					new_scale("color") +
					geom_smooth(method = "lm", inherit.aes = FALSE, aes(x = Genotype_Numeric, y = Mean), color = "black", size = 0.5) +
					scale_x_continuous(breaks  = c(0,0.5), labels=c("A/A", "A/G"))  +
					geom_text(data = t_test_norm_df_cryo, aes(label=paste0("beta = ", Beta_lm)), x = -Inf, y = Inf,  hjust = -0.1, vjust = 1.4,inherit.aes = FALSE, size = 2.5)
save_figs(pCHCHD2_counts_eQTL_mean_norm_cryo, paste0(outdir,"CHCHD2/counts_eQTL_mean_norm_cryo"), width = 10, height = 7)



##### Quintile Data #####
datadir_quint <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/RNAvelocity/nb_partitioning_rep_separate/data/"

quint_files <- list.files(datadir_quint)


seurat_quint_list <- lapply(quint_files, function(x){
	readRDS(paste0(datadir_quint,x))
})
names(seurat_quint_list) <- gsub("_seurat_1pct.rds", "", quint_files)


seurat_quint_list <- lapply(seurat_quint_list, function(x){
	x$Location <- ifelse(x$Time %in% c("Thawed Village Day 0", "Thawed Village Day 7"), "Sydney_Cryopreserved", x$Location)
	x$Time <- ifelse(x$Time %in% c("Village Day 4", "Thawed Village Day 7", "Village"), "Village", "Baseline")
	x$Time
	for (location in names(site_updates)){
		x$Location <- gsub(location, site_updates[location], x$Location)
	}
	return(x)
})




##### RARRES2 eQTL #####
RARRES2 <- expression_df(seurat_quint_list, "ENSG00000106538", c("Location", "Time", "Final_Assignment", "Site_rep", "Quintile"))
RARRES2$Location <-factor(gsub("cryo", " Cryo", RARRES2$Location), levels = c("Site 1", "Site 2", "Site 3", "Site 3 Cryopreserved"))

### 7:56049019 SNP genotype (genotyped SNP in LD block)
SNP_vcf <- fread("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/rs2108851_RARRES2.vcf") ### ref = T, alt = C
colnames(SNP_vcf) <- gsub("22_FSA", "FSA0006", colnames(SNP_vcf)) %>% gsub("29_MBE", "MBE1006", .) %>% gsub("36_TOB00421_i_E8", "TOB0421", .)

SNP_vcf$FSA0006 <- gsub(":\\d\\.\\d\\d\\d,\\d\\.\\d\\d\\d,\\d", "", SNP_vcf$FSA0006) %>% gsub(".+:", "", .)
SNP_vcf$MBE1006 <- gsub(":\\d\\,\\d\\,\\d", "", SNP_vcf$MBE1006) %>% gsub(".+:", "", .)
SNP_vcf$TOB0421 <- gsub(":\\d\\,\\d\\,\\d", "", SNP_vcf$TOB0421) %>% gsub(".+:", "", .)

SNP_vcf_long = melt(SNP_vcf, id.vars = c("#CHROM", "POS", "ID", "REF", "ALT"), measure.vars = c("FSA0006", "MBE1006", "TOB0421"),  variable.name = "Final_Assignment")

SNP_vcf_long$Genotype <- ifelse(SNP_vcf_long$value < 0.5, paste0(SNP_vcf_long$REF, "/", SNP_vcf_long$REF), 
							ifelse(SNP_vcf_long$value > 1.5, paste0(SNP_vcf_long$ALT, "/", SNP_vcf_long$ALT), paste0(SNP_vcf_long$REF, "/", SNP_vcf_long$ALT)))

SNP_vcf_long$Genotype_Numeric <- ifelse(SNP_vcf_long$value < 0.5, 0, 
									ifelse(SNP_vcf_long$value > 1.5, 2, 1))


RARRES2 <- data.table(RARRES2)
RARRES2 <- RARRES2[SNP_vcf_long, on = "Final_Assignment"]

RARRES2_means <- RARRES2[,.(Mean=mean(Counts)),.(Final_Assignment, Genotype, Genotype_Numeric, Time, Site_rep, Location, Quintile)]

RARRES2_means_norm <- RARRES2[,.(Mean=mean(Normalized)),.(Final_Assignment, Genotype, Genotype_Numeric, Time, Site_rep, Location, Quintile)]


dir.create(paste0(outdir,"RARRES2"))


### Calculate t test significance ###
RARRES2_t_test_df <- unique(data.table(Location = RARRES2$Location, Time = RARRES2$Time, Quintile = RARRES2$Quintile))
RARRES2_t_test_df$Beta_nb <- 0
RARRES2_t_test_df$Beta_lm <- 0

RARRES2_t_test_norm_df <- unique(data.table(Location = RARRES2$Location, Time = RARRES2$Time, Quintile = RARRES2$Quintile))
RARRES2_t_test_norm_df$Beta_lm <- 0

model_nb <- list()
model_lm <- list()
model_nb_norm <- list()
model_lm_norm <- list()

for (row in 1:nrow(RARRES2_t_test_df)){
	model_nb[[paste0(RARRES2_t_test_df$Location[row], "_", RARRES2_t_test_df$Time[row])]] <- summary(glmmTMB(Mean ~ Genotype_Numeric + (1|Site_rep), data = RARRES2_means[Location == RARRES2_t_test_df$Location[row] & Time == RARRES2_t_test_df$Time[row] & Quintile == RARRES2_t_test_df$Quintile[row]], family = nbinom2))
	model_lm[[paste0(RARRES2_t_test_df$Location[row], "_", RARRES2_t_test_df$Time[row])]] <- summary(lm(Mean ~ Genotype_Numeric, data = RARRES2_means[Location == RARRES2_t_test_df$Location[row] & Time == RARRES2_t_test_df$Time[row] & Quintile == RARRES2_t_test_df$Quintile[row]]))
	model_lm_norm[[paste0(RARRES2_t_test_df$Location[row], "_", RARRES2_t_test_df$Time[row])]] <- summary(lm(Mean ~ Genotype_Numeric, data = RARRES2_means_norm[Location == RARRES2_t_test_df$Location[row] & Time == RARRES2_t_test_df$Time[row] & Quintile == RARRES2_t_test_df$Quintile[row]]))
	RARRES2_t_test_df[row]$Beta_nb <- round(model_nb[[paste0(RARRES2_t_test_df$Location[row], "_", RARRES2_t_test_df$Time[row])]]$coefficients$cond[2,"Estimate"], 1)
	RARRES2_t_test_df[row]$Beta_lm <- round(model_lm[[paste0(RARRES2_t_test_df$Location[row], "_", RARRES2_t_test_df$Time[row])]]$coefficients[2,"Estimate"], 1)
	RARRES2_t_test_norm_df[row]$Beta_lm <- round(model_lm_norm[[paste0(RARRES2_t_test_norm_df$Location[row], "_", RARRES2_t_test_norm_df$Time[row])]]$coefficients[2,"Estimate"], 1)
}

RARRES2_t_test_df[RARRES2_t_test_df == 0] <- NA
RARRES2_t_test_norm_df[RARRES2_t_test_norm_df == 0] <- NA

pRARRES2_counts_eQTL_mean <- list()

for (location in levels(RARRES2_means$Location)){
	pRARRES2_counts_eQTL_mean[[location]] <- ggplot(RARRES2_means[which(RARRES2_means$Location == location),], aes(x = Genotype_Numeric, y = Mean, color = Final_Assignment)) +
						geom_point(size = 0.2, position = position_jitter(width = 0.05)) +
						theme_classic() +
						facet_grid(Quintile ~ Time) +
						ylab("CHCHD2 Counts")+
						xlab("rs2304376 Genotype") +
						scale_color_manual(values = line_colors) +
						new_scale("color") +
						geom_smooth(method = "lm", inherit.aes = FALSE, aes(x = Genotype_Numeric, y = Mean), color = "black", size = 0.5) +
						scale_x_continuous(breaks  = c(0,1), labels=c("T/T", "T/C"))  +
						geom_text(data = RARRES2_t_test_df[which(RARRES2_t_test_df$Location == location),], aes(label=paste0("Beta = ",Beta_lm)), x = -Inf, y = Inf,  hjust = -0.1, vjust = 1.4, inherit.aes = FALSE, size = 2) +
						ggtitle(location) +
						theme(legend.position = "none",
							plot.title = element_text(hjust = 0.5)) 
	save_figs(pRARRES2_counts_eQTL_mean[[location]], paste0(outdir,"RARRES2/",location,"_counts_eQTL_mean"), width = 9, height = 12)
}

pCombined_RARRES2 <- ggarrange(plotlist = pRARRES2_counts_eQTL_mean, nrow = 1, align = "h")
save_figs(pCombined_RARRES2, paste0(outdir, "RARRES2/counts_eQTL_mean_combined"), width = 22, height = 12)


pRARRES2_counts_eQTL_mean_norm <- list()

for (location in levels(RARRES2_means$Location)){
	pRARRES2_counts_eQTL_mean_norm[[location]] <- ggplot(RARRES2_means_norm[which(RARRES2_means_norm$Location == location),], aes(x = Genotype_Numeric, y = Mean, color = Final_Assignment)) +
						geom_point(size = 1, position = position_jitter(width = 0.05)) +
						theme_classic() +
						facet_grid(Quintile ~ Time) +
						ylab("CHCHD2 Normalized Counts")+
						xlab("rs2304376 Genotype") +
						scale_color_manual(values = line_colors) +
						new_scale("color") +
						geom_smooth(method = "lm", inherit.aes = FALSE, aes(x = Genotype_Numeric, y = Mean), color = "black", size = 0.5) +
						scale_x_continuous(breaks  = c(0,1), labels=c("T/T", "T/C"))  +
						geom_text(data = RARRES2_t_test_norm_df[which(RARRES2_t_test_norm_df$Location == location),], aes(label=paste0("beta = ", Beta_lm)), x = -Inf, y = Inf,  hjust = -0.1, vjust = 1.4,inherit.aes = FALSE, size = 2.5) +
						ggtitle(location) +
						theme(legend.position = "none",
							plot.title = element_text(hjust = 0.5)) 
	save_figs(pRARRES2_counts_eQTL_mean_norm[[location]], paste0(outdir,"RARRES2/",location,"_counts_eQTL_mean_norm"), width = 6, height = 15)
}


pCombined_RARRES2_norm <- ggarrange(plotlist = pRARRES2_counts_eQTL_mean_norm, nrow = 1, align = "h")
save_figs(pCombined_RARRES2_norm, paste0(outdir, "RARRES2/counts_eQTL_mean_norm_combined"), width = 20, height = 12)


CHCHD2_cryo <- CHCHD2[grepl("Site 3", CHCHD2$Location),]
CHCHD2_cryo$Location <- gsub("Site 3_Cryopreserved", "Cryopreserved", CHCHD2_cryo$Location) %>% gsub("Site 3", "Fresh", .)


pCHCHD2_counts_cryo_eQTL <- ggplot(CHCHD2_cryo, aes(x = Genotype, y = Counts)) +
					geom_boxplot(outlier.size = 0.3) +
					theme_classic() +
					facet_grid(Location ~ Time) +
					ylab("CHCHD2 Counts") +
					xlab("rs2304376 Genotype") +
					geom_signif(comparisons = list(c("A/A", "A/G")), 
								map_signif_level=TRUE, y = max(CHCHD2[which(CHCHD2$Location != "Site 3_Cryopreserved"),]$Counts) + 5,
								test = "t.test")+
					ylim(0,max(CHCHD2[which(CHCHD2$Location != "Site 3_Cryopreserved"),]$Counts) + 25)
save_figs(pCHCHD2_counts_cryo_eQTL, paste0(outdir,"CHCHD2/counts_cryo_eQTL"), width = 5, height = 7)


CHCHD2_means_norm_cryo <- CHCHD2_means_norm[grepl("Site 3", CHCHD2_means_norm$Location),]
CHCHD2_means_norm_cryo$Location <- gsub("Site 3_Cryopreserved", "Cryopreserved", CHCHD2_means_norm_cryo$Location) %>% gsub("Site 3", "Fresh", .)
CHCHD2_means_norm_cryo$Location <- factor(CHCHD2_means_norm_cryo$Location, levels = c("Fresh", "Cryopreserved"))

t_test_norm_df_cryo <- t_test_norm_df[grepl("Site 3", t_test_norm_df$Location),]
t_test_norm_df_cryo$Location <- gsub("Site 3_Cryopreserved", "Cryopreserved", t_test_norm_df_cryo$Location) %>% gsub("Site 3", "Fresh", .)



pCHCHD2_counts_eQTL_mean_norm_cryo <- ggplot(CHCHD2_means_norm_cryo, aes(x = Genotype_Numeric, y = Mean, color = Final_Assignment)) +
					geom_point(size = 1, position = position_jitter(width = 0.05)) +
					theme_classic() +
					facet_grid(factor(Location, levels = c("Fresh", "Cryopreserved")) ~ Time) +
					ylab("CHCHD2 Normalized Counts")+
					xlab("rs2304376 Genotype") +
					scale_color_manual(values = line_colors) +
					new_scale("color") +
					geom_smooth(method = "lm", inherit.aes = FALSE, aes(x = Genotype_Numeric, y = Mean), color = "black", size = 0.5) +
					scale_x_continuous(breaks  = c(0,0.5), labels=c("A/A", "A/G"))  +
					geom_text(data = t_test_norm_df_cryo, aes(label=paste0("beta = ", Beta_lm)), x = -Inf, y = Inf,  hjust = -0.1, vjust = 1.4,inherit.aes = FALSE, size = 2.5)
save_figs(pCHCHD2_counts_eQTL_mean_norm_cryo, paste0(outdir,"CHCHD2/counts_eQTL_mean_norm_cryo"), width = 10, height = 7)






##### Pluri Genes #####
pluri_genes <- data.frame(Gene = c("MYC", "NANOG", "POU5F1", "SOX2"), ENSG = c("ENSG00000136997", "ENSG00000111704", "ENSG00000204531", "ENSG00000181449"))

df_list <- list()
for (gene in pluri_genes$Gene){
	df_list[[gene]] <- expression_df(seurat_list, pluri_genes[which(pluri_genes$Gene == gene),"ENSG"], c("Location", "Time", "Final_Assignment"))
	df_list[[gene]]$Gene <- gene
}
df <- data.table(do.call(rbind, df_list))

dir.create(paste0(outdir,"pluri_genes/"))




##### Read in significance #####
pluri_deg <- readRDS("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Expression_Boxplots/pluri_degs/LR_DEGs_4pluri_genes.rds")


pluri_deg <- lapply(names(pluri_deg), function(x){
	strings <- unlist(str_split(x, "_"))
	pluri_deg[[x]]$Location <- gsub("Brisbane", "Site 1", strings[1]) %>% gsub("Melbourne", "Site 2", .) %>% gsub("Sydney", "Site 3", .)
	pluri_deg[[x]]$Cryopreservation <- strings[2]
	pluri_deg[[x]]$Final_Assignment <- strings[3]
	return(pluri_deg[[x]])
})

pluri_deg_dt <- do.call(rbind,pluri_deg)
colnames(pluri_deg_dt) <- gsub("GeneID", "Gene", colnames(pluri_deg_dt))
pluri_deg_dt$Symbol <- "*"
pluri_deg_dt$counts_position <- ifelse(pluri_deg_dt$Gene == "MYC", 26,
									ifelse(pluri_deg_dt$Gene == "NANOG", 22, 
										ifelse(pluri_deg_dt$Gene == "POU5F1", 115, 60)))

df$Cryopreservation <- ifelse(df$Location != "3_Cryopreserved", "Fresh", "Cryopreserved")
df$Location <- gsub("_Cryopreserved", "", df$Location)


p_counts <- ggplot(df[Cryopreservation != "Cryopreserved"], aes(x = Final_Assignment, y = Counts, color = Time)) +
					geom_boxplot(aes(fill = Time), outlier.size = 0.5) +
					theme_classic() +
					scale_color_manual(values = village_colors) +
					scale_fill_manual(values = alpha(village_colors, 0.3)) +
					theme(legend.position = "none",
						axis.title.x = element_blank()) +
					geom_text(
						data = pluri_deg_dt[Cryopreservation != "Cryopreserved"],
						aes(x = Final_Assignment, y = counts_position,label = Symbol), 
						color = "black",
						size = 4) +
					facet_grid(Gene ~ Location, scales = "free_y") +
					theme(axis.text.x = element_text(angle = 45, hjust = 1))
					
save_figs(p_counts, paste0(outdir,"pluri_genes/pluri_counts"), width = 16, height = 16)



pluri_deg_dt$counts_position <- ifelse(pluri_deg_dt$Gene == "MYC", 15.5,
									ifelse(pluri_deg_dt$Gene == "NANOG", 10, 
										ifelse(pluri_deg_dt$Gene == "POU5F1", 9, 12)))

p_scaled <- ggplot(df[Location != "Cryopreserved"], aes(x = Final_Assignment, y = Normalized, color = Time)) +
					geom_boxplot(aes(fill = Time),outlier.size = 0.15, lwd=0.3) +
					theme_classic() +
					facet_grid(Gene ~ Location, scales = "free_y") +
					scale_color_manual(values = village_colors) +
					scale_fill_manual(values = alpha(village_colors, 0.4)) +
					theme(axis.title.x = element_blank()) +
					ylab("Normalized Expression") +
					theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
					geom_text(
						data = pluri_deg_dt[Cryopreservation != "Cryopreserved"],
						aes(x = Final_Assignment, y = counts_position,label = Symbol), 
						color = "black",
						size = 3) 
save_figs(p_scaled, paste0(outdir,"pluri_genes/pluri_normalized"), width = 15, height = 10)


df_cryo <- df[grepl("Site 3", df$Location),]
df_cryo$Location <- factor(gsub("Site 3_Cryopreserved", "Cryopreserved", df_cryo$Location) %>% gsub("Site 3", "Fresh", .), levels = c("Fresh", "Cryopreserved"))


p_counts_cryo <- ggplot(df_cryo, aes(x = Final_Assignment, y = Counts, color = Time)) +
					geom_boxplot(outlier.size = 0.5) +
					theme_classic() +
					facet_grid(Gene ~ Location, scales = "free_y") +
					scale_color_manual(values = village_colors) +
					theme(legend.position = "none",
						axis.title.x = element_blank()) +
					theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_figs(p_counts_cryo, paste0(outdir,"pluri_genes/pluri_counts_cryo"), width = 10, height = 8)

p_scaled_cryo <- ggplot(df_cryo, aes(x = Final_Assignment, y = Normalized, color = Time)) +
					geom_boxplot(outlier.size = 0.5) +
					theme_classic() +
					facet_grid(Gene ~ Location, scales = "free_y") +
					scale_color_manual(values = village_colors) +
					theme(legend.position = "none",
						axis.title.x = element_blank()) +
					theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
					ylab("Normalized Expression")
save_figs(p_scaled_cryo, paste0(outdir,"pluri_genes/pluri_normalized_cryo"), width = 10, height = 9)



