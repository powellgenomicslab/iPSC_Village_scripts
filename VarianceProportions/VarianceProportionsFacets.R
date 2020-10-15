library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(jcolors)
library(ggridges)

##### Set up directories #####
dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
var_dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/VarianceProportions/"
outdir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/VarianceProportions/Compare_SeparatedSite/"

dir.create(outdir)

##### Read in files #####
variance_files <- list.files(paste0(var_dir,"scaterVariance_SeparatedSite"), pattern = ".txt")

variance <- lapply(variance_files, function(x){
    read.table(paste0(var_dir,"scaterVariance_SeparatedSite/",x), sep= "\t")
})
names(variance) <- gsub("_VarianceContributionsResults.txt","", variance_files)

variance_MtRbReg <- lapply(variance_files, function(x){
    read.table(paste0(var_dir,"scaterVariance_SeparatedSite_MtRb_regressed/",x), sep= "\t")
})
names(variance_MtRbReg) <- gsub("_VarianceContributionsResults.txt","", variance_files)

##### Add in a column to indicate the input #####
variance <- lapply(names(variance), function(x){
    variance[[x]]$Model <- paste0(x,"_No_Regression")
    variance[[x]]$Gene <- rownames(variance[[x]])
    return(variance[[x]])
})
names(variance) <- gsub("_VarianceContributionsResults.txt","", variance_files)

variance_MtRbReg <- lapply(names(variance_MtRbReg), function(x){
    variance_MtRbReg[[x]]$Model <- paste0(x,"_MtRb_Regression")
    variance_MtRbReg[[x]]$Gene <- rownames(variance_MtRbReg[[x]])
    return(variance_MtRbReg[[x]])
})
names(variance_MtRbReg) <- gsub("_VarianceContributionsResults.txt","", variance_files)


##### Join all the dataframes together for 
variance_df <- do.call(rbind, variance)
variance_MtRbReg_df <- do.call(rbind, variance_MtRbReg)

variance_combined <- rbind(variance_df, variance_MtRbReg_df)
colnames(variance_combined) <- gsub("FinalAssignment","Individual",colnames(variance_combined)) %>% gsub("MULTI_ID", "Pool",.)
variance_combined_long <- pivot_longer(variance_combined, cols = c("Individual","Time","Pool","percent.mt","percent.rb"),names_to = "Variable", values_to = "Variance")
variance_combined_long$Model <- factor(variance_combined_long$Model, levels = c("Brisbane_No_Regression","Sydney_No_Regression","Melbourne_No_Regression","Brisbane_MtRb_Regression","Sydney_MtRb_Regression","Melbourne_MtRb_Regression"))

plot <- ggplot(variance_combined_long, aes(log10(Variance), color = Variable)) +
    geom_line(stat="density", size=1, alpha = 0.5) +
    facet_wrap(vars(Model)) +
    theme_classic() +
    scale_color_jcolors() +
    scale_x_continuous(name="% Variance Explained", labels = c(0.001, 0.01, 0.1, 1, 10, 100), limits=c(log10(0.001),log10(100)))
ggsave(plot, filename = paste0(outdir,"facet_density_variance.png"))


########## Make a new figure with the replicate pairs as variance ##########
##### Read in files #####
variance_sample_files_Day0 <- list.files(paste0(var_dir,"scaterVariance_SeparatedReplicate"), pattern = "Day_0.+\\.txt")
variance_sample_files_Day4 <- list.files(paste0(var_dir,"scaterVariance_SeparatedReplicate"), pattern = "Day_4.+\\.txt")

variance_sample_Day0 <- lapply(variance_sample_files_Day0, function(x){
    read.table(paste0(var_dir,"scaterVariance_SeparatedReplicate/",x), sep= "\t")
})
names(variance_sample_Day0) <- gsub("_VarianceContributionsResults.txt","", variance_sample_files_Day0)

variance_sample_Day4 <- lapply(variance_sample_files_Day4, function(x){
    read.table(paste0(var_dir,"scaterVariance_SeparatedReplicate/",x), sep= "\t")
})
names(variance_sample_Day4) <- gsub("_VarianceContributionsResults.txt","", variance_sample_files_Day4)


##### Add in a column to indicate the input #####
variance_sample_Day0 <- lapply(names(variance_sample_Day0), function(x){
    variance_sample_Day0[[x]]$Sample <- x
    variance_sample_Day0[[x]] <- separate(variance_sample_Day0[[x]], col = "Sample", into = c("Replicate","Garbage","Day"), sep = "_", remove = FALSE)
    variance_sample_Day0[[x]]$Garbage <- NULL
    variance_sample_Day0[[x]]$Site  <- variance_sample_Day0[[x]]$Replicate
    variance_sample_Day0[[x]]$Site <- gsub("\\d","", variance_sample_Day0[[x]]$Site)
    variance_sample_Day0[[x]]$Day <- paste0("Day_",variance_sample_Day0[[x]]$Day)
    variance_sample_Day0[[x]]$Replicate  <- variance_sample_Day0[[x]]$Replicate
    variance_sample_Day0[[x]]$Replicate <- gsub("Sydney","Replicate_", variance_sample_Day0[[x]]$Replicate) %>% gsub("Melbourne","Replicate_", .) %>% gsub("Brisbane","Replicate_",.)
    variance_sample_Day0[[x]]$Gene <- rownames(variance_sample_Day0[[x]])
    return(variance_sample_Day0[[x]])
})
names(variance_sample_Day0) <- gsub("_VarianceContributionsResults.txt","", variance_sample_files_Day0) %>% gsub("_Day_\\d","",.)

variance_sample_Day4 <- lapply(names(variance_sample_Day4), function(x){
    variance_sample_Day4[[x]]$Sample <- x
    variance_sample_Day4[[x]] <- separate(variance_sample_Day4[[x]], col = "Sample", into = c("Replicate","Garbage","Day"), by = "_", remove = FALSE)
    variance_sample_Day4[[x]]$Garbage <- NULL
    variance_sample_Day4[[x]]$Site  <- variance_sample_Day4[[x]]$Replicate
    variance_sample_Day4[[x]]$Site <- gsub("\\d","", variance_sample_Day4[[x]]$Site)
    variance_sample_Day4[[x]]$Day <- paste0("Day_",variance_sample_Day4[[x]]$Day)
    variance_sample_Day4[[x]]$Replicate  <- variance_sample_Day4[[x]]$Replicate
    variance_sample_Day4[[x]]$Replicate <- gsub("Sydney","Replicate_", variance_sample_Day4[[x]]$Replicate) %>% gsub("Melbourne","Replicate_", .) %>% gsub("Brisbane","Replicate_",.)
    variance_sample_Day4[[x]]$Gene <- rownames(variance_sample_Day4[[x]])
    return(variance_sample_Day4[[x]])
})
names(variance_sample_Day4) <- gsub("_VarianceContributionsResults.txt","", variance_sample_files_Day4)%>% gsub("_Day_\\d","",.)


##### Order the dataframes for the timepoints to all be the same for rows and columns #####
gene_list <- lapply(names(variance_sample_Day0), function(x){
    gene1 <- as.data.frame(rownames(variance_sample_Day0[[x]]))
    colnames(gene1) <- c("Gene")
    gene2 <- as.data.frame(rownames(variance_sample_Day4[[x]]))
    colnames(gene2) <- c("Gene")
    genes <- inner_join(gene1, gene2)
    return(genes)
})
names(gene_list) <- names(variance_sample_Day0)

variance_sample_Day4 <- lapply(names(variance_sample_Day4), function(x){
    left_join(gene_list[[x]], variance_sample_Day4[[x]], by = c("Gene"))
})
names(variance_sample_Day4) <- gsub("_VarianceContributionsResults.txt","", variance_sample_files_Day4)%>% gsub("_Day_\\d","",.)

variance_sample_Day0 <- lapply(names(variance_sample_Day0), function(x){
    left_join(gene_list[[x]], variance_sample_Day0[[x]], by = c("Gene"))
})
names(variance_sample_Day0) <- gsub("_VarianceContributionsResults.txt","", variance_sample_files_Day0)%>% gsub("_Day_\\d","",.)

##### Check that all the row names match #####
lapply(names(variance_sample_Day0), function(x){
    print(all(rownames(variance_sample_Day0[[x]]) == rownames(variance_sample_Day4[[x]])))
})
### All TRUE

### Remove the Time column (all NAs since they are separated by time) ###
variance_sample_Day4 <- lapply(variance_sample_Day4, function(x){
    x$Time <- NULL
    return(x)
})

variance_sample_Day0 <- lapply(variance_sample_Day0, function(x){
    x$Time <- NULL
    return(x)
})

### Rearange the columns to be in the same order ###
variance_sample_Day4 <- lapply(names(variance_sample_Day4), function(x){
   variance_sample_Day4[[x]] <- variance_sample_Day4[[x]][,colnames(variance_sample_Day0[[x]])]
   return(variance_sample_Day4[[x]])
})
names(variance_sample_Day4) <- gsub("_VarianceContributionsResults.txt","", variance_sample_files_Day4)%>% gsub("_Day_\\d","",.)

### Subtract the 4 week from the baseline ###
variance_sample_diff <- lapply(names(variance_sample_Day4), function(x){
    M <- merge(variance_sample_Day0[[x]][,1:7],variance_sample_Day4[[x]][,1:7],by="Gene")
    S <- M[,grepl("*\\.x$",names(M))] - M[,grepl("*\\.y$",names(M))]
    df <- cbind(M[,1,drop=FALSE],S)
    return(df)
})
names(variance_sample_diff) <- names(variance_sample_Day4)

variance_sample_diff <- lapply(names(variance_sample_diff), function(x){
    variance_sample_diff[[x]]$Replicate <- x
    return(variance_sample_diff[[x]])
})
names(variance_sample_diff) <- names(variance_sample_Day4)

variance_sample_diff_df <- do.call(rbind, variance_sample_diff)
colnames(variance_sample_diff_df) <- gsub("\\.x","",colnames(variance_sample_diff_df))

variance_sample_diff_df_long <- pivot_longer(variance_sample_diff_df, cols = c("FinalAssignment",  "percent.rb",   "nCount_RNA", "nFeature_RNA", "percent.mt",   "nCount_HTO"), names_to = "Variable", values_to = "Variance_Difference")

### ggplot ###
var_dif_plot <- ggplot(variance_sample_diff_df_long, aes(Variable, Variance_Difference, fill = Variable)) +
    geom_violin() +
    geom_boxplot(width=0.1) +
    facet_wrap(vars(Replicate))
ggsave(var_dif_plot, filename = paste0(outdir,"violin_variance_diff.png"))

var_dif_ridge_plot <- ggplot(variance_sample_diff_df_long, aes(Variance_Difference, Variable, fill = Variable)) +
    geom_density_ridges() +
    facet_wrap(vars(Replicate))
ggsave(var_dif_ridge_plot, filename = paste0(outdir,"ridge_variance_diff.png"))





















##### Pivot longer for ggplot #####
variance_sample_df <- do.call(rbind, variance_paired)
colnames(variance_sample_df) <- gsub("FinalAssignment","Individual", colnames(variance_sample_df))
variance_sample_df_long <- pivot_longer(variance_sample_df, cols = c("Individual","Time","percent.mt","percent.rb","nCount_HTO","nCount_RNA","nFeature_RNA","nFeature_HTO"),names_to = "Variable", values_to = "Variance")
variance_sample_df_long$Site <- factor(variance_sample_df_long$Site, levels = c("Brisbane","Sydney","Melbourne"))
variance_sample_df_long$Replicate <- factor(variance_sample_df_long$Replicate, levels = c("1","2","3"))

colors <- c(jcolors("default"),"tan4")
names(colors)[6] <- "Brown"

plot <- ggplot(variance_sample_df_long, aes(log10(Variance), color = Variable)) +
    geom_line(stat="density", size=1, alpha = 0.5) +
    facet_grid(cols = vars(Replicate), rows = vars(Site)) +
    theme_classic() +
    scale_color_manual(values = c("#29BF12", "#00A5CF", "#DE1A1A", "#574AE2", "#FFBF00", "tan4")) +
    scale_x_continuous(name="% Variance Explained", labels = c(0.001, 0.01, 0.1, 1, 10, 100), limits=c(log10(0.001),log10(100)))
ggsave(plot, filename = paste0(outdir,"facet_density_variance_paired.png"))



########## Calculate the difference in the variance at baseline and 4 days ##########
##### Read in files #####
variance_paired <- lapply(names(variance_paired), function(x){
    variance_paired[[x]]$Site <- x
    variance_paired[[x]]$Site <- gsub("\\d","", variance_paired[[x]]$Site)
    variance_paired[[x]]$Replicate <- x
    variance_paired[[x]]$Replicate <- gsub("Sydney","", variance_paired[[x]]$Replicate) %>% gsub("Melbourne","", .) %>% gsub("Brisbane","",.)
    variance_paired[[x]]$Gene <- rownames(variance_paired[[x]])
    return(variance_paired[[x]])
})
names(variance_paired) <- gsub("_VarianceContributionsResults.txt","", variance_paired_files)







