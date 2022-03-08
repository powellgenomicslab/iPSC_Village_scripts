##### Reason: want to test nb interaction model before running on large number of genes

library(haven)
library(ggplot2)
library(lme4)
library(glmmTMB)
library(Seurat)
library(tidyverse)
library(lmtest)



dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
outdir <- paste0(dir,"output/test_interaction_nb_model")
dir.create(outdir, recursive = TRUE)


##### Read in seurat with genes #####
seurat <- readRDS(paste0(dir,"output/Distribution_tests/seurat_integrated_all_times_clustered_1pct_expressing.rds"))

seurat@meta.data$Location <- gsub("_Baseline", "", seurat@meta.data$Location) %>% gsub("_Village.+", "", .) %>% gsub("Thawed", "Cryopreserved",.)
seurat@meta.data$Time <- gsub("Thawed Village Day 0", "Baseline", seurat@meta.data$Time) %>% gsub("Thawed Village Day 7", "Village", .)




### Make dataframes for modeling ###
gene <- "ENSG00000106153"

df <- data.frame("Expression" = data.frame(seurat[["SCT"]]@counts[gene,]), "Line" = seurat@meta.data$Final_Assignment, "Village" = ifelse(seurat@meta.data$Time == "Baseline", 0, 1), "Replicate" = gsub("[A-Z][a-z]+", "", seurat@meta.data$MULTI_ID), "Location" = seurat@meta.data$Location) ### The results are not dependent on the order or the covariates included so don't need to include extra covariates
colnames(df)[1] <- "Expression"





### Calculate R2 for each separately + remaining (1 - model with all in ) & use proporiton of these in figure
mode <- list()
model[["Base"]] <- glmmTMB(Expression ~ 1, data = df, family = nbinom2)
model[["Line"]] <- glmmTMB(Expression ~ 1 + Line , data = df, family = nbinom2)
model[["Replicate"]] <- glmmTMB(Expression ~ 1 + Replicate , data = df, family = nbinom2)
model[["Location"]] <- glmmTMB(Expression ~ 1 + Location , data = df, family = nbinom2)
model[["Village"]] <- glmmTMB(Expression ~ 1 + Village, data = df, family = nbinom2)
model[["Line_Replicate_Location_Village"]] <- glmmTMB(Expression ~ 1 + Village, data = df, family = nbinom2)
model[["Line_Replicate_Location_Village"]]



### Calculate model together ###
pseudoR2[[1-logLik(model_line)/logLik(base_model)


1-logLik(model_rep)/logLik(model)


1-logLik(model_site)/logLik(model)



1-logLik(model_line_village)/logLik(model)
1-logLik(model_village_line)/logLik(model)


### fit model with interaction - lrt
anova(model_village_line_rep_site2, model_village_line_rep_site_interaction2)
pchisq(137.96, df=1, lower.tail=FALSE)/2