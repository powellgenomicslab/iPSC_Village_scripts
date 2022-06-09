library(haven)
library(ggplot2)
library(glmmTMB)
library(Seurat)
library(tidyverse)
library(specr)
library(data.table)
library(dsLib)
library(pkgcond)
library(texreg)


inicio("Starting Analysis")



##### Define functions #####
icc_glmmtmb <- function(model, percent = TRUE) {
    tmp <- VarCorr(model)
    var <- do.call(rbind, lapply(names(tmp$cond), function(x) data.table("grp" = x, "vcov" = attr(tmp$cond[[x]], "stddev")^2)))
    var <- rbind(var, data.table("grp" = "Residual", "vcov" = sigma(model)^2))
    sum_var <- sum(var$vcov)
    var <- var %>% dplyr::mutate(icc = vcov/sum_var)
    if (isTRUE(percent)) {
        var <- var %>% dplyr::mutate(percent = .data$icc * 100)
    }
    return(var)
}



##### Bring in variables #####
### Bring in arguments
args <- commandArgs(trailingOnly = TRUE)
icc_interaction_outdir <- paste0(args[1])
icc_outdir <- paste0(args[2])
model_interaction_outdir <- paste0(args[3])
model_outdir <- paste0(args[4])
resid_outdir <- paste0(args[5])
gene <- as.character(args[6])

print(icc_outdir)
print(icc_outdir)
print(model_outdir)
print(resid_outdir)
print(gene)



##### Read in data #####
### Seurat object with normalized data and covariates needed ###
seurat <- readRDS("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/RNAvelocity/post_review/data/seurat_integrated_all_times_clustered_1pct_expressing_pseudotime.rds")

### Dataframe of icc summaries so know what variables need to be fit for each gene ###
icc_summary <- fread("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/variance_partitioning_all_cells/combined/sig_results.tsv.gz", sep = "\t")
colnames(icc_summary) <- gsub("gene", "ensg", colnames(icc_summary))


### Make DF for modeling ###
df_hier_unscale <- data.frame("Expression" = seurat[["SCT"]]@scale.data[gene,], "Village" = as.factor(ifelse(seurat@meta.data$Time == "Baseline", 0, 1)), "Line" = seurat@meta.data$Final_Assignment, "Replicate" = as.factor(gsub("[A-Z][a-z]+", "", seurat@meta.data$MULTI_ID)),  "Cryopreserved" = seurat$Cryopreserved, "Site" = seurat$Location, "Pseudotime" = round(seurat$latent_time, 2))
colnames(df_hier_unscale)[1] <- "Expression"



##### Get list of variables to fit before testing pseudotime effect #####
variables <- icc_summary[ensg == gene & grp != "Residual"]$grp



##### Test pseudotime impact #####
### Fit the known variables to get residuals ###
model_all <- as.formula(paste0("Expression ~ (1|", paste0(variables, collapse = ") + (1|"), ")"))
model_glmmtmb <- suppress_warnings(glmmTMB(formula = noquote(model_all), data = df_hier_unscale, REML = TRUE), "giveCsparse")


### Test pseudotime on residuals ###
df_hier_unscale$Residuals <- resid(model_glmmtmb)

model_pseudotime <- suppress_warnings(glmmTMB(Residuals ~ Pseudotime, data = df_hier_unscale, REML = TRUE), "giveCsparse")
model_pseudotime2 <- suppress_warnings(glmmTMB(Residuals ~ 1, data = df_hier_unscale, REML = TRUE), "giveCsparse")

### Test with Anova ###
P_value <- anova(model_pseudotime2, model_pseudotime)$`Pr(>Chisq)`[2]
P_value <- anova(model_glmmtmb, model_pseudotime2)$`Pr(>Chisq)`[2]


test_plot <- ggplot(df_hier_unscale, aes(Pseudotime, Residuals, color = Site)) +
    geom_point()  +
    facet_grid(~Line) +
    geom_smooth(method = "lm", se = FALSE)

ggsave(test_plot, filename = "/directflow/SCCGGroupShare/projects/DrewNeavin/test.png")


if (P_value < 0.05/(length(variables) + 1)){
    ### Test for amount of variance explained ###
    model_pseudotime2 <- suppress_warnings(glmmTMB(Residuals ~ 1 + Pseudotime, data = df_hier_unscale, REML = TRUE), "giveCsparse")



    ### Test for interactions

}
