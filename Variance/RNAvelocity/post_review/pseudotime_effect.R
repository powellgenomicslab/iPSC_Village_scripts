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
plot_outdir <- paste0(args[3])
gene <- as.character(args[4])
effects_outdir <- paste0(args[5])


print(icc_interaction_outdir)
print(icc_outdir)
print(plot_outdir)
print(gene)

line_colors <- c(FSA0006 = "#F79E29", MBE1006 = "#9B2C99", TOB0421 = "#35369C")



##### Read in data #####
### Seurat object with normalized data and covariates needed ###
seurat <- readRDS("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/RNAvelocity/post_review/data/seurat_integrated_all_times_clustered_1pct_expressing_pseudotime.rds")

### Dataframe of icc summaries so know what variables need to be fit for each gene ###
icc_summary <- fread("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/variance_partitioning_all_cells/combined/sig_results.tsv.gz", sep = "\t")
colnames(icc_summary) <- gsub("gene", "ensg", colnames(icc_summary))


### Make DF for modeling ###
df_hier_unscale <- data.frame("Expression" = seurat[["SCT"]]@scale.data[gene,],"Normalized Counts" = seurat[["SCT"]]@counts[gene,], "Log Expression" = seurat[["SCT"]]@data[gene,],  "Village" = as.factor(ifelse(seurat@meta.data$Time == "Baseline", 0, 1)), "Line" = seurat@meta.data$Final_Assignment, "Replicate" = as.factor(gsub("[A-Z][a-z]+", "", seurat@meta.data$MULTI_ID)),  "Cryopreserved" = seurat$Cryopreserved, "Site" = seurat$Location, "Pseudotime" = round(seurat$latent_time, 2))
colnames(df_hier_unscale)[1] <- "Expression"



##### Get list of variables to fit before testing pseudotime effect #####
variables <- c(icc_summary[ensg == gene & grp != "Residual"]$grp, "Pseudotime")
variables <- c(variables[!grepl(":",variables)], variables[grepl(":", variables)])

model_all <- as.formula(paste0("Expression ~ (1|", paste0(variables, collapse = ") + (1|"), ")"))

if (length(variables) > 1){
    model <- as.formula(paste0("Expression ~ (1|", paste0(variables[!variables %in% "Pseudotime"], collapse = ") + (1|"), ")"))
} else {
    model <- as.formula(paste0("Expression ~ 1"))
}


boolFalse<-F
while(boolFalse==F & length(variables) > 0){
  tryCatch({
    print(variables)
    model_glmmtmb <- suppress_warnings(glmmTMB(formula = noquote(model_all), data = df_hier_unscale, REML = TRUE), "giveCsparse")
    model_glmmtmb_loo <- suppress_warnings(glmmTMB(formula = noquote(model), data = df_hier_unscale, REML = TRUE), "giveCsparse")
    boolFalse<-T
  },error=function(e){
    if (length(variables) > 1){
        variables <- variables[1:(length(variables) -1)]
    } else {
        variables <- c()
    }
  })
}


if ("Pseudotime" %in% variables){

    ### Deal with singular fits by removing last variable until a fit can be found - ordered in variables buy importance
    while ((!model_glmmtmb$sdr$pdHess | !model_glmmtmb_loo$sdr$pdHess) & length(variables) > 0){
        print("Singular fit: removing last variable and rerunning with one less covariate.")
        if (length(variables) > 1){
            variables <- variables[1:(length(variables) -1)]
            print(variables)
            model_all <- as.formula(paste0("Expression ~ (1|", paste0(variables, collapse = ") + (1|"), ")"))
            model_glmmtmb <- suppress_warnings(glmmTMB(formula = noquote(model_all), data = df_hier_unscale, REML = TRUE), "giveCsparse")
            model <- as.formula(paste0("Expression ~ (1|", paste0(variables[!variables %in% "Pseudotime"], collapse = ") + (1|"), ")"))
            model_glmmtmb_loo <- suppress_warnings(glmmTMB(formula = noquote(model), data = df_hier_unscale, REML = TRUE), "giveCsparse")
        } else {
            variables <- c()
        }
    }

    print(variables)
    
    if ("Pseudotime" %in% variables){

        icc <- data.table(grp = c("Pseudotime", "Residual"), P = as.numeric(NA))

        icc[grp == "Pseudotime"]$P <- anova(model_glmmtmb_loo, model_glmmtmb)$`Pr(>Chisq)`[2]


        if (icc[grp == "Pseudotime"]$P < 0.05/length(variables)){

            ##### Calculate full model #####
            updated_model <- as.formula(paste0("Expression ~ 1 + (1|", paste0(variables, collapse = ") + (1|"), ")"))


            model_loo_updated <- list()
            model_loo_updated[["all"]] <- suppress_warnings(glmmTMB(formula = noquote(updated_model), data = df_hier_unscale, REML = TRUE), "giveCsparse")

            ### Calculate the variance explained by each of the included variables ###
            icc <- icc_glmmtmb(model_loo_updated[["all"]])

            ### Recalfulate significance ###
            icc$P <- as.numeric(NA)
            icc$gene <- gene

            for (variable in variables){
                print(variable)
                if (length(variables) > 1){
                    model <- as.formula(paste0("Expression ~ 1 + (1|", paste0(variables[!variables %in% variable], collapse = ") + (1|"), ")"))
                } else {
                    model <- as.formula(paste0("Expression ~ 1"))
                }
                model_loo_updated[[variable]] <- suppress_warnings(glmmTMB(formula = noquote(model), data = df_hier_unscale, REML = TRUE), "giveCsparse")
                icc[grp == variable]$P <- anova(model_loo_updated[[variable]], model_loo_updated[["all"]])$`Pr(>Chisq)`[2]
            }
        }

        interaction_variables <- c()

        if ("Line" %in% variables & "Pseudotime" %in% variables){
            interaction_variables <- c(interaction_variables, "Line:Pseudotime")

            model_all_interaction <- as.formula(paste0("Expression ~ (1|", paste0(c(variables, interaction_variables), collapse = ") + (1|"), ")"))

            boolFalse<-F
            while(boolFalse==F & length(interaction_variables) > 0){
            tryCatch({
                print(c(variables, interaction_variables))
                model_glmmtmb_interaction <- suppress_warnings(glmmTMB(formula = noquote(model_all_interaction), data = df_hier_unscale, REML = TRUE), "giveCsparse")
                boolFalse<-T
            },error=function(e){
                if (length(interaction_variables) > 1){
                    interaction_variables <- interaction_variables[1:(length(interaction_variables) -1)]
                } else {
                    interaction_variables <- c()
                }
            })
            }

            if ("Line:Pseudotime" %in% interaction_variables){
                ### Deal with singular fits by removing last variable until a fit can be found - ordered in variables buy importance
                while (!model_glmmtmb_interaction$sdr$pdHess & length(interaction_variables) > 0 ){
                    print("Singular fit: removing last variable and rerunning with one less covariate.")
                    if (length(interaction_variables) > 1){
                        interaction_variables <- interaction_variables[1:(length(interaction_variables) -1)]
                        print(c(interaction_variables, variables))
                        model_all_interaction <- as.formula(paste0("Expression ~ (1|", paste0(c(variables, interaction_variables), collapse = ") + (1|"), ")"))
                        model_glmmtmb_interaction <- suppress_warnings(glmmTMB(formula = noquote(model_all_interaction), data = df_hier_unscale, REML = TRUE), "giveCsparse")
                    } else {
                        interaction_variables <- c()
                    }
                }

                if ("Line:Pseudotime" %in% interaction_variables){
                    
                    model_loo_interaction <- list()

                    icc_interaction <- data.table(grp = interaction_variables, P = as.numeric(NA))

                    for (variable in c(interaction_variables)){
                        print(variable)
                        if (length(interaction_variables) > 1){
                            model_interaction <- as.formula(paste0("Expression ~ (1|", paste0(c(variables, interaction_variables)[!c(variables, interaction_variables) %in% variable], collapse = ") + (1|"), ")"))
                        } else {
                            model_interaction <- as.formula(paste0("Expression ~ (1|", paste0(variables[!variables %in% variable], collapse = ") + (1|"), ")"))
                        }
                        model_loo_interaction[[variable]] <- suppress_warnings(glmmTMB(formula = noquote(model_interaction), data = df_hier_unscale, REML = TRUE), "giveCsparse")
                        icc_interaction[grp == variable]$P <- anova(model_loo_interaction[[variable]], model_glmmtmb_interaction)$`Pr(>Chisq)`[2]
                    }

                    if (icc_interaction[grp == "Line:Pseudotime"]$P < 0.05/length(c(variables, interaction_variables)) & !is.na(icc_interaction[grp == "Line:Pseudotime"]$P)){
                        print("interaction significant:")
                        print(icc_interaction)

                        model_loo_interaction_updated <- model_loo_interaction

                        updated_model_interaction <- as.formula(paste0("Expression ~ 1 + (1|", paste0(c(variables, interaction_variables), collapse = ") + (1|"), ")"))

                        model_loo_interaction_updated[["all"]] <- suppress_warnings(glmmTMB(formula = noquote(updated_model_interaction), data = df_hier_unscale, REML = TRUE), "giveCsparse")

                        ### Calculate the variance explained by each of the included variables ###
                        icc_interaction <- icc_glmmtmb(model_loo_interaction_updated[["all"]])


                        ### Recalculate significance ###
                        icc_interaction$P <- as.numeric(NA)
                        icc_interaction$gene <- gene

                        for (variable in c(variables, interaction_variables)){
                            print(variable)
                            # if (length(c(interaction_variables)) > 1){
                                model <- as.formula(paste0("Expression ~ 1 + (1|", paste0(c(variables, interaction_variables)[!c(variables, interaction_variables) %in% variable], collapse = ") + (1|"), ")"))
                            # } else {
                            #     model <- as.formula(paste0("Expression ~ (1|", paste0(variables[!variables %in% variable], collapse = ") + (1|"), ")"))
                            # }
                            model_loo_interaction_updated[[variable]] <- suppress_warnings(glmmTMB(formula = noquote(model), data = df_hier_unscale, REML = TRUE), "giveCsparse")
                            icc_interaction[grp == variable]$P <- anova(model_loo_interaction_updated[[variable]], model_loo_interaction_updated[["all"]])$`Pr(>Chisq)`[2]
                        }


                        plot <- ggplot(df_hier_unscale, aes(Pseudotime, Expression, color = Line)) +
                                    geom_point(alpha = 0.3)  +
                                    # geom_smooth(method = "lm", se = TRUE) +
                                    geom_smooth(se = TRUE) +
                                    theme_classic() +
                                    scale_color_manual(values = line_colors)

                        ggsave(plot, filename = paste0(plot_outdir, gene,".png"), width = 5, height = 4)

                        plot_facet <- ggplot(df_hier_unscale, aes(Pseudotime, Expression, color = Line)) +
                                    geom_point(alpha = 0.3)  +
                                    facet_wrap(vars(paste0(Site, " ", Cryopreserved)), ncol = 1, scales = "free_y") +
                                    geom_smooth(se = TRUE) +
                                    # geom_smooth(method = "lm", se = TRUE) +
                                    theme_classic() +
                                    scale_color_manual(values = line_colors)

                        ggsave(plot_facet, filename = paste0(plot_outdir, gene,"_site_facet.png"), width = 4, height = 7)


                        saveRDS(icc_interaction, paste0(icc_interaction_outdir, gene, "_icc.rds"), compress = TRUE)


                        model_interaction_effects_formula <- as.formula(paste0("Expression ~ ", paste0(c(variables, interaction_variables), collapse = " + ")))
                        model_interaction_effects <- suppress_warnings(glmmTMB(formula = noquote(model_interaction_effects_formula), data = df_hier_unscale, REML = TRUE), "giveCsparse")

                        dt <- data.table(grp = rownames(coef(summary(model_interaction_effects))$cond), Effect = coef(summary(model_interaction_effects))$cond[,"Estimate"], P = coef(summary(model_interaction_effects))$cond[,"Pr(>|z|)"])

                        saveRDS(dt, paste0(effects_outdir, gene, "_effects.rds"), compress = TRUE)
                    }
                } 
            } 
        } 
    } else {
        icc <- data.table(grp=character(), vcov=numeric(), icc=numeric(), percent=numeric(), P=numeric(), gene=character())
    }
}


saveRDS(icc, paste0(icc_outdir, gene, "_icc.rds"), compress = TRUE)


