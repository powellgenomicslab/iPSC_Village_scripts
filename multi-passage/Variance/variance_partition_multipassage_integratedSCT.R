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
icc_action_outdir <- paste0(args[1])
icc_outdir <- paste0(args[2])
model_action_outdir <- paste0(args[3])
model_outdir <- paste0(args[4])
resid_outdir <- paste0(args[5])
gene <- as.character(args[6])

print(icc_outdir)
print(icc_outdir)
print(model_outdir)
print(resid_outdir)
print(gene)



##### Read in seurat with genes #####
seurat <- readRDS("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/preQC/time-integrated_filtered_seurat_1pct_expressing.rds")



### Make DF for modeling ###
df_hier_unscale <- data.frame("Expression" = seurat[["SCT"]]@scale.data[gene,], "Line" = seurat@meta.data$Assignment, "Time" = as.factor(seurat@meta.data$Pool))
colnames(df_hier_unscale)[1] <- "Expression"



##### Leave one out method #####
variables <- c("Line", "Time")


model_all <- as.formula(paste0("Expression ~ (1|", paste0(variables, collapse = ") + (1|"), ")"))


boolFalse<-F
while(boolFalse==F & length(variables) > 0){
  tryCatch({
    print(variables)
    model_glmmtmb <- suppress_warnings(glmmTMB(formula = noquote(model_all), data = df_hier_unscale, REML = TRUE), "giveCsparse")
    boolFalse<-T
  },error=function(e){
    if (length(variables) > 1){
        variables <- variables[1:(length(variables) -1)]
    } else {
        variables <- c()
    }
  })
}


if (!length(variables) == 0){


    ### Deal with singular fits by removing last variable until a fit can be found - ordered in variables buy importance
    while (!model_glmmtmb$sdr$pdHess & length(variables) > 0 ){
        print("Singular fit: removing last variable and rerunning with one less covariate.")
        if (length(variables) > 1){
            variables <- variables[1:(length(variables) -1)]
            print(variables)
            model_all <- as.formula(paste0("Expression ~ (1|", paste0(variables, collapse = ") + (1|"), ")"))
            model_glmmtmb <- suppress_warnings(glmmTMB(formula = noquote(model_all), data = df_hier_unscale, REML = TRUE), "giveCsparse")
        } else {
            variables <- c()
        }
    }

    print(variables)

    if (length(variables) > 0){

        model_loo <- list()

        icc <- data.table(grp = variables, P = as.numeric(NA))

        for (variable in variables){
            print(variable)
            if (length(variables) > 1){
                model <- as.formula(paste0("Expression ~ (1|", paste0(variables[!variables %in% variable], collapse = ") + (1|"), ")"))
            } else {
                model <- as.formula(paste0("Expression ~ 1"))
            }
            model_loo[[variable]] <- suppress_warnings(glmmTMB(formula = noquote(model), data = df_hier_unscale, REML = TRUE), "giveCsparse")
            icc[grp == variable]$P <- anova(model_loo[[variable]], model_glmmtmb)$`Pr(>Chisq)`[2]
        }


        if (!(any(icc[grp != "Residual"]$P > 0.05/length(variables)) | any(is.na(icc[grp != "Residual"]$P)))){
            model_loo_updated <- model_loo

            updated_model <- as.formula(paste0("Expression ~ 1 + (1|", paste0(variables, collapse = ") + (1|"), ")"))

            model_loo_updated[["all"]] <- suppress_warnings(glmmTMB(formula = noquote(updated_model), data = df_hier_unscale, REML = TRUE), "giveCsparse")

            ### Calculate the variance explained by each of the included variables ###
            icc <- icc_glmmtmb(model_loo_updated[["all"]])


            ### Recalculate significance ###
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


        while((any(icc[grp != "Residual"]$P > 0.05/length(variables)) | any(is.na(icc[grp != "Residual"]$P)))){

            print("Removing non-significant vartiables and retesting signficance")

            ##### Identify variables to keep #####
            variables <- icc[P < 0.05/length(variables)]$grp

            if (length(variables) > 0){
                
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


        
            } else {
                icc <- data.table(grp=character(), vcov=numeric(), icc=numeric(), percent=numeric(), P=numeric(), gene=character())
                model_loo_updated <- list()
            }
        }

        action_variables <- c()

        if (length(variables) > 1){
            ### Add in actions of the significant variables
            if ("Line" %in% variables & "Time" %in% variables){
                action_variables <- c(action_variables, "Line:Time")
            }

            model_all_action <- as.formula(paste0("Expression ~ (1|", paste0(c(variables, action_variables), collapse = ") + (1|"), ")"))


            boolFalse<-F
            while(boolFalse==F & length(action_variables) > 0){
                tryCatch({
                    print(c(variables, action_variables))
                    model_glmmtmb_action <- suppress_warnings(glmmTMB(formula = noquote(model_all_action), data = df_hier_unscale, REML = TRUE), "giveCsparse")
                    boolFalse<-T
                },error=function(e){
                    if (length(action_variables) > 1){
                        action_variables <- action_variables[1:(length(action_variables) -1)]
                    } else {
                        action_variables <- c()
                    }   
                })
            }

            ### Deal with singular fits by removing last variable until a fit can be found - ordered in variables by importance
            while (!model_glmmtmb_action$sdr$pdHess & length(action_variables) > 0 ){
                print("Singular fit: removing last variable and rerunning with one less covariate.")
                if (length(action_variables) > 1){
                    action_variables <- action_variables[1:(length(action_variables) -1)]
                    print(c(action_variables, variables))
                    model_all_action <- as.formula(paste0("Expression ~ (1|", paste0(c(variables, action_variables), collapse = ") + (1|"), ")"))
                    model_glmmtmb_action <- suppress_warnings(glmmTMB(formula = noquote(model_all_action), data = df_hier_unscale, REML = TRUE), "giveCsparse")
                } else {
                    action_variables <- c()
                }
            }

            if (length(action_variables) > 0){

                model_loo_action <- list()

                icc_action <- data.table(grp = action_variables, P = as.numeric(NA))

                for (variable in c(action_variables)){
                    print(variable)
                    if (length(action_variables) > 1){
                        model_action <- as.formula(paste0("Expression ~ (1|", paste0(c(variables, action_variables)[!c(variables, action_variables) %in% variable], collapse = ") + (1|"), ")"))
                    } else {
                        model_action <- as.formula(paste0("Expression ~ (1|", paste0(variables[!variables %in% variable], collapse = ") + (1|"), ")"))
                    }
                    model_loo_action[[variable]] <- suppress_warnings(glmmTMB(formula = noquote(model_action), data = df_hier_unscale, REML = TRUE), "giveCsparse")
                    icc_action[grp == variable]$P <- anova(model_loo_action[[variable]], model_glmmtmb_action)$`Pr(>Chisq)`[2]
                }


                if (!(any(icc_action[grp != "Residual"]$P > 0.05/length(c(variables, action_variables))) | any(is.na(icc_action[grp != "Residual"]$P)))){
                    model_loo_action_updated <- model_loo_action

                    updated_model_action <- as.formula(paste0("Expression ~ 1 + (1|", paste0(c(variables, action_variables), collapse = ") + (1|"), ")"))

                    model_loo_action_updated[["all"]] <- suppress_warnings(glmmTMB(formula = noquote(updated_model_action), data = df_hier_unscale, REML = TRUE), "giveCsparse")

                    ### Calculate the variance explained by each of the included variables ###
                    icc_action <- icc_glmmtmb(model_loo_action_updated[["all"]])


                    ### Recalculate significance ###
                    icc_action$P <- as.numeric(NA)
                    icc_action$gene <- gene

                    for (variable in c(variables, action_variables)){
                        print(variable)
                        if (length(c(action_variables)) > 1){
                            model <- as.formula(paste0("Expression ~ 1 + (1|", paste0(c(variables, action_variables)[!c(variables, action_variables) %in% variable], collapse = ") + (1|"), ")"))
                        } else {
                            model <- as.formula(paste0("Expression ~ (1|", paste0(variables[!variables %in% variable], collapse = ") + (1|"), ")"))
                        }
                        model_loo_action_updated[[variable]] <- suppress_warnings(glmmTMB(formula = noquote(model), data = df_hier_unscale, REML = TRUE), "giveCsparse")
                        icc_action[grp == variable]$P <- anova(model_loo_action_updated[[variable]], model_loo_action_updated[["all"]])$`Pr(>Chisq)`[2]
                    }
                }


                while((any(icc_action[!(grp %in% c("Residual", variables))]$P > 0.05/length(c(variables, action_variables))) | any(is.na(icc_action[!(grp %in%c("Residual", variables))]$P)))){

                    print("Removing non-significant vartiables and retesting signficance")

                    ##### Identify variables to keep #####
                    action_variables <- icc_action[!(grp %in% c("Residual", variables)) & P < 0.05/length(c(variables, action_variables))]$grp

                    if (length(action_variables) > 0){
                        
                        ##### Calculate full model #####
                        updated_model_action <- as.formula(paste0("Expression ~ 1 + (1|", paste0(c(variables, action_variables), collapse = ") + (1|"), ")"))


                        model_loo_action_updated <- list()
                        model_loo_action_updated[["all"]] <- suppress_warnings(glmmTMB(formula = noquote(updated_model_action), data = df_hier_unscale, REML = TRUE), "giveCsparse")



                        ### Calculate the variance explained by each of the included variables ###
                        icc_action <- icc_glmmtmb(model_loo_action_updated[["all"]])



                        ### Recalfulate significance ###
                        icc_action$P <- as.numeric(NA)
                        icc_action$gene <- gene

                        for (variable in c(variables, action_variables)){
                            print(variable)
                            model_action <- as.formula(paste0("Expression ~ 1 + (1|", paste0(c(variables, action_variables)[!(c(variables, action_variables) %in% variable)], collapse = ") + (1|"), ")"))
                            model_loo_action_updated[[variable]] <- suppress_warnings(glmmTMB(formula = noquote(model_action), data = df_hier_unscale, REML = TRUE), "giveCsparse")
                            icc_action[grp == variable]$P <- anova(model_loo_action_updated[[variable]], model_loo_action_updated[["all"]])$`Pr(>Chisq)`[2]
                        }
                    } else {
                        icc_action <- data.table(grp=character(), vcov=numeric(), icc=numeric(), percent=numeric(), P=numeric(), gene=character())
                        model_loo_action_updated <- list()
                    }
                }
            } else {
                icc_action <- data.table(grp=character(), vcov=numeric(), icc=numeric(), percent=numeric(), P=numeric(), gene=character())
                model_loo_action_updated <- list()
            }
            if (nrow(icc_action) > nrow(icc)){
                saveRDS(icc_action, paste0(icc_action_outdir, gene, "_icc.rds"), compress = TRUE)
                saveRDS(model_loo_action_updated, paste0(model_action_outdir, gene, "_fitted_models.rds"), compress = TRUE)
            }


            ### If line is significant, then get residuals for downstream qtl checks ###
            if ("Line" %in% variables){
                print("Making residuals for qtl detection")
                if (length(variables) > 1){
                        if (length(action_variables) > 0){
                            model_no_line <- as.formula(paste0("Expression ~ (1|", paste0(c(variables, action_variables)[-grep("Line", c(variables, action_variables))], collapse = ") + (1|"), ")"))
                        } else {
                            model_no_line <- as.formula(paste0("Expression ~ (1|", paste0(variables[!variables %in% "Line"], collapse = ") + (1|"), ")"))
                        }
                } else {
                    model_no_line <- as.formula(paste0("Expression ~ 1"))
                }
                fit_no_line <- glmmTMB(formula = noquote(model_no_line), data = df_hier_unscale, REML = TRUE)
                residuals <- resid(fit_no_line)
                saveRDS(residuals, paste0(resid_outdir, gene, "_residuals4qtl.rds"), compress = TRUE)
            }
        
        } else if (length(variables) ==0) {
            icc <- data.table(grp=character(), vcov=numeric(), icc=numeric(), percent=numeric(), P=numeric(), gene=character())
            model_loo_updated <- list()
        }
    } else {
        icc <- data.table(grp=character(), vcov=numeric(), icc=numeric(), percent=numeric(), P=numeric(), gene=character())
        model_loo_updated <- list()
    }
}


saveRDS(icc, paste0(icc_outdir, gene, "_icc.rds"), compress = TRUE)
saveRDS(model_loo_updated, paste0(model_outdir, gene, "_fitted_models.rds"), compress = TRUE)


fin()