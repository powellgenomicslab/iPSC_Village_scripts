##### Reason: want to test nb interaction model before running on large number of genes

library(haven)
library(ggplot2)
library(lme4)
library(glmmTMB)
library(Seurat)
library(tidyverse)
library(lmtest)
library(hier.part)
library(MuMIn)
library(specr)
library(lmerTest)
library(tictoc)
library(Rfast)
library(sjstats)



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


# ### Fit models
# model <- glmmTMB(Expression ~ 1 + (1|Line) + (1|Village) + (1|Replicate) + (1|Location), data = df, family = nbinom2) ### The results are not dependent on the order or the covariates included so don't need to include extra covariates
# model_small <- glmmTMB(Expression ~ 1 + (1|Line) + (1|Village) + (1|Location), data = df, family = nbinom2) ### The results are not dependent on the order or the covariates included so don't need to include extra covariates
# model_reorder <- glmmTMB(Expression ~ 1 + (1|Village) + (1|Replicate) + (1|Location) + (1|Line) , data = df, family = nbinom2) ### The results are not dependent on the order or the covariates included so don't need to include extra covariates
# model_nonrandom <- glmmTMB(Expression ~ 1 + Line + Village + Replicate + Location, data = df, family = nbinom2) ### The results are not dependent on the order or the covariates included so don't need to include extra covariates
# model_nonrandom_nested <- glmmTMB(Expression ~ 1 + Village + Replicate + Location, data = df, family = nbinom2)
# model_nonrandom_no_rep_nested <- glmmTMB(Expression ~ 1 + Line + Village + Location, data = df, family = nbinom2) ### The results are not dependent on the order or the covariates included so don't need to include extra covariates

# lrtest(model_nonrandom_nested, model_nonrandom)
# lrtest(model_nonrandom_no_rep_nested, model_nonrandom)


# model_nonrandom_noone <- glmmTMB(Expression ~ Line + Village + Replicate + Location, data = df, family = nbinom2) ### The results are not dependent on the order or the covariates included so don't need to include extra covariates
# model_nonrandom_reorder <- glmmTMB(Expression ~ 1 + Village + Replicate + Location + Line, data = df, family = nbinom2) ### The results are not dependent on the order or the covariates included so don't need to include extra covariates
# model_nonrandom_reorder_noone <- glmmTMB(Expression ~ Village + Replicate + Location + Line, data = df, family = nbinom2) ### The results are not dependent on the order or the covariates included so don't need to include extra covariates
# lm <- lm(Expression ~ Line + Village + Replicate + Location, data = df) ### The results are not dependent on the order or the covariates included so don't need to include extra covariates
# model_no1 <- glmmTMB(Expression ~ (1|Line) + (1|Village) + (1|Replicate) + (1|Location), data = df, family = nbinom2) ### The results are not dependent on the order or the covariates included so don't need to include extra covariates

# summary(model)
# summary(model_small)
# summary(model_nonrandom)
# summary(model_nonrandom_noone)
# summary(model_nonrandom_reorder)
# summary(model_nonrandom_reorder_noone)
# summary(model_no1)
# summary(lm)
# summary(model_reorder)




# model_line_x_location <- glmmTMB(Expression ~ 1 + (1|Line) + (1|Location), data = df, family = nbinom2) ### The results are not dependent on the order or the covariates included so don't need to include extra covariates
# model_line_x_location_rev <- glmmTMB(Expression ~ 1 + (1|Location) * (1|Line), data = df, family = nbinom2) ### The results are not dependent on the order or the covariates included so don't need to include extra covariates
# model_line_x_location_norandom <- glmmTMB(Expression ~ 1 + Line * Location, data = df, family = nbinom2) ### The results are not dependent on the order or the covariates included so don't need to include extra covariates
# model_line_x_location_norandom <- glmmTMB(Expression ~ 1 + Line * Location, data = df, family = nbinom2) ### The results are not dependent on the order or the covariates included so don't need to include extra covariates

# summary(model_line_x_location)
# summary(model_line_x_location_rev)
# summary(model_line_x_location_norandom)

# model_w_interaction <- glmmTMB(Expression ~ 1 + (1|Line) + (1|Village) + (1|Replicate) + (1|Location) + (1|Village:Line), data = df, family = nbinom2) ### The results are not dependent on the order or the covariates included so don't need to include extra covariates
# model_w_interaction2 <- glmmTMB(Expression ~ 1 + (1|Line) + (1|Village) + (1|Replicate) + (1|Location) + (1|Village:Location) + (1|Village:Line), data = df, family = nbinom2) ### The results are not dependent on the order or the covariates included so don't need to include extra covariates

# summary(model)
# summary(model_w_interaction)
# summary(model_w_interaction2)$coefficients$cond

# anova(model, model_w_interaction, test="LRT")
# lrtest(model, model_w_interaction)
# lrtest(model_w_interaction, model_w_interaction2)

# var_model <- variance_partition(model)
# rowSums(var_model[2:ncol(var_model)])

# var_model_w_interaction <- variance_partition(model_w_interaction)
# rowSums(var_model_w_interaction[2:ncol(var_model_w_interaction)])

# var_model_w_interaction2 <- variance_partition(model_w_interaction2)
# rowSums(var_model_w_interaction2[2:ncol(var_model_w_interaction2)])

# model_line_only_var_part <- variance_partition(glmmTMB(Expression ~ 1 + (1|Line), data = df, family = nbinom2))
# model_village_only_var_part <- variance_partition(glmmTMB(Expression ~ 1 + (1|Village), data = df, family = nbinom2))
# model_location_only_var_part <- variance_partition(glmmTMB(Expression ~ 1 + (1|Location), data = df, family = nbinom2))
# model_replicate_only_var_part <- variance_partition(glmmTMB(Expression ~ 1 + (1|Replicate), data = df, family = nbinom2))



nbGLM <- glm.nb(Expression ~ 1, data=df)
model <- glmmTMB(Expression ~ 1, data = df, family = nbinom2)
model_line <- glmmTMB(Expression ~ 1 + (1|Line) , data = df, family = nbinom2)
model_rep <- glmmTMB(Expression ~ 1 + (1|Replicate) , data = df, family = nbinom2)
model_site <- glmmTMB(Expression ~ 1 + (1|Location) , data = df, family = nbinom2)
model_village <- glmmTMB(Expression ~ 1 + (1|Village), data = df, family = nbinom2)
model_line_village <- glmmTMB(Expression ~ 1 + (1|Line) + (1|Village), data = df, family = nbinom2)
model_village_line <- glmmTMB(Expression ~ 1 + (1|Village) + (1|Line), data = df, family = nbinom2)
model_village_line_lm <- lm(Expression ~ 1 + Village + Line, data = df)

head(resid(model_village_line))
sum <- summary(model_village_line)
mean(resid(model_village_line)) + mean(resid(model_village_line))^2 * (1/sum$sigma)


summary(nbGLM)
summary(model)

mean(df$Expression) + mean(df$Expression)^2 * (1/nbGLM$theta)

model_sum <- summary(model)
summary(model_village_line)

1-logLik(model_line)/logLik(model)
1-logLik(model_village_line)/logLik(model_village)

1-logLik(model_village)/logLik(model)
1-logLik(model_line_village)/logLik(model_line)

1-logLik(model_rep)/logLik(model)


1-logLik(model_site)/logLik(model)



1-logLik(model_line_village)/logLik(model)
1-logLik(model_village_line)/logLik(model)



model2 <- glmmTMB(Expression ~ 1, data = df, family = nbinom2)
model_line2 <- glmmTMB(Expression ~ 1 + Line , data = df, family = nbinom2)
model_rep2 <- glmmTMB(Expression ~ 1 + Replicate , data = df, family = nbinom2)
model_site2 <- glmmTMB(Expression ~ 1 + Location , data = df, family = nbinom2)
model_village2 <- glmmTMB(Expression ~ 1 + Village, data = df, family = nbinom2)
model_interaction2 <- glmmTMB(Expression ~ 1 + Line + Village + Line:Village, data = df, family = nbinom2)
model_line_village2 <- glmmTMB(Expression ~ 1 +Line + Village, data = df, family = nbinom2)
model_village_line2 <- glmmTMB(Expression ~ 1 + Village + Line, data = df, family = nbinom2)
model_village_line_rep_site2 <- glmmTMB(Expression ~ 1 + Village + Line + Replicate + Location, data = df, family = nbinom2)
model_village_line_rep_site_interaction2 <- glmmTMB(Expression ~ 1 + Village + Line + Replicate + Location + Line:Village, data = df, family = nbinom2)

sum_vil_line <- summary(model_village_line2)


1-logLik(model_line2)/logLik(model2)
1-logLik(model_village_line2)/logLik(model_village2)

1-logLik(model_village2)/logLik(model2)
1-logLik(model_line_village2)/logLik(model_line2)

1-logLik(model_rep2)/logLik(model2)

1-logLik(model_interaction2)/logLik(model2)


1-logLik(model_site2)/logLik(model2)



1-logLik(model_line_village2)/logLik(model2)
1-logLik(model_village_line2)/logLik(model2)
1-logLik(model_village_line_rep_site2)/logLik(model2)
1-logLik(model_village_line_rep_site_interaction2)/logLik(model2)
1-logLik(model_village_line_rep_site_interaction2)/logLik(model_village_line_rep_site2)

anova(model_village_line_rep_site2, model_village_line_rep_site_interaction2)
pchisq(137.96, df=1, lower.tail=FALSE)/2




gene_mt <- "ENSG00000198899"

df_mt <- data.frame("Expression" = data.frame(seurat[["SCT"]]@counts[gene_mt,]), "Line" = seurat@meta.data$Final_Assignment, "Village" = ifelse(seurat@meta.data$Time == "Baseline", 0, 1), "Replicate" = gsub("[A-Z][a-z]+", "", seurat@meta.data$MULTI_ID), "Location" = seurat@meta.data$Location) ### The results are not dependent on the order or the covariates included so don't need to include extra covariates
colnames(df_mt)[1] <- "Expression"

model_mt <- glmmTMB(Expression ~ 1, data = df_mt, family = nbinom2)
model_line_mt <- glmmTMB(Expression ~ 1 + Line , data = df_mt, family = nbinom2)
model_rep_mt <- glmmTMB(Expression ~ 1 + Replicate , data = df_mt, family = nbinom2)
model_site_mt <- glmmTMB(Expression ~ 1 + Location , data = df_mt, family = nbinom2)
model_village_mt <- glmmTMB(Expression ~ 1 + Village, data = df_mt, family = nbinom2)
model_line_village_mt <- glmmTMB(Expression ~ 1 +Line + Village, data = df_mt, family = nbinom2)
model_village_line_mt <- glmmTMB(Expression ~ 1 + Village + Line, data = df_mt, family = nbinom2)
model_village_line_rep_site_mt <- glmmTMB(Expression ~ 1 + Village + Line + Replicate + Location, data = df_mt, family = nbinom2)
model_village_line_rep_site_interaction_mt <- glmmTMB(Expression ~ 1 + Village + Line + Replicate + Location + Line:Village, data = df_mt, family = nbinom2)



1-logLik(model_line_mt)/logLik(model_mt)
1-logLik(model_village_line_mt)/logLik(model_village_mt)

1-logLik(model_village_mt)/logLik(model_mt)
1-logLik(model_line_village_mt)/logLik(model_line_mt)

1-logLik(model_rep_mt)/logLik(model_mt)


1-logLik(model_site_mt)/logLik(model_mt)



1-logLik(model_line_village_mt)/logLik(model_mt)
1-logLik(model_village_line_mt)/logLik(model_mt)
1-logLik(model_village_line_rep_site_mt)/logLik(model_mt)
1-logLik(model_village_line_rep_site_interaction_mt)/logLik(model_mt)

anova(model_village_line_rep_site_mt, model_village_line_rep_site_interaction_mt)





df_lm <- data.frame("Expression" = data.frame(seurat[["SCT"]]@data[gene,]), "Line" = seurat@meta.data$Final_Assignment, "Village" = ifelse(seurat@meta.data$Time == "Baseline", 0, 1), "Replicate" = gsub("[A-Z][a-z]+", "", seurat@meta.data$MULTI_ID), "Location" = seurat@meta.data$Location) ### The results are not dependent on the order or the covariates included so don't need to include extra covariates
colnames(df_lm)[1] <- "Expression"



model <- lm(Expression ~ 1, data = df_lm)
model_line <- lm(Expression ~ 1 + (1|Line) , data = df)
model_rep <- lm(Expression ~ 1 + (1|Replicate) , data = df)
model_site <- lm(Expression ~ 1 + (1|Location) , data = df)
model_village <- lm(Expression ~ 1 + (1|Village), data = df)
model_line_village <- lm(Expression ~ 1 + (1|Line) + (1|Village), data = df)
model_village_line <- lm(Expression ~ 1 + (1|Village) + (1|Line), data = df)
model_village_line_lm <- lmer(Expression ~ 1 + (1|Village) + (1|Line) + (1|Replicate) + (1|Location), data = df_lm)
model_village_line_lm_inter <- lmer(Expression ~ 1 + (1|Village) + (1|Line) + (1|Replicate) + (1|Location) + (1|Replicate:Village), data = df_lm)
anova(model_village_line_lm, model_village_line_lm_inter)

model_village_line_lm2 <- lm(Expression ~ 1 + Village  + Replicate + Location +Line, data = df_lm)
model_village_line_lm_inter2 <- lm(Expression ~ 1 + Village + Line + Replicate + Location + Line:Village, data = df_lm)
model_village_line_lm_inter2b <- lm(Expression ~ 1 + Village + Line + Replicate + Location + Replicate:Village, data = df_lm)
model_village_line_lm_inter3 <- lm(Expression ~ 1 + Village + Line + Replicate + Location + Line:Village + Replicate:Village, data = df_lm)
model_village_line_lm_inter3b <- lm(Expression ~ 1 + Village + Line + Replicate + Location + Replicate:Village + Line:Village, data = df_lm)

anova(model_village_line_lm2, model_village_line_lm_inter2)
anova(model_village_line_lm2, model_village_line_lm_inter2b)
anova(model_village_line_lm_inter2, model_village_line_lm_inter3)
anova(model_village_line_lm_inter2b, model_village_line_lm_inter3b)


df_hier_unscale <- data.frame("Expression" = seurat[["SCT"]]@data[gene,], "Village" = as.factor(ifelse(seurat@meta.data$Time == "Baseline", 0, 1)), "Line" = seurat@meta.data$Final_Assignment, "Replicate" = as.factor(gsub("[A-Z][a-z]+", "", seurat@meta.data$MULTI_ID)), "Site" = seurat$Location)
colnames(df_hier_unscale)[1] <- "Expression"

# df_hier <- data.frame("Expression" = seurat[["SCT"]]@scale.data[gene,], "Village" = ifelse(seurat@meta.data$Time == "Baseline", 0, 1), "Line" = seurat@meta.data$Final_Assignment, "Replicate" = gsub("[A-Z][a-z]+", "", seurat@meta.data$MULTI_ID))
# colnames(df_hier)[1] <- "Expression"


model_rand_unscale <- lmer(Expression ~ (1 | Replicate) + (1 | Village) + (1 | Line) + (1 | Site) + (1|Replicate:Village) + (1|Replicate:Line) + (1|Replicate:Site) + (1|Village:Line) + (1|Village:Site) + (1|Line:Site),  data = df_hier_unscale)
# model_rand <- lmer(Expression ~ (1 | Replicate) + (1 | Village) + (1 | Line), data = df_hier)

hier_model_unscale <- hier.part(df_hier_unscale$Expression, df_hier_unscale[,2:ncol(df_hier_unscale)])
# hier_model <- hier.part(df_hier$Expression, df_hier[,2:ncol(df_hier)])


icc_specs(model_rand_unscale)
hier_model_unscale
# icc_specs(model_rand)






df_hier_unscale_mt <- data.frame("Expression" = seurat[["SCT"]]@data[gene_mt,], "Village" = as.factor(ifelse(seurat@meta.data$Time == "Baseline", 0, 1)), "Line" = seurat@meta.data$Final_Assignment, "Replicate" = as.factor(gsub("[A-Z][a-z]+", "", seurat@meta.data$MULTI_ID)), "Site" = seurat$Location)
colnames(df_hier_unscale_mt)[1] <- "Expression"

# df_hier_mt <- data.frame("Expression" = seurat[["SCT"]]@scale.data[gene_mt,], "Village" = as.factor(ifelse(seurat@meta.data$Time == "Baseline", 0, 1)), "Line" = seurat@meta.data$Final_Assignment, "Replicate" = as.factor(gsub("[A-Z][a-z]+", "", seurat@meta.data$MULTI_ID)))
# colnames(df_hier_mt)[1] <- "Expression"



model_rand_unscale_mt <- lmer(Expression ~ (1 | Replicate) + (1 | Village) + (1 | Line) + (1 | Site) + (1|Replicate:Village) + (1|Replicate:Line) + (1|Replicate:Site) + (1|Village:Line) + (1|Village:Site) + (1|Line:Site), data = df_hier_unscale_mt)
# model_rand_mt <- lmer(Expression ~ (1 | Replicate) + (1 | Village) + (1 | Line), data = df_hier_mt)
summary(model_rand_unscale_mt)
# summary(model_rand_mt)

hier_model_unscale_mt <- hier.part(df_hier_unscale_mt$Expression, df_hier_unscale_mt[,2:ncol(df_hier_unscale_mt)])
# hier_model_mt <- hier.part(df_hier_mt$Expression, df_hier_mt[,2:ncol(df_hier_mt)])



icc_specs(model_rand_unscale_mt)
hier_model_unscale_mt
# icc_specs(model_rand_mt)



df <- df_hier_unscale_mt

### To use for variance explained:
model <- lmer(Expression ~ (1 | Replicate) + (1 | Village) + (1 | Line) + (1 | Site), data = df)

base_model <- model

base_model <- lmer(Expression ~ (1 | Replicate) + (1 | Village) + (1 | Line) + (1 | Site) + (1|Replicate:Village) + (1|Replicate:Line) + (1|Replicate:Site) + (1|Village:Site) + (1|Line:Site), data = df)
model_2 <- lmer(Expression ~ (1 | Replicate) + (1 | Village) + (1 | Line) + (1 | Site) + (1|Replicate:Village) + (1|Replicate:Line) + (1|Replicate:Site) + (1|Village:Line) + (1|Village:Site) + (1|Line:Site), data = df)

anova_tmp <- anova(base_model, model_2)




dfb <- df_hier_unscale

### To use for variance explained:
base_modelb <- lmer(Expression ~ (1 | Replicate) + (1 | Village) + (1 | Line) + (1 | Site), data = dfb)

modelb <- lmer(Expression ~ (1 | Replicate) + (1 | Village) + (1 | Line) + (1 | Site) + (1|Replicate:Village) + (1|Replicate:Line) + (1|Replicate:Site) + (1|Village:Site) + (1|Line:Site), data = dfb)
model_all <- lmer(Expression ~ (1 | Replicate) + (1 | Village) + (1 | Line) + (1 | Site) + (1|Replicate:Village) + (1|Replicate:Line) + (1|Replicate:Site) + (1|Village:Line) + (1|Village:Site) + (1|Line:Site), data = dfb)


model_all_fixed <- lmer(Expression ~ 1 + (1 | Replicate) + (1 | Village) + (1 | Line) + (1 | Site) + (1|Replicate:Village) + (1|Replicate:Line) + (1|Replicate:Site) + (1|Village:Line) + (1|Village:Site) + (1|Line:Site), data = df_hier_unscale)

model_all_fast <- rint.regs(df_hier_unscale_update$Expression, data.matrix(df_hier_unscale_update[,c("Village", "Replicate", "Site")]), id = as.numeric(gsub("MBE1006", 3, gsub("TOB0421", 2, gsub("FSA0006", 1, df_hier_unscale_update$Line)))))


icc_specs(base_modelb)
icc_specs(modelb)
icc_specs(model_all)
icc_specs(model_all_fixed)

anova_tmpb <- anova(modelb, model_all)


ranova(model_all)
drop1(model_all)

model_all_step <- lmerTest::step(model_all)



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


##### Test time of different packages ####
tic("lme4 lmer time test, high line effect")
model_lme4 <- lme4::lmer(Expression ~ (1 | Replicate) + (1 | Village) + (1 | Line) + (1 | Site) + (1|Replicate:Village) + (1|Replicate:Line) + (1|Replicate:Site) + (1|Village:Line) + (1|Village:Site) + (1|Line:Site), data = df_hier_unscale)
toc()

tic("lmerTest lmer time test, high line effect")
model_lmertest <- lmerTest::lmer(Expression ~ (1 | Replicate) + (1 | Village) + (1 | Line) + (1 | Site) + (1|Replicate:Village) + (1|Replicate:Line) + (1|Replicate:Site) + (1|Village:Line) + (1|Village:Site) + (1|Line:Site), data = df_hier_unscale)
toc()

tic("glmmTMB time test, high line effect")
model_glmmtmb <- glmmTMB(Expression ~ (1 | Replicate) + (1 | Village) + (1 | Line) + (1 | Site) + (1|Replicate:Village) + (1|Replicate:Line) + (1|Replicate:Site) + (1|Village:Line) + (1|Village:Site) + (1|Line:Site), data = df_hier_unscale, REML = TRUE)
toc()

icc_specs(model_lme4)
icc_specs(model_lmertest)
icc_glmmtmb(model_glmmtmb)

model_glmmtmb_nested <- glmmTMB(Expression ~ (1 | Replicate) + (1 | Village) + (1 | Line) + (1 | Site) + (1|Replicate:Village) + (1|Replicate:Line) + (1|Replicate:Site) + (1|Village:Site) + (1|Line:Site), data = df_hier_unscale, REML = TRUE)

anova(model_glmmtmb, model_glmmtmb_nested)$`Pr(>Chisq)`



tic("lme4 lmer time test, high village effect")
model_lme4_mt <- lme4::lmer(Expression ~ (1 | Replicate) + (1 | Village) + (1 | Line) + (1 | Site) + (1|Replicate:Village) + (1|Replicate:Line) + (1|Replicate:Site) + (1|Village:Line) + (1|Village:Site) + (1|Line:Site), data = df_hier_unscale_mt)
toc()

tic("lmerTest lmer time test, high village effect")
model_lmertest_mt <- lmerTest::lmer(Expression ~ (1 | Replicate) + (1 | Village) + (1 | Line) + (1 | Site) + (1|Replicate:Village) + (1|Replicate:Line) + (1|Replicate:Site) + (1|Village:Line) + (1|Village:Site) + (1|Line:Site), data = df_hier_unscale_mt)
toc()

tic("glmmTMB time test, high line effect")
model_glmmtmb_mt <- glmmTMB(Expression ~ (1 | Replicate) + (1 | Village) + (1 | Line) + (1 | Site) + (1|Replicate:Village) + (1|Replicate:Line) + (1|Replicate:Site) + (1|Village:Line) + (1|Village:Site) + (1|Line:Site), data = df_hier_unscale_mt, REML = TRUE)
toc()

icc_specs(model_lme4_mt)
icc_specs(model_lmertest)
icc_glmmtmb(model_glmmtmb_mt)
sum(icc_glmmtmb(model_glmmtmb_mt)$vcov)



##### Test comparisons with different order/number of covariates #####
model_glmmtmb_short <- glmmTMB(Expression ~ (1 | Replicate) + (1 | Village) + (1 | Line) + (1 | Site), data = df_hier_unscale, REML = TRUE)
model_glmmtmb_short_interaction <- glmmTMB(Expression ~ (1 | Replicate) + (1 | Village) + (1 | Line) + (1 | Site) + (1 | Village:Line), data = df_hier_unscale, REML = TRUE)

anova(model_glmmtmb_short, model_glmmtmb_short_interaction)$`Pr(>Chisq)`

ranova_mt <- ranova(model_glmmtmb_mt)
