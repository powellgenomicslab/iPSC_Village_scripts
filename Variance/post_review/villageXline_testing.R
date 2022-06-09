library(data.table)
library(Seurat)
library(pkgcond)
library(glmmTMB)



### Set up directories and variables ###
args <- commandArgs(trailingOnly = TRUE)
icc_outdir <- paste0(args[1])
model_outdir <- paste0(args[2])
resid_outdir <- paste0(args[3])
ensg <- as.character(args[4])

print(icc_outdir)
print(model_outdir)
print(resid_outdir)
print(ensg)



### Read in data ###
seurat <- readRDS("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Distribution_tests/seurat_integrated_all_times_clustered_1pct_expressing.rds")
icc_dt <- fread("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/variance_partition_post_review/combined/sig_results.tsv.gz", sep = "\t")



### Process counts ###
## Subset non-cryo
seurat@meta.data$Location <- gsub("_Baseline", "", seurat@meta.data$Location) %>% gsub("_Village.+", "", .) %>% gsub("Thawed", "Cryopreserved",.)
seurat@meta.data$Time <- gsub("Thawed Village Day 0", "Baseline", seurat@meta.data$Time) %>% gsub("Thawed Village Day 7", "Village", .) %>% gsub(" Day 4", "", .)

seurat <- subset(seurat, subset = Location != "Sydney_Cryopreserved")





## Subset village and non-village
seurat_village <- subset(seurat, subset = Time == "Village")
seurat_nonvillage <- subset(seurat, subset = Time == "Baseline")

vars <- icc_dt[gene == ensg]$grp
vars <- vars[!(vars %in% c("Village:Line", "Residual"))]
vars <- gsub("Village:","", vars) %>% gsub(":Village", "", .)


### Make DF for modeling ###
df_village <- data.frame("Expression" = seurat_village[["SCT"]]@data[gene,], "Line" = seurat_village@meta.data$Final_Assignment, "Replicate" = as.factor(gsub("[A-Z][a-z]+", "", seurat_village@meta.data$MULTI_ID)), "Site" = seurat_village$Location)
colnames(df_village)[1] <- "Expression"

df_nonvillage <- data.frame("Expression" = seurat_nonvillage[["SCT"]]@data[gene,], "Line" = seurat_nonvillage@meta.data$Final_Assignment, "Replicate" = as.factor(gsub("[A-Z][a-z]+", "", seurat_nonvillage@meta.data$MULTI_ID)), "Site" = seurat_nonvillage$Location)
colnames(df_nonvillage)[1] <- "Expression"



### 1. Fit significant variables except line effect ###
model1 <- as.formula(paste0("Expression ~ ", paste0(vars, collapse = " + ")))

# model_all <- as.formula(paste0("Expression ~ (1|", paste0(vars, collapse = ") + (1|"), ")"))
# model_test <- suppress_warnings(glmmTMB(formula = noquote(model_all), data = df_village, REML = TRUE), "giveCsparse")
# model_test2 <- suppress_warnings(glmmTMB(formula = noquote(model_all), data = df_nonvillage, REML = TRUE), "giveCsparse")

# icc_glmmtmb(model_test)
# icc_glmmtmb(model_test2)


model_village <- suppress_warnings(glmmTMB(formula = noquote(model1), data = df_village, REML = TRUE), "giveCsparse")
model_nonvillage <- suppress_warnings(glmmTMB(formula = noquote(model1), data = df_nonvillage, REML = TRUE), "giveCsparse")



### 2. Fit residuals with Line ###
df_village$resid <- resid(model_village)
df_nonvillage$resid <- resid(model_nonvillage)

model_resid_village <- suppress_warnings(glmmTMB(resid ~ Line, data = df_village, REML = TRUE), "giveCsparse")
model_resid_nonvillage <- suppress_warnings(glmmTMB(resid ~ Line, data = df_nonvillage, REML = TRUE), "giveCsparse")



### 3. Capture beta of village and non-village in matrix ###
results <- data.table(Gene = gene,
                        beta_Village = ,
                        beta_Nonvillage = ,
                        LRT_Village = ,
                        LRT_Nonvillage =)


### LRT between no line and line effect and add to table ###



### Write out table ###
