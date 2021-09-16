#library(lme4)
library(haven)
library(ggplot2)
library(lme4)
library(glmmTMB)
library(Seurat)



##### Bring in variables #####
### Bring in arguments
args <- commandArgs(trailingOnly = TRUE)
outdir <- paste0(args[1])
number <- as.numeric(args[2])




dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"




##### Read in seurat with genes #####
seurat <- readRDS(paste0(dir,"output/Distribution_tests/seurat_integrated_all_times_clustered_1pct_expressing.rds"))

seurat_sub <- subset(seurat, subset = Time != "Thawed Village Day 0")
seurat_sub <- subset(seurat_sub, subset = Time != "Thawed Village Day 7")
gene <- rownames(seurat[["SCT"]]@counts)[number]

seurat_sub@meta.data$Location <- gsub("_.+", "", seurat_sub@meta.data$Location)

seurat_sub_list  <- list()

for (location in unique(seurat_sub@meta.data$Location)){
	seurat_sub_list[[location]] <- subset(seurat_sub, subset = Location == location)
}


df_list <- lapply(seurat_sub_list, function(x){
	# tmp <- data.frame("Expression" = data.frame(x[["SCT"]]@counts[gene,]), "Line" = x@meta.data$Final_Assignment, "Village" = ifelse(x@meta.data$Time == "Baseline", 0, 1), "Replicate" = gsub("[A-Z][a-z]+", "", x@meta.data$MULTI_ID),"Phase" = x@meta.data$phases,  "MT" = x@meta.data$percent.mt, "RB" = x@meta.data$percent.rb)
	tmp <- data.frame("Expression" = data.frame(x[["SCT"]]@counts[gene,]), "Line" = x@meta.data$Final_Assignment, "Village" = ifelse(x@meta.data$Time == "Baseline", 0, 1), "Replicate" = gsub("[A-Z][a-z]+", "", x@meta.data$MULTI_ID)) ### The results are not dependent on the order or the covariates included so don't need to include extra covariates
	colnames(tmp)[1] <- "Expression"
	return(tmp)
})


# Fit model
models <- lapply(df_list, function(x){
	# glmmTMB(Expression ~ 1 + (1|Line) + (1|Village) + (1|Replicate) + (1|Phase) + (1|MT) + (1|RB), data = x, family = nbinom2)
	glmmTMB(Expression ~ 1 + (1|Line) + (1|Village) + (1|Replicate), data = x, family = nbinom2) ### The results are not dependent on the order or the covariates included so don't need to include extra covariates
})

summaries <- lapply(models, function(x){
	summary(x)
})



beta0 <- list()
sigma2u <- list()
alpha <- list()
expectation <- list()
variance <- list()
variance2 <- list()
variance1 <- list()
vpc_variable <- list()
vpc_resid <- list()

for (location in names(models)){
	# Intercept
	beta0[[location]] <- summary(models[[location]])$coefficients$cond[1,1]


	vpc_variable[[location]] <- data.frame(matrix(nrow=1,ncol=(length(names(summary(models[[location]])$varcor$cond))) + 1))
	colnames(vpc_variable[[location]]) <- c("Gene", names(summary(models[[location]])$varcor$cond))
	vpc_variable[[location]]$Gene <- gene
	vpc_resid[[location]] <- data.frame(matrix(nrow=1,ncol=(length(names(summary(models[[location]])$varcor$cond))) + 1))
	colnames(vpc_resid[[location]]) <- c("Gene", names(summary(models[[location]])$varcor$cond))
	vpc_resid[[location]]$Gene <- gene


	for (variable in names(summary(models[[location]])$varcor$cond)){
		# Cluster variance
		sigma2u[[location]][[variable]] <- summary(models[[location]])$varcor$cond[[variable]][1,1]


		# Overdispersion parameter
		alpha[[location]] <- 1/(summary(models[[location]])$sigma)

		# Marginal expectation
		expectation[[location]][[variable]] <- exp(beta0[[location]] + sigma2u[[location]][[variable]]/2)


		# Marginal variance
		variance[[location]][[variable]] <- expectation[[location]][[variable]] + expectation[[location]][[variable]]^2*(exp(sigma2u[[location]][[variable]])*(1 + alpha[[location]]) - 1)

		# Marginal variance: Level-2 component
		variance2[[location]][[variable]] <- expectation[[location]][[variable]]^2*(exp(sigma2u[[location]][[variable]]) - 1)

		# Marginal variance: Level-1 component
		variance1[[location]][[variable]] <- expectation[[location]][[variable]] + expectation[[location]][[variable]]^2*exp(sigma2u[[location]][[variable]])*alpha[[location]]

		# Level-2 VPC
		vpc_variable[[location]][1,variable] <- variance2[[location]][[variable]]/(variance2[[location]][[variable]] + variance1[[location]][[variable]])

		# Level-1 VPC
		vpc_resid[[location]][1,variable] <- variance1[[location]][[variable]]/(variance2[[location]][[variable]] + variance1[[location]][[variable]])

	}

}

saveRDS(vpc_variable, paste0(outdir, gene, "_variable_variance_partition_coefficients.rds"))

