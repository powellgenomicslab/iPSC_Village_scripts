library(haven)
library(ggplot2)
library(lme4)
library(glmmTMB)
library(Seurat)



##### Bring in variables #####
### Bring in arguments
args <- commandArgs(trailingOnly = TRUE)
outdir <- as.character(args[1])
number <- as.numeric(args[2])
rep <- as.character(args[3])
datadir <- as.character(args[4])


##### Make variance partitioning function #####
variance_partition <- function(model_out){
	sigma2u <- list()
	expectation <- list()
	variance <- list()
	variance2 <- list()
	variance1 <- list()

	beta0 <- summary(model_out)$coefficients$cond[1,1]


	vpc_variable <- data.frame(matrix(nrow=1,ncol=(length(names(summary(model_out)$varcor$cond))) + 1))
	colnames(vpc_variable) <- c("Gene", names(summary(model_out)$varcor$cond))
	vpc_variable$Gene <- gene
	vpc_resid <- data.frame(matrix(nrow=1,ncol=(length(names(summary(model_out)$varcor$cond))) + 1))
	colnames(vpc_resid) <- c("Gene", names(summary(model_out)$varcor$cond))
	vpc_resid$Gene <- gene


	for (variable in names(summary(model_out)$varcor$cond)){
		# Cluster variance
		sigma2u[[variable]] <- summary(model_out)$varcor$cond[[variable]][1,1]


		# Overdispersion parameter
		alpha <- 1/(summary(model_out)$sigma)

		# Marginal expectation
		expectation[[variable]] <- exp(beta0 + sigma2u[[variable]]/2)


		# Marginal variance
		variance[[variable]] <- expectation[[variable]] + expectation[[variable]]^2*(exp(sigma2u[[variable]])*(1 + alpha) - 1)

		# Marginal variance: Level-2 component
		variance2[[variable]] <- expectation[[variable]]^2*(exp(sigma2u[[variable]]) - 1)

		# Marginal variance: Level-1 component
		variance1[[variable]] <- expectation[[variable]] + expectation[[variable]]^2*exp(sigma2u[[variable]])*alpha

		# Level-2 VPC
		vpc_variable[1,variable] <- variance2[[variable]]/(variance2[[variable]] + variance1[[variable]])

		# Level-1 VPC
		vpc_resid[1,variable] <- variance1[[variable]]/(variance2[[variable]] + variance1[[variable]])
	}
	return(vpc_variable)
}



dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"


SCT <- readRDS(paste0(datadir, rep, "_seurat_1pct.rds"))


if (number <= nrow(SCT[["SCT"]]@counts)){
	gene <- rownames(SCT[["SCT"]]@counts)[number]

	### Make dataframes for modeling ###
	df <- data.frame("Expression" = data.frame(SCT[["SCT"]]@counts[gene,]), "Line" = SCT@meta.data$Final_Assignment, "Village" = ifelse(SCT@meta.data$Time == "Baseline", 0, 1)) ### The results are not dependent on the order or the covariates included so don't need to include extra covariates
	colnames(df)[1] <- "Expression"


	### Fit models
	model <- glmmTMB(Expression ~ 1 + (1|Line) + (1|Village), data = df, family = nbinom2) ### The results are not dependent on the order or the covariates included so don't need to include extra covariates


	##### Get Variance Partitions #####
	variance_partition <- variance_partition(model)


	saveRDS(variance_partition, paste0(outdir, rep, "_", gene, "_variable_variance_partition_coefficients.rds"))

}





q()
