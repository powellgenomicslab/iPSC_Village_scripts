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
location <- as.character(args[3])
village <- as.character(args[4])
datadir <- as.character(args[5])


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


SCT <- readRDS(paste0(datadir,  village, "_", location,"_seurat.rds"))

# for (rep in uniuqe(SCT_sub@meta.data$MULTI_ID)){
# 	SCT_sub <- subset(SCT, subset = MULTI_ID == rep)
# }

if (number <= nrow(SCT[["SCT"]]@counts)){
	gene <- rownames(SCT[["SCT"]]@counts)[number]

	### Make dataframes for modeling ###
	df <- data.frame("Expression" = data.frame(SCT[["SCT"]]@counts[gene,]), "Line" = SCT@meta.data$Final_Assignment, "Replicate" = gsub("[A-Z][a-z]+", "", SCT@meta.data$MULTI_ID)) ### The results are not dependent on the order or the covariates included so don't need to include extra covariates
	colnames(df)[1] <- "Expression"

	# df_list <- lapply(SCT_sub, function(x){
	# 	tmp <- data.frame("Expression" = data.frame(x[["SCT"]]@counts[gene,]), "Line" = x@meta.data$Final_Assignment) ### The results are not dependent on the order or the covariates included so don't need to include extra covariates
	# 	colnames(tmp)[1] <- "Expression"
	# 	return(tmp)
	# })


	### Fit models
	model <- glmmTMB(Expression ~ 1 + (1|Line) + (1|Replicate), data = df, family = nbinom2) ### The results are not dependent on the order or the covariates included so don't need to include extra covariates


	# model_list <- lapply(df_list, function(x){
	# 	glmmTMB(Expression ~ 1 + (1|Line), data = y, family = nbinom2) ### The results are not dependent on the order or the covariates included so don't need to include extra covariates
	# })


	##### Get Variance Partitions #####
	variance_partition <- variance_partition(model)

	# variance_partition_list <- lapply(model_list, variance_partition)

	saveRDS(variance_partition, paste0(outdir, village, "_", location, "_", gene, "_variable_variance_partition_coefficients.rds"))

	# saveRDS(variance_partition_list, paste0(outdir, gene, "_variable_variance_partition_coefficients.rds"))
}





