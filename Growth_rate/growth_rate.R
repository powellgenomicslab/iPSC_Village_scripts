library(data.table)
library(tidyverse)
library(gtools)
library(growthcurver)
library("ggpubr")


dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/"
data <- fread(paste0(dir, "Growth_rate/growth_rate_measurements.csv"))
outdir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/growth_rate/"


dir.create(outdir, recursive = TRUE)


models <- lapply(unique(data$CellLine), function(line){
	SummarizeGrowth(data[CellLine == line]$Time, data[CellLine == line]$Count)
})
names(models) <- unique(data$CellLine)

predict <- lapply(names(models), function(line){
	tmp <- data.table(Time = data[CellLine == line]$Time, pred.wt = predict(models[[line]]$model))
	tmp$CellLine <- line
	return(tmp)
})

predict_df <- do.call(rbind, predict)




### Make Plot ###
pBasic <- ggplot(data, aes(Time, Count, color = CellLine)) +
  geom_point() +
  theme_classic() +
  geom_line(data=predict_df, aes(y=pred.wt))

ggsave(pBasic, filename = paste0(outdir, "/Basic_growth_rates_plot.png"))




##### Correlate growth rate with confluency #####
pCorr <- ggplot(data, aes(Count,Confluency, color = CellLine)) +
			geom_point() +
			theme_classic()

ggsave(pCorr, filename = paste0(outdir, "Count_confluency_corr.png"))


pCorr <- ggscatter(data, x = "Count", y = "Confluency", color = "CellLine",
			add = "reg.line", conf.int = TRUE, 
			cor.method = "pearson") +
			stat_cor(aes(color = CellLine), label.x = 90)
ggsave(pCorr, filename = paste0(outdir, "Count_confluency_corr.png"))
