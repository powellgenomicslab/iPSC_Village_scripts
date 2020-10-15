library(tidyr)
library(tidyverse)
library(dplyr)
library(ggplot2)

dir <- "/Volumes/ScratchRoot/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Seurat/"
harmony_dir <- paste0(dir,"Multi_Resolution_Harmony_DayPatient/")
seurat_dir <- paste0(dir,"Multi_Resolution_Integration_DayPatient/")

harmony <- read_delim(paste0(harmony_dir,"Identities_Resolution_0.14.txt"), delim = "\t")
seurat <- read_delim(paste0(seurat_dir,"Identities_Resolution_0.02.txt"), delim = "\t")

df <- cbind(harmony,seurat)
colnames(df) <- c("harmony","seurat")
df$test <- rownames(df)

ggplot(df, aes(harmony, fill = as.factor(seurat))) +
  geom_bar(position = "fill") +
  theme_classic()
