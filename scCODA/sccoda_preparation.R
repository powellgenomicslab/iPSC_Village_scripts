library(data.table)
library(tidyverse)
library(Seurat)




outdir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/scCODA/"
dir.create(outdir, recursive = TRUE)



seurat <- readRDS("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/All_data_integrated_remove_bad/seurat_integrated_all_times_clustered.rds")


head(seurat@meta.data)

seurat@meta.data$Location <- gsub("_.+", "", seurat@meta.data$Location)
seurat@meta.data$Cryopreserved <- ifelse(grepl("Thawed", seurat@meta.data$Location_Time), "Cryopreserved", "Fresh")
seurat@meta.data$Location_Cryopreserved_Line <- paste0(seurat@meta.data$Location, "_", seurat@meta.data$Cryopreserved, "_", seurat@meta.data$Final_Assignment)
seurat@meta.data$Village <- gsub("Baseline", "Uni-Culture", seurat@meta.data$Time) %>% 
        gsub("Thawed Village Day 0", "Uni-Culture", .) %>%
        gsub("Thawed Village Day 7", "Village", .) %>%
        gsub("Village Day 4", "Village", .)
seurat@meta.data$Replicate <- gsub("Brisbane", "Replicate", seurat@meta.data$MULTI_ID) %>%
        gsub("Melbourne", "Replicate", .) %>% 
        gsub("Sydney", "Replicate", .)


seurat@meta.data$Group <- paste0(seurat@meta.data$Location, "_", seurat@meta.data$Village, "_",seurat@meta.data$Replicate, "_", seurat@meta.data$Cryopreserved)


cell_prop_long <- data.table(table(seurat@meta.data$Group, seurat@meta.data$Final_Assignment))


cell_prop_dt <- dcast(cell_prop_long, V1 ~ V2, value.var = "N")
colnames(cell_prop_dt)[1] <- c("Group")


fwrite(cell_prop_dt, paste0(outdir, "cell_line_numbers.tsv"), sep = "\t")
