library(data.table)
library(ggplot2)


##### Set up Directories #####
data <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Growth_rate/18-line_village/growth_measurements.tsv"
outdir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/growth_rate/18_line_village/"

dir.create(outdir, recursive = TRUE)


##### Read in data #####
data <- fread(data, sep = "\t")
colors_original <- readRDS("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Nona_multiome/line_colors")



##### Pivot data longer (reps) #####
data_long <- melt(data, id.vars = c("Plate", "Cell_Line", "Day", "Sector"),
                measure.vars = c("rep1", "rep2", "rep3"),variable.name = "Replicate", value.name = "Confluence")




##### Plot #####
plot <- ggplot(data_long, aes(Day, Confluence, color = Cell_Line, shape = Replicate)) +
    geom_point(alpha = 0.6) +
    geom_smooth(se=FALSE, size = 0.5, alpha = 0.6) +
    scale_color_manual(values = colors_original) 
    # scale_linetype_manual(values=c("longdash", "dashed", "dotted"))

ggsave(plot, filename = paste0(outdir, "growth_plot.png"))
# , aes(linetype=Replicate)


unique(data_long$Cell_Line)[!unique(data_long$Cell_Line) %in% names(colors_original)]



##### Normalize confluency #####
data$rep1_norm <- as.numeric(NA)
data$rep2_norm <- as.numeric(NA)
data$rep3_norm <- as.numeric(NA)

for (line in unique(data$Cell_Line)){
    for (quad in unique(data$Sector)){
        for (day in unique(data[Cell_Line == line & Sector == quad]$Day)){
            data[Cell_Line == line & Sector == quad & Day == day]$rep1_norm <- data[Cell_Line == line & Sector == quad & Day == day]$rep1/data[Cell_Line == line & Sector == quad & Day == 1]$rep1
            data[Cell_Line == line & Sector == quad & Day == day]$rep2_norm <- data[Cell_Line == line & Sector == quad & Day == day]$rep2/data[Cell_Line == line & Sector == quad & Day == 1]$rep2
            data[Cell_Line == line & Sector == quad & Day == day]$rep3_norm <- data[Cell_Line == line & Sector == quad & Day == day]$rep3/data[Cell_Line == line & Sector == quad & Day == 1]$rep3
        }
    }
}


data_norm_long <-  melt(data, id.vars = c("Plate", "Cell_Line", "Day", "Sector"),
                measure.vars = c("rep1_norm", "rep2_norm", "rep3_norm"),variable.name = "Replicate", value.name = "Confluence")



##### Plot #####
plot_norm <- ggplot(data_norm_long[Day < 6], aes(Day, Confluence, color = Cell_Line, shape = Replicate)) +
    geom_point(alpha = 0.6) +
    geom_smooth(se=FALSE, size = 0.5, alpha = 0.6) +
    scale_color_manual(values = colors_original) 
    # scale_linetype_manual(values=c("longdash", "dashed", "dotted"))

ggsave(plot_norm, filename = paste0(outdir, "growth_plot_norm.png"))





##### make files for running ratrack #####
## Different file for each cell line and replicate
ratrack_outdir <- paste0(outdir, "/ratrack/data/")
dir.create(ratrack_outdir)

ratrack_dt <- list()

for (line in unique(data_long$Cell_Line)){
    # for (rep in unique(data_long$Replicate)){
        tmp <- data_long[Cell_Line == line]
        ratrack_dt[[line]] <- data.table(name = paste0(tmp$Cell_Line, "_", tmp$Replicate), time = tmp$Day, count = tmp$Confluence*10000, sample1 = 1/9)

        fwrite(ratrack_dt[[line]], paste0(ratrack_outdir, line, "_ratrack.csv"), sep = ",")
    # }
}




ratrack_outdir_n <- paste0(outdir, "/ratrack/data/limited_times/")
dir.create(ratrack_outdir_n)

cell_n <- data.table(Cell_Line = c("166", "180N", "FSA0001", "FSA0004", "IST1877", "IST3323", "MBE0953", "MBE2817", "MBE2900", "TOB0198", "TOB0199", "TOB0205","TOB0220","TOB0421","TOB0435", "WAB0004", "WAB0038", "WAB0103"),
            N = c(4, 4, 3, 3, 4, 4, 3, 4, 3, 4,3,3,3,3,3,4,8,3))

for (line in unique(data_long$Cell_Line)){
    tmp <- data_long[Cell_Line == line]
    ratrack_dt[[line]] <- data.table(name = paste0(tmp$Cell_Line, "_", tmp$Replicate), time = tmp$Day, count = tmp$Confluence*10000, sample1 = 1/9)
    ratrack_dt[[line]] <- ratrack_dt[[line]][time <= cell_n[Cell_Line == line]$N]

    fwrite(ratrack_dt[[line]], paste0(ratrack_outdir_n, line, "_ratrack.csv"), sep = ",")
}





# 166 - days 0-4
# 180N - days 0-4 (or maybe 3 if that seems off)
# FSA0001 - days 0-3
# FSA0004 - days 0-3
# IST1877 - days 0-4 (or maybe 3 if that seems off)
# IST3323 - days 0-4 
# MBE0953 - days 0-3
# MBE2817 - days 0-4, remove 
# MBE2900 - days 0-3
# TOB0198 - days 0-4 (or maybe 3 if that seems off)
# TOB0199 - days 0-3
# TOB0205 - days 0-3 (or maybe 2 if that seems off)
# TOB0220 - days 0-3 (or maybe 2 if that seems off)
# TOB0421 - days 0-3 (or maybe 2 if that seems off)
# TOB0435 - days 0-3 
# WAB0004 - days 0-4
# WAB0038 - days 0-8
# WAB0103 - days 0-3 (or maybe 2 if that seems off)




ratrack_outdir_4days <- paste0(outdir, "/ratrack/data/days4/")
dir.create(ratrack_outdir_4days)

cell_n <- data.table(Cell_Line = c("166", "180N", "FSA0001", "FSA0004", "IST1877", "IST3323", "MBE0953", "MBE2817", "MBE2900", "TOB0198", "TOB0199", "TOB0205","TOB0220","TOB0421","TOB0435", "WAB0004", "WAB0038", "WAB0103"),
            N = c(4, 4, 3, 3, 4, 4, 3, 4, 3, 4,3,3,3,3,3,4,8,3))

for (line in unique(data_long$Cell_Line)){
    tmp <- data_long[Cell_Line == line]
    ratrack_dt[[line]] <- data.table(name = paste0(tmp$Cell_Line, "_", tmp$Replicate), time = tmp$Day, count = tmp$Confluence*10000, sample1 = 1/9)
    ratrack_dt[[line]] <- ratrack_dt[[line]][time <= 4]

    fwrite(ratrack_dt[[line]], paste0(ratrack_outdir_4days, line, "_ratrack.csv"), sep = ",")
}

