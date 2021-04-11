# === tool ========
.libPaths(c(Sys.getenv('R_LIBS_1'), 
            Sys.getenv('R_LIBS_2'),
            Sys.getenv('R_LIBS_3')))
.libPaths()

suppressPackageStartupMessages({
    library(dplyr)
    library(stringr)
    library(ggplot2)
})

R.Version()
# ==== arguments from batch R=====================
argsList <- (commandArgs(trailingOnly = TRUE))
argsList
for (i in seq_along(argsList)) {
    if(length(argsList[[i]])>1){
        assign(names(argsList)[i],argsList[[i]])
    }else{eval(parse(text = argsList[[i]]))}
}


## === load data from data preparation step ======
print(inPath)
print(outPath)

#main <- "/Volumes/fiona/phd/treeclimbR_article/timeBenchmark/output/RData/Analysis"

dataFiles <- file.path(inPath, list.files(path = inPath))
sapply(dataFiles, function(file) {load(file, envir=globalenv())}) 

# function
getDF <- function(method_time, method_name) {
    method_time <- lapply(method_time, FUN = function(x){x["elapsed"]/3600}) 
    scene <- strsplit(names(method_time), "_")
    data.frame(sample = sapply(scene, FUN = function(x){x[[1]]}),
               tree_size = sapply(scene, FUN = function(x){x[[2]]}),
               time_min = unlist(method_time),
               Method = method_name) %>%
        mutate(sample = as.numeric(gsub(pattern = "sample", "", sample))*2,
               tree_size = as.numeric(gsub(pattern = "leaf", "", tree_size)))
}


lefse <- lefse[-grep("leaf10000", names(lefse))]
res <- list("miLineage1" = miLineage1,
            "miLineage2" = miLineage2,
            "treeclimbR" = treeclimbR_all,
            "lasso" = Lasso_all,
            "BH" = BH,
            "HFDR" = HFDR,
            "minP" = minP_all,
            "StructFDR" = structFDR,
            "LEfSe" = lefse)
df <- lapply(seq_along(res), FUN = function(x) {
  getDF(method_time = res[[x]], method_name = names(res)[[x]])
})

DF <- do.call(rbind, df) 


head(DF)

# Generate the figure
vcolor <- c("treeclimbR" = "#E41A1C", "BH" = "#377EB8",
            "StructFDR" = "#4DAF4A", "HFDR" = "#984EA3",
            "lasso" = "#FF7F00", "minP" = "#A65628",
            "miLineage1" = "#999999", "miLineage2" = "#666666",
            "LEfSe" = "#E7298A")

vline <- c("treeclimbR" = "solid", "BH" = "dashed",
            "StructFDR" = "dashed", "HFDR" = "dashed",
            "lasso" = "dashed", "minP" = "dashed",
            "miLineage1" = "dashed", "miLineage2" = "dashed",
            "LEfSe" = "dashed")

prettify <- theme_bw(base_size = 10) + theme(
    aspect.ratio = 1,
    #plot.margin = unit(c(0, 0, 0.2, 0), "cm"),
    panel.grid = element_blank(),
    panel.spacing = unit(0, "lines"),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.key.size= unit(1.5, "mm"),
    legend.spacing.x = unit(0.5, "mm"),
    plot.title = element_text(hjust = 0.5),
    legend.text = element_text(size = 7),
    legend.position="right",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin=margin(-5,-5,-10,-5),
    strip.background = element_rect(colour = "black", fill = "gray90"),
    strip.text.x = element_text(color = "black", size = 8),
    strip.text.y = element_text(color = "black", size = 8))


fig <- ggplot(DF, aes(x = tree_size, y = time_min, color = Method)) +
    geom_point() +
    geom_line(aes(linetype = Method)) +
    facet_grid(cols = vars(sample)) +
    scale_linetype_manual("Method", values = vline) +
    scale_color_manual("Method", values = vcolor) +
    scale_x_sqrt(breaks = c(100, 500, 1000, 3000, 6000, 10000)) +
    scale_y_sqrt(breaks = c(0.5, 1, 2, 5, 10)) +
    labs(x = "Tree sizes", y = "Time (h)") +
    theme(axis.text = element_text(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(outPath, fig, units = "in", width = 7, height = 3.5,
       dpi = 300)

sessionInfo()