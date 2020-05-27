suppressPackageStartupMessages({
    library(TreeHeatmap)
    library(ggnewscale)
    library(ggplot2)
    library(ggnewscale)
    library(dplyr)
    library(ggtree)
    library(cowplot)
    library(scales)
})

# s <- 5
# load functions
source("summary/drawF.R")   

if (reso == "high") {
source("summary/rm_ancestor.R")
}

# tree heatmaps: BS, US, SS
pathD <- "output/RData/DataPrep"
pathA <- "output/RData/Analysis"
ss <- c("BS_sp25_sim5_r4_L2",
        "US_sp25_sim5_r4_L2",
        "SS_sp25_sim5_r4_pct0.4_L2")
mm <- c("treeclimbR", "StructFDR", "HFDR", "minP",
        "BH", "miLineage", "LassoGLM")
scene <- c("BS", "US", "SS")

# load results of lefse
load("lefse/output/nodes/out_lefse.RData")
# --------------------------------- result -------------------------------------
loc_list <-  tse_list <- vector("list", length(ss))

for (i in 1:3) {
    message(i)
    si <- ss[i]
    
    # load data
    dataFile <- file.path(pathD, paste0(si, ".RData"))
    resultFile <- file.path(pathA, si, paste0(mm, ".RData"))
    lapply(c(dataFile, resultFile), load, .GlobalEnv)
    
    #
    tse_list[[i]] <- tse
    loc_da <- list(
        lasso = loc.Lasso[[s]],
        HFDR = loc.HFDR_0.05[[s]],
        StructFDR = loc.StructFDR_0.05[[s]],
        BH = loc.bh_0.05[[s]],
        miLineage1 = loc1_0.05.MLA[[s]],
        miLineage2 = loc2_0.05.MLA[[s]],
        minP = outMin[[s]]$node[outMin[[s]]$keep_0.05],
        lefse = out_lefse_low[[si]]$`0.05`[[s]],
        treeclimbR = outsel_0.05[[s]]$node[outsel_0.05[[s]]$signal.node])
    loc_list[[i]] <- lapply(loc_da, FUN = function(x){
        if(length(x)) {x} else{ NULL }
        if (reso == "high") {
            rm_ancestor(node = x, tree = rowTree(tse))
        } else{x}
    })
}




fig_list <- lapply(seq_along(loc_list), FUN = function(x, r) {
    loc_da <- loc_list[[x]]
    drawF(tse = tse_list[[x]], s = r, 
          anno_node = loc_da$treeclimbR, 
          loc_da = loc_da, scene = scene[x])
}, r = s)

# legends
legend <- get_legend(fig_list[[1]])
fig_bottom <- plot_grid(fig_list[[1]] +
                        theme(legend.position = "none",
                              plot.margin = unit(c(0, 2, 0, 0),"mm")),
                    fig_list[[2]] +
                        theme(legend.position = "none",
                              plot.margin = unit(c(0, 2, 0, 0),"mm")),
                    fig_list[[3]] +
                        theme(legend.position = "none",
                              plot.margin = unit(c(0, 2, 0, 0),"mm")),
                    legend, rel_widths = c(0.3, 0.3, 0.3, 0.1),
                    nrow = 1)


fig_bottom
