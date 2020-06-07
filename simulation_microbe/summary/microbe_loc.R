suppressPackageStartupMessages({
    library(TreeSummarizedExperiment)
    library(treeclimbR)
    library(ggtree)
    library(cowplot)
    library(ggplot2)
    
})

dir <- c("LassoGLM.RData", "HFDR.RData", "StructFDR.RData",
         "BH.RData", "miLineage.RData", "minP.RData",
         "treeclimbR.RData")
mdir <- "output/RData"

scene <- "BS_sp25_sim5_r4_L2"
load(file.path("output/RData/DataPrep", paste0(scene, ".RData")))
lapply(file.path(mdir, "Analysis", scene, dir),
       load,.GlobalEnv)

# load results of lefse
load("output/RData/lefse/out_lefse.RData")
# tree
phy_tree <- rowTree(tse)

br <- metadata(tse)$branch
br <- c(br$A, br$B)


# detected nodes
s <- 5
loc.treeclimbR <- lapply(outsel_0.05, FUN = function(x){
    x$node[x$signal.node]
}) 
loc.minP <- lapply(outMin, FUN = function(x){
    x$node[x$keep_0.05]
}) 
loc <- list("lasso" = loc.Lasso[[s]], 
            "HFDR" = loc.HFDR_0.05[[s]],
            "StructFDR" = loc.StructFDR_0.05[[s]],
            "BH" = loc.bh_0.05[[s]],
            "miLineage1" = loc1_0.05.MLA[[s]], 
            "miLineage2" = loc2_0.05.MLA[[s]],
            "minP" = loc.minP[[s]],
            "LEfSe" = out_lefse_low[[scene]]$`0.05`[[s]],
            "treeclimbR" = loc.treeclimbR[[s]])


p <- ggtree(phy_tree, branch.length = "none")
p <- groupClade(p, br) + aes(color = group) +
    scale_color_manual(values=c("grey", "orange", "blue")) 
p <- scaleClade(p, node = br[1], scale = 20)
p <- scaleClade(p, node = br[2], scale = 20)
p




fig_loc <- plot_grid(
    p + 
        geom_point2(aes(subset = (node %in% loc[[1]])), 
                    color = "red", size = 0.9) +
        theme(legend.position = "none") +
        expand_limits(y = 1770), 
    p + 
        geom_point2(aes(subset = (node %in% loc[[2]])), 
                    color = "red", size = 0.9) +
        theme(legend.position = "none") +
        expand_limits(y = 1770),
    p + 
        geom_point2(aes(subset = (node %in% loc[[3]])), 
                    color = "red", size = 0.6) +
        theme(legend.position = "none")+
        expand_limits(y = 1770),
    p + 
        geom_point2(aes(subset = (node %in% loc[[4]])), 
                    color = "red", size = 0.9) +
        theme(legend.position = "none"),
    p + 
        geom_point2(aes(subset = (node %in% loc[[5]])), 
                    color = "red", size = 0.9) +
        theme(legend.position = "none")+
        expand_limits(y = 1770),
    p + 
        geom_point2(aes(subset = (node %in% loc[[6]])), 
                    color = "red", size = 0.9) +
        theme(legend.position = "none")+
        expand_limits(y = 1770),
    p + 
        geom_point2(aes(subset = (node %in% loc[[7]])), 
                    color = "red", size = 0.9) +
        theme(legend.position = "none")+
        expand_limits(y = 1770),
    p + 
        geom_point2(aes(subset = (node %in% loc[[8]])), 
                    color = "red", size = 0.9) +
        theme(legend.position = "none") +
        expand_limits(y = 1770), 
    p + 
        geom_point2(aes(subset = (node %in% loc[[9]])), 
                    color = "red", size = 0.9) +
        theme(legend.position = "none") +
        expand_limits(y = 1770), 
    labels = names(loc),
    label_size = 8, 
    label_x = c(0, 0, 0, 0, -0.02, -0.02, 0, -0.08),
    nrow = 3)

