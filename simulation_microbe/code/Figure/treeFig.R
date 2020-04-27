
# === tool ========
source("/home/fiona/phd/microbes/simulation/Tool/argsR_compile.R")
.libPaths(c(Sys.getenv('R_LIBS_1'), 
            Sys.getenv('R_LIBS_2'),
            .libPaths()))
.libPaths()

suppressPackageStartupMessages({
    library(treeclimbR)
    library(ggtree)
    library(ggplot2)
    library(cowplot)
})



# ==== arguments from batch R=====================
args <- (commandArgs(trailingOnly = TRUE))
args
argsList <- argRun(args, grp.pattern = "dirGRP")
argsList
for (i in seq_along(argsList)) {
    if(length(argsList[[i]])>1){
        assign(names(argsList)[i],argsList[[i]])
    }else{eval(parse(text = argsList[[i]]))}
}

## === load data from data preparation step ======
lapply(dirGRP, load, .GlobalEnv)
load(inRDat)
print(outFig)

# truth
fc <- metadata(tse)$FC
truth <- names(fc[fc != 1])
truth <- transNode(tree = rowTree(tse), node = truth)

br <- c(metadata(tse)$br$A, metadata(tse)$br$B)
leafBr <- findOS(tree = rowTree(tse), 
                 node = br)

# truth in branch A
truthA <- intersect(truth, leafBr[[1]])
truthB <- intersect(truth, leafBr[[2]])


# the number of simulation
n <- length(assays(tse))

# ----------------------- 0.05 ------------------------------------
# minP
loc.minP_0.05 <- lapply(outMin, 
                        FUN = function(x) {
                            x$node[x$keep_0.05]
                        })
# aggP
loc.aggP_0.05 <- lapply(outsel_0.05,
                        FUN = function(x){
                            xx <- x$node[x$signal.node]
                            xx
                        })

# structFDR
loc.str_0.05

# BH
loc.bh_0.05

# hcFDR
loc.hcFDR_0.05

# miLineage
loc1_0.05.MLA 
loc2_0.05.MLA

# LASSO
loc.Lasso

loc_0.05 <- list(loc.aggP_0.05, loc.minP_0.05, loc.str_0.05, loc.bh_0.05,
                 loc.hcFDR_0.05, loc1_0.05.MLA, loc2_0.05.MLA, loc.Lasso)
methods <- c("aggFDR", "minP", "strFDR", "BH", "hcFDR",
             "miL1", "miL2", "Lasso")
names(loc_0.05) <- methods

hf <- c("orange", "blue")[seq_along(br)]
edSize <- 0.2
psize <- 2
zscale <- 50
for (i in seq_len(n)){
    figList <- vector("list", length(methods))
    for (j in seq_along(methods)) {
      cat(j, "\n")
        figList[[j]] <- viewBranch(tree = rowTree(tse), 
                                   hlight_node = br, 
                                   hlight_fill = hf,
                                   hlight_alpha = 0.2,
                                   group_leaf = list(A = truthA, B = truthB),
                                   group_color = c("A" = "red", "B" = "red", "0" = "black"),
                                   point_node = loc_0.05[[j]][[i]], 
                                   point_color = "#1B9E77", 
                                   point_size = psize, 
                                   layout = "rectangular",
                                   edge_size = edSize, 
                                   zoom_node = br, 
                                   zoom_scale = zscale) +
            ggtitle(methods[j]) + 
            theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15)) 
    }
    
    otf <- gsub(pattern = "Fig.txt",
                replacement = paste0("Fig", i,"_0.05.png"), outFig)
    png(otf, height = 600, bg = NA)
    p <- plot_grid(plotlist = figList, nrow = 2)
    print(p)
    dev.off()
}

# ----------------------- 0.01 ------------------------------------
# minP
loc.minP_0.01 <- lapply(outMin, 
                        FUN = function(x) {
                            x$node[x$keep_0.01]
                        })
# aggP
loc.aggP_0.01 <- lapply(outsel_0.01,
                        FUN = function(x){
                            xx <- x$node[x$signal.node]
                            xx
                        })

# structFDR
loc.str_0.01

# BH
loc.bh_0.01

# hcFDR
loc.hcFDR_0.01

# miLineage
loc1_0.01.MLA 
loc2_0.01.MLA

# LASSO
loc.Lasso

loc_0.01 <- list(loc.aggP_0.01, loc.minP_0.01, loc.str_0.01, loc.bh_0.01,
                 loc.hcFDR_0.01, loc1_0.01.MLA, loc2_0.01.MLA, loc.Lasso)
methods <- c("aggFDR", "minP", "strFDR", "BH", "hcFDR",
             "miL1", "miL2", "Lasso")
names(loc_0.01) <- methods

hf <- c("orange", "blue")[seq_along(br)]
edSize <- 0.2
psize <- 2
zscale <- 50
for (i in seq_len(n)){
    figList <- vector("list", length(methods))
    for (j in seq_along(methods)) {
        figList[[j]] <- viewBranch(tree = rowTree(tse), 
                                   hlight_node = br, 
                                   hlight_fill = hf,
                                   hlight_alpha = 0.2,
                                   group_leaf = list(A = truthA, B = truthB),
                                   group_color = c("A" = "red", "B" = "red", "0" = "black"),
                                   point_node = loc_0.01[[j]][[i]], 
                                   point_color = "#1B9E77", 
                                   point_size = psize, 
                                   layout = "rectangular",
                                   edge_size = edSize, 
                                   zoom_node = br, 
                                   zoom_scale = zscale) +
            ggtitle(methods[j]) + 
            theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15)) 
    }
    
    otf <- gsub(pattern = "Fig.txt", 
                replacement = paste0("Fig", i,"_0.01.png"), outFig)
    png(otf, height = 600, bg = NA)
    p <- plot_grid(plotlist = figList, nrow = 2)
    print(p)
    dev.off()  
}


# ----------------------- 0.1 ------------------------------------
# minP
loc.minP_0.1 <- lapply(outMin, 
                        FUN = function(x) {
                            x$node[x$keep_0.1]
                        })
# aggP
loc.aggP_0.1 <- lapply(outsel_0.1,
                        FUN = function(x){
                            xx <- x$node[x$signal.node]
                            xx
                        })

# structFDR
loc.str_0.1

# BH
loc.bh_0.1

# hcFDR
loc.hcFDR_0.1

# miLineage
loc1_0.1.MLA 
loc2_0.1.MLA

# LASSO
loc.Lasso

loc_0.1 <- list(loc.aggP_0.1, loc.minP_0.1, loc.str_0.1, loc.bh_0.1,
                 loc.hcFDR_0.1, loc1_0.1.MLA, loc2_0.1.MLA, loc.Lasso)
methods <- c("aggFDR", "minP", "strFDR", "BH", "hcFDR",
             "miL1", "miL2", "Lasso")
names(loc_0.1) <- methods

hf <- c("orange", "blue")[seq_along(br)]
edSize <- 0.2
psize <- 3
zscale <- 50
for (i in seq_len(n)){
    figList <- vector("list", length(methods))
    for (j in seq_along(methods)) {
        figList[[j]] <- viewBranch(tree = rowTree(tse), 
                                   hlight_node = br, 
                                   hlight_fill = hf,
                                   hlight_alpha = 0.2,
                                   group_leaf = list(A = truthA, B = truthB),
                                   group_color = c("A" = "red", "B" = "red", "0" = "black"),
                                   point_node = loc_0.1[[j]][[i]], 
                                   point_color = "#1B9E77", 
                                   point_size = psize, 
                                   layout = "rectangular",
                                   edge_size = edSize, 
                                   zoom_node = br, 
                                   zoom_scale = zscale) +
            ggtitle(methods[j]) + 
            theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15)) 
    }
    
    otf <- gsub(pattern = "Fig.txt", 
                replacement = paste0("Fig", i,"_0.1.png"), outFig)
    png(otf, height = 600, bg = NA)
    p <- plot_grid(plotlist = figList, nrow = 2)
    print(p)
    dev.off()  
}
capture.output(writeLines("success"), file = outFig)

