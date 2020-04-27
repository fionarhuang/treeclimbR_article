# === tool ========
source("/home/fiona/phd/microbes/simulation/Tool/argsR_compile.R")
.libPaths(c(Sys.getenv('R_LIBS_1'), 
            Sys.getenv('R_LIBS_2'),
            .libPaths()))
.libPaths()
suppressPackageStartupMessages({
    library(TreeSummarizedExperiment)
    library(treeclimbR)
    library(ggtree)
    library(ggplot2)
    library(cowplot)
    library(dplyr)
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
load(parSet)
outFig
ls()
# ====================== TPR/FDR function ===========================
rateFun <- function(loc0.01, loc0.05, loc0.1, nSim, method) {
    res <- vector("list", nSim)
    for (i in seq_len(nSim)) {
        res[[i]] <- list(loc_0.01 = loc0.01[[i]],
                         loc_0.05 = loc0.05[[i]],
                         loc_0.1 = loc0.1[[i]])
    }
    rate <- lapply(res, FUN = function(x) {
        xx <- lapply(x, function(y) {
            fdr.y <- fdr(tree = rowTree(tse), truth = truth, 
                         found = y, only.leaf = TRUE)
            tpr.y <- tpr(tree = rowTree(tse), truth = truth, 
                         found = y, only.leaf = TRUE)
            c(fdr.y, tpr.y)
        })
        do.call(rbind, xx)
    })
    
    rate <- lapply(seq_along(rate), FUN = function(x) {
        cbind.data.frame(rate[[x]], a = c(0.01, 0.05, 0.1), sim = x, method = method)
    })
    rate <- do.call(rbind, rate)
    rownames(rate) <- NULL
    return(rate)
}

# ==================================================================
# the number of simulation
nSim <- length(assays(tse))
rejLimit <- 0.05
# truth
fc <- metadata(tse)$FC
truth <- names(fc[fc != 1])
truth <- signalNode(tree = rowTree(tse), node = truth)

# minP
loc.minP_0.01 <- lapply(outMin, FUN = function(x) {
    x$node[x$keep_0.01]
})
loc.minP_0.05 <- lapply(outMin, FUN = function(x) {
    x$node[x$keep_0.05]
})
loc.minP_0.1 <- lapply(outMin, FUN = function(x) {
    x$node[x$keep_0.1]
})
rate.minP <- rateFun(loc.minP_0.01, loc.minP_0.05,
                     loc.minP_0.1, nSim, method = "minP")

# treeclimbR

loc.aggP_0.01 <- lapply(outsel_0.01, FUN = function(x) {
    x$node[x$signal.node]
})
loc.aggP_0.05 <- lapply(outsel_0.05, FUN = function(x) {
    x$node[x$signal.node]
})
loc.aggP_0.1 <- lapply(outsel_0.1, FUN = function(x) {
    x$node[x$signal.node]
})
rate.aggP <- rateFun(loc.aggP_0.01, loc.aggP_0.05, 
                     loc.aggP_0.1, nSim, method = "aggFDR")

# miLineage part I
rate.ml1 <- rateFun(loc1_0.01.MLA, loc1_0.05.MLA, 
                    loc1_0.1.MLA, nSim, method = "miL1")


# miLineage part II
rate.ml2 <- rateFun(loc2_0.01.MLA, loc2_0.05.MLA, 
                    loc2_0.1.MLA, nSim, method = "miL2")


# Lasso has no level (e.g. 0.01, 0.05, 0.1) to specified

# rate.Lasso <- lapply(loc.Lasso, FUN = function(x) {
#     x <- as.numeric(x)
#     fdr.x <- fdr(tree = rowTree(tse), truth = truth, 
#                  found = x, only.leaf = TRUE)
#     tpr.x <- tpr(tree = rowTree(tse), truth = truth,
#                  found = x, only.leaf = TRUE)
#     c(fdr.x, tpr.x)
# })
# rate.Lasso <- do.call(rbind, rate.Lasso)
# rate.Lasso <- rate.Lasso %>% 
#     data.frame %>%
#     mutate(sim = 1:5)
rate.Lasso <- rateFun(loc.Lasso, loc.Lasso, loc.Lasso, 
                      nSim, method = "Lasso")
# structFDR
rate.str <- rateFun(loc.str_0.01, loc.str_0.05,
                    loc.str_0.1, nSim, method = "strFDR")

# BH
rate.bh <- rateFun(loc.bh_0.01, loc.bh_0.05, 
                   loc.bh_0.1, nSim, method = "BH")

# hcFDR
rate.hc <- rateFun(loc.hcFDR_0.01, loc.hcFDR_0.05,
                   loc.hcFDR_0.1, nSim, method = "hcFDR")




rateDat <- rbind(rate.aggP, rate.minP, rate.ml1, rate.ml2,
                 rate.Lasso, rate.hc, rate.bh, rate.str)
avDat <- rateDat %>%
    group_by(method, a) %>% 
    summarise(fdr = mean(fdr), tpr = mean(tpr)) 

# title <- paste(scene, paste("lib", muNB, sep = "-"), 
#                paste("sp-", nSam[1], sep = ""), 
#                paste("ratio-", metadata(dat)$branch$ratio, sep = ""), 
#                sep = ",")


vshape <- c(4, 5, 16, 16, 17, 2, 10, 1)

cbbPalette <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                "#66A61E", "#E6AB02", "#A6761D", "#666666")
# vcolor <- cbbPalette[c(1, 2, 3, 7, 8)]
vcolor <- cbbPalette
names(vshape) <- names(vcolor) <- unique(rateDat$Methods)

# if (ncol(tse) == 100) {
#     yl <- c(0.7, 1.01)
#     xl <- c(-0.01, 1.01)
# } else {
    yl <- c(0, 1)
    xl <- c(0, 1)
# }

# p1 <- ggplot(rateDat) +
#     geom_jitter(aes(x = fdr, y = tpr,
#                     color = method, shape = method),
#                 size = 3, alpha = 0.8, stroke = 1.5,
#                 width = 0.007, height = 0.007) +
#     xlab("FDR") +
#     ylab("TPR") +
#     theme(strip.text = element_text(size = 25),
#           axis.text = element_text(size = 18),
#           axis.title = element_text(size = 25),
#           legend.title = element_text(size = 18),
#           legend.text = element_text(size = 18),
#           legend.position="bottom") +
#     scale_shape_manual(values = vshape) +
#     scale_color_manual(values = vcolor) +
#     guides(shape = guide_legend(nrow = 2)) +
#     #    guides(color=guide_legend(nrow=2,byrow=TRUE)) +
#     geom_vline(xintercept = 0.05, linetype = "dashed", color = "darkviolet") +
#     ylim(yl[1], yl[2]) +
#     xlim(xl[1], xl[2])

p2 <- ggplot(avDat) +
    geom_point(aes(x = fdr, y = tpr,
                    color = method, shape = method),
                size = 3, alpha = 0.8, stroke = 1.5) +
    geom_path(aes(x = fdr, y = tpr,
                  color = method, group = method), size = 3, alpha = 0.3) +
    xlab("FDR") +
    ylab("TPR") +
    theme(strip.text = element_text(size = 15),
          axis.text = element_text(size = 18),
          axis.title = element_text(size = 18),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.position="right") +
    scale_shape_manual(values = vshape) +
    scale_color_manual(values = vcolor) +
    guides(shape = guide_legend(ncol = 1)) +
    #    guides(color=guide_legend(nrow=2,byrow=TRUE)) +
    geom_vline(xintercept = c(0.01, 0.05, 0.1),
               linetype = "dashed", color = "darkviolet") +
    ylim(yl[1], yl[2]) +
    xlim(xl[1], xl[2]) + 
    ggtitle(paste0(scene, "\n (", nSam[1], " replicates in each group)")) +
    theme(plot.title = element_text(hjust = 0.5))

# p2 <- ggplot(rateDat) +
#     geom_boxplot(aes(x = factor(Methods), y = signalBr)) +
#     geom_jitter(aes(x = Methods, y = signalBr, color = Methods)) +
#     xlab("") +
#     ylab("The number of found branches") +
#     theme(strip.text = element_text(size = 25),
#           axis.text.y = element_text(size = 15),
#           axis.title.y = element_text(size = 15),
#           axis.text.x = element_text(angle = 90, size = 13)) +
#     scale_y_log10(limits = c(-1, 100)) +
#     geom_hline(yintercept = length(truth), color = "blue", linetype = "dashed") +
#     guides(color = FALSE)
#plot_grid(p1, p2, nrow = 2)
png(outFig)
p2
dev.off()

# outAv <- gsub(pattern = ".png", replacement = ".pdf", out)
# pdf(outAv)
# p2
# dev.off()
# =============for later publication ==========================
# bitmap("Plot.tiff", height = 4, width = 4, units = 'in', type="tifflzw", res=300)
# plot(x, y)
# dev.off()


# pdf(out)
# par(mar = c(5, 5, 4, 1))
# plot(rr.minP[, 1], rr.minP[, 2], xlim = c(0, 1), ylim = c(0, 1), 
#      col = rgb(0, 0.5, 0.1, 0.5), xlab = "FDR", ylab = "TPR",
#      type = "n", cex.axis = 2, cex.lab = 2, main = title)
# points(rr.minP[, 1], rr.minP[, 2], col = rgb(0, 0.5, 0.1, 0.7), pch = 15, cex = 2)
# points(rr.MLA[, 1], rr.MLA[, 2], col = rgb(0, 0, 1, 0.7), pch = 16, cex = 2)
# points(rr100.MLA[, 1], rr100.MLA[, 2], col = rgb(0, 0.2, 1, 0.7), pch = 17, cex = 2)
# points(rr.Lasso[, 1], rr.Lasso[, 2], col = rgb(1, 0, 1, 0.7), pch = 18, cex = 2)
# legend("topright", pch = 15:18, 
#        col = c(rgb(0, 0.5, 0.1, 0.7), 
#                rgb(0, 0, 1, 0.7), 
#                rgb(0, 0.2, 1, 0.7),
#                rgb(1, 0, 1, 0.7)), pt.cex=2,
#        legend = c("minP", "miLineage", "miLineage100", "Lasso"))
# dev.off()


