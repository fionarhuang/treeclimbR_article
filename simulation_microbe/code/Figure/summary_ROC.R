

library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(treeAGG2)



sim <- c("S1_lib1E5_sp10_sim5_r4_L2", 
         "S1_lib1E5_sp25_sim5_r4_L2", 
         "S1_lib1E5_sp50_sim5_r4_L2",
         "S2_lib1E5_sp10_sim5_r4_L2", 
         "S2_lib1E5_sp25_sim5_r4_L2", 
         "S2_lib1E5_sp50_sim5_r4_L2",
         "S3_lib1E5_sp10_sim5_r4_pct0.2_L2", 
         "S3_lib1E5_sp25_sim5_r4_pct0.2_L2", 
         "S3_lib1E5_sp50_sim5_r4_pct0.2_L2")
# labOutput <- c("S1, 10 replicates", "S1, 25 replicates", "S1, 50 replicates",
#                "S2, 10 replicates", "S2, 25 replicates", "S2, 50 replicates",
#                "S3, 10 replicates", "S3, 25 replicates", "S3, 50 replicates") 
scenario <- rep(c("S1", "S2", "S3"), each = 3)
repl <- rep(c(10, 25, 50), 3)
mdir <- "/Volumes/fiona/phd/microbes/simulation/wflow_microbes_new/output/RData"
avDat <- vector("list", length(sim))
for (i in seq_along(avDat)) {
    cat(i, "\n")
    cat(scenario[i], "\n")
    cat(repl[i], "\n")
    si <- sim[i]

load(file.path(mdir, "DataPrep", paste0(si, ".RData")))
load(file.path(mdir, "parameter", paste0(si, ".RData")))
load(file.path(mdir, "Analysis", si, "BH.RData"))
load(file.path(mdir, "Analysis", si, "edgeR_aggP.RData"))
load(file.path(mdir, "Analysis", si, "edgeR_minP.RData"))
load(file.path(mdir, "Analysis", si, "hcFDR.RData"))
load(file.path(mdir, "Analysis", si, "LassoGLM.RData"))
load(file.path(mdir, "Analysis", si, "miLineage.RData"))
load(file.path(mdir, "Analysis", si, "structFDR.RData"))


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
    x$nodeNum[x$keep_0.01]
})
loc.minP_0.05 <- lapply(outMin, FUN = function(x) {
    x$nodeNum[x$keep_0.05]
})
loc.minP_0.1 <- lapply(outMin, FUN = function(x) {
    x$nodeNum[x$keep_0.1]
})
rate.minP <- rateFun(loc.minP_0.01, loc.minP_0.05,
                     loc.minP_0.1, nSim, method = "minP")

# aggP

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


avDat[[i]] <- rateDat %>%
    group_by(method, a) %>% 
    summarise(fdr = mean(fdr), tpr = mean(tpr)) %>%
    mutate(scenario = scenario[i]) %>% 
    mutate(replicates = repl[i])

}


Dat <- do.call(rbind, avDat)
#Dat$output <- paste(Dat$output, "replicates", sep = " ")
vshape <- c(1, 2, 4)
names(vshape) <- unique(Dat$a)
cbbPalette <-  c("#377EB8", "#F781BF", "#B2DF8A", "#33A02C", 
                 "#A65628", "#984EA3", "#FF7F00", "#E41A1C")
# vcolor <- cbbPalette[c(1, 2, 3, 7, 8)]
vcolor <- cbbPalette
names(vcolor) <- unique(Dat$method)


prettify <- theme_bw(base_size = 9) + theme(
    aspect.ratio = 1,
    #plot.margin = unit(c(0, 0, 0.2, 0), "cm"),
    panel.grid = element_blank(),
    panel.spacing = unit(0, "lines"),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.key.size= unit(2, "mm"),
    plot.title = element_text(hjust = 0.5),
    legend.text = element_text(size = 8),
    legend.position="bottom",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin=margin(-1,-1,-1,-1),
    strip.background = element_rect(colour = "black", fill = "gray90"),
    strip.text.x = element_text(color = "black", size = 8),
    strip.text.y = element_text(color = "black", size = 8, angle = 0)) 

p <- ggplot(Dat) +
    geom_vline(xintercept = c(0.01, 0.05, 0.1),
               linetype = "longdash", color = "black",
               size = 0.2) +
    geom_point(aes(x = fdr, y = tpr,
                    color = method)) +
    geom_path(aes(x = fdr, y = tpr,
                  color = method, group = method), 
              show.legend = FALSE) +
#    facet_wrap(~output) +
    facet_grid(rows = vars(scenario), cols = vars(replicates)) +
    xlab("FDR") +
    ylab("TPR") +
    scale_shape_manual(values = vshape) +
    scale_color_manual(values = vcolor) +
    guides(color=guide_legend(NULL, nrow = 1, byrow = TRUE,
                              override.aes = list(size = 3))) +
           #shape = guide_legend(title = "FDR level", ncol = 1)) +
    #    guides(color=guide_legend(nrow=2,byrow=TRUE)) +
    
    # ylim(yl[1], yl[2]) +
    scale_x_continuous(trans = "sqrt", 
                       breaks = c(0.01, 0.05, 0.1, seq(0, 1, 0.2)),
                       limits = c(0, 1), labels = function(u) formatC(u)) +
    scale_y_continuous(breaks = c(seq(0, 1, 0.2)), limits = c(0, 1)) +
    prettify 
p
ggsave("~/Documents/ms/microbes/microbes_roc.eps", p, 
       units = "cm", width = 15, height = 15,
       dpi = 300)
       #, useDingbats = FALSE)

