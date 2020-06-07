
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(RColorBrewer)
    library(treeclimbR)
    library(dplyr)
})


sim <- c("BS_sp10_sim5_r4_L2", 
         "BS_sp25_sim5_r4_L2", 
         "BS_sp50_sim5_r4_L2",
         "US_sp10_sim5_r4_L2", 
         "US_sp25_sim5_r4_L2", 
         "US_sp50_sim5_r4_L2",
         "SS_sp10_sim5_r4_pct0.4_L2", 
         "SS_sp25_sim5_r4_pct0.4_L2", 
         "SS_sp50_sim5_r4_pct0.4_L2")

# reso: high or low
# reso <- "high"
scenario <- rep(c("BS", "US", "SS"), each = 3)
repl <- rep(c(10, 25, 50), 3)
mdir <- "output/RData"
dir <- c("BH.RData", "treeclimbR.RData", "minP.RData",
         "HFDR.RData", "LassoGLM.RData", "miLineage.RData",
         "StructFDR.RData")

# remove ancestor
source("summary/rm_ancestor.R")

# calculate TPR/FDR function 
source("summary/rateFun.R")

# load results of lefse
load("output/RData/lefse/out_lefse.RData")

avDat <- vector("list", length(sim))
for (i in seq_along(avDat)) {
    cat(i, "\n")
    cat(scenario[i], "\n")
    cat(repl[i], "\n")
    si <- sim[i]
    
    
    load(file.path(mdir, "DataPrep", paste0(si, ".RData")))
    load(file.path(mdir, "parameter", paste0(si, ".RData")))
    lapply(file.path(mdir, "Analysis", si, dir),
           load,.GlobalEnv)
    
    # ==================================================================
    # the number of simulation
    nSim <- length(assays(tse))
    
    
    # truth
    fc <- metadata(tse)$FC
    truth <- names(fc[fc != 1])
    truth <- signalNode(tree = rowTree(tse), node = truth)
    
    # lefse
    if (reso == "high") {
        out_lefse <- out_lefse_high
    } else {
        out_lefse <- out_lefse_low
    }
    loc.lefse <- out_lefse[[si]]
    rate.lefse <- rateFun(loc.lefse, nSim, method = "LEfSe")
    
    # minP
    loc.minP <- lapply(c(0.01, 0.05, 0.1), FUN = function(x){
        xl <- lapply(outMin, FUN = function(y) {
            yy <- y$node[y[[paste0("keep_", x)]]]
            if (reso == "high") {
               yy <- rm_ancestor(node = yy, tree = rowTree(tse))
            }
            yy
        })
    })
    names(loc.minP) <- c(0.01, 0.05, 0.1)
    rate.minP <- rateFun(loc.minP, nSim, method = "minP")
    
    # treeclimbR
    loc.treeclimbR <- lapply(list(outsel_0.01, outsel_0.05, outsel_0.1),
                             FUN = function(x) {
        xl <- lapply(x, FUN = function(y) {
            yy <- y$node[y$signal.node]
            if (reso == "high") {
                yy <- rm_ancestor(node = yy, tree = rowTree(tse))
            }
            yy
        })
    })
    names(loc.treeclimbR) <- c(0.01, 0.05, 0.1)
    rate.treeclimbR <- rateFun(loc.treeclimbR, nSim, method = "treeclimbR")
    
    # miLineage part I
    loc1.MLA <- lapply(list(loc1_0.01.MLA, loc1_0.05.MLA, 
                            loc1_0.1.MLA), 
                       FUN = function(x) {
                           xl <- lapply(x, FUN = function(y) {
                               if (reso == "high") {
                                   y <- rm_ancestor(node = y, tree = rowTree(tse))
                               } 
                               return(y)
                           })
                        })
    names(loc1.MLA) <- c(0.01, 0.05, 0.1)
    rate.ml1 <- rateFun(loc1.MLA, nSim, method = "miLineage1")
    
    
    # miLineage part II
    loc2.MLA <- lapply(list(loc2_0.01.MLA, loc2_0.05.MLA, 
                            loc2_0.1.MLA), 
                       FUN = function(x) {
                           xl <- lapply(x, FUN = function(y) {
                               if (reso == "high") {
                                   y <- rm_ancestor(node = y, tree = rowTree(tse))
                               } 
                               return(y)
                           })
                       })
    names(loc2.MLA) <- c(0.01, 0.05, 0.1)
    rate.ml2 <- rateFun(loc2.MLA, nSim, method = "miLineage2")
    
    
    # lasso
    loc.Lasso <- lapply(loc.Lasso, FUN = function(x) {
        if (reso == "high") {
        x <- rm_ancestor(node = x, tree = rowTree(tse))
        } 
        return(x)
    })
    loc.ll <- list(loc.Lasso, loc.Lasso, loc.Lasso)
    names(loc.ll) <- c(0.01, 0.05, 0.1)
    rate.Lasso <- rateFun(loc.ll, 
                          nSim, method = "lasso")
    # StructFDR
    loc.StructFDR <- lapply(list(loc.StructFDR_0.01, loc.StructFDR_0.05, 
                          loc.StructFDR_0.1), FUN = function(x) {
                              xl <- lapply(x, FUN = function(y) {
                                  if (reso == "high") {
                                      y <- rm_ancestor(node = y, tree = rowTree(tse))
                                  } 
                                  return(y)
                              })
                          })
    names(loc.StructFDR) <- c(0.01, 0.05, 0.1)
    rate.StructFDR <- rateFun(loc.StructFDR, nSim, method = "StructFDR")
    
    # BH
    loc.bh <- lapply(list(loc.bh_0.01, loc.bh_0.05, 
                          loc.bh_0.1), FUN = function(x) {
        xl <- lapply(x, FUN = function(y) {
            if (reso == "high") {
                y <- rm_ancestor(node = y, tree = rowTree(tse))
            } 
            return(y)
        })
    })
    names(loc.bh) <- c(0.01, 0.05, 0.1)
    rate.bh <- rateFun(loc.bh, nSim, method = "BH")
    
    # HFDR
    loc.HFDR <- lapply(list(loc.HFDR_0.01, loc.HFDR_0.05, 
                          loc.HFDR_0.1), FUN = function(x) {
                              xl <- lapply(x, FUN = function(y) {
                                  if (reso == "high") {
                                      y <- rm_ancestor(node = y, tree = rowTree(tse))
                                  } 
                                  return(y)
                              })
                          })
    names(loc.HFDR) <- c(0.01, 0.05, 0.1)
    rate.HFDR <- rateFun(loc.HFDR, nSim, method = "HFDR")
    
    
    
    
    rateDat <- rbind(rate.treeclimbR, rate.minP, 
                     rate.ml1, rate.ml2,
                     rate.Lasso, rate.HFDR, 
                     rate.bh, rate.StructFDR,
                     rate.lefse)
    
    
    avDat[[i]] <- rateDat %>%
        group_by(method, alpha) %>% 
        summarise(fdr = mean(fdr), tpr = mean(tpr)) %>%
        mutate(scenario = scenario[i]) %>% 
        mutate(replicates = repl[i])
    
}


Dat <- do.call(rbind, avDat)

Dat <- Dat %>%
    ungroup() %>%
    mutate(pc = factor(scenario, levels = c("BS", "US", "SS"),
                       labels = c("BS", "US", "SS"))) %>%
    mutate(method = factor(method, 
                           levels = c("lasso", "HFDR",
                                      "StructFDR", "BH", 
                                      "miLineage1", 
                                      "miLineage2", "minP",
                                      "LEfSe",
                                      "treeclimbR"),
                           labels = c("lasso", "HFDR",
                                      "StructFDR", "BH", 
                                      "miLineage1", 
                                      "miLineage2", "minP", 
                                      "LEfSe",
                                      "treeclimbR"))) %>%
    arrange(method) %>%
    filter(!(method == "lasso" & alpha != 0.05)) 

vcolor <- c("treeclimbR" = "#E41A1C", "BH" = "#377EB8",
            "StructFDR" = "#4DAF4A", "HFDR" = "#984EA3",
            "lasso" = "#FF7F00", "minP" = "#A65628",
            "miLineage1" = "#999999", "miLineage2" = "#666666",
            "LEfSe" = "#E7298A")
vshape <- c("treeclimbR" = 5, "BH" = 16,
            "StructFDR" = 16, "HFDR" = 16,
            "lasso" = 16, "minP" = 16,
            "miLineage1" = 16, "miLineage2" = 16,
            "LEfSe" = 16)

vsize <- c("treeclimbR" = 3, "BH" = 2,
           "StructFDR" = 2, "HFDR" = 2,
           "lasso" = 2, "minP" = 2,
           "miLineage1" = 2, "miLineage2" = 2,
           "LEfSe" = 2)
prettify <- theme_bw(base_size = 8) + theme(
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
    legend.position="bottom",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin=margin(-5,-5,-10,-5),
    strip.background = element_rect(colour = "black", fill = "gray90"),
    strip.text.x = element_text(color = "black", size = 8),
    strip.text.y = element_text(color = "black", size = 8))



p_roc <- ggplot(Dat) +
    geom_vline(xintercept = c(0.01, 0.05, 0.1),
               linetype = "longdash", color = "black",
               size = 0.2) +
    geom_point(aes(x = fdr, y = tpr, color = method, 
                   shape = method, size = method)) +
    geom_path(aes(x = fdr, y = tpr,
                  color = method, group = method), 
              size = 0.3, show.legend = FALSE) +
    #    facet_wrap(~output) +
    xlab("FDR") +
    ylab("TPR") +
    facet_grid(rows = vars(pc), cols = vars(replicates)) +
    scale_color_manual("Method", values = vcolor) +
    scale_shape_manual("Method", values = vshape) +
    scale_size_manual("Method", values = vsize) +
    guides(color=guide_legend(nrow = 1, byrow = TRUE,
                              override.aes = list(size = 2.5))) +
    scale_x_continuous(trans = "sqrt", 
                       breaks = c(0.01, 0.05, 0.1, seq(0, 1, 0.2)),
                       limits = c(0, 1), labels = function(u) formatC(u)) +
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    prettify +
    theme(legend.box.margin=margin(-5,-5,0,-5))

p_roc
#save(p_roc, Dat , file = "summary/roc_L2.RData")


