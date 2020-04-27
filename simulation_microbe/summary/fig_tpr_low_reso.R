
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(RColorBrewer)
    library(treeclimbR)
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

scenario <- rep(c("BS", "US", "SS"), each = 3)
repl <- rep(c(10, 25, 50), 3)
mdir <- "output/RData"
dir <- c("BH.RData", "treeclimbR.RData", "minP.RData",
         "HFDR.RData", "LassoGLM.RData", "miLineage.RData",
         "StructFDR.RData")

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
            cbind.data.frame(rate[[x]], alpha = c(0.01, 0.05, 0.1), 
                             sim = x, method = method)
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
    
    loc.treeclimbR_0.01 <- lapply(outsel_0.01, FUN = function(x) {
        x$node[x$signal.node]
    })
    loc.treeclimbR_0.05 <- lapply(outsel_0.05, FUN = function(x) {
        x$node[x$signal.node]
    })
    loc.treeclimbR_0.1 <- lapply(outsel_0.1, FUN = function(x) {
        x$node[x$signal.node]
    })
    rate.treeclimbR <- rateFun(loc.treeclimbR_0.01, loc.treeclimbR_0.05, 
                         loc.treeclimbR_0.1, nSim, method = "treeclimbR")
    
    # miLineage part I
    rate.ml1 <- rateFun(loc1_0.01.MLA, loc1_0.05.MLA,
                        loc1_0.1.MLA, nSim, method = "miLineage1")


    # miLineage part II
    rate.ml2 <- rateFun(loc2_0.01.MLA, loc2_0.05.MLA,
                        loc2_0.1.MLA, nSim, method = "miLineage2")
    
    
    rate.Lasso <- rateFun(loc.Lasso, loc.Lasso, loc.Lasso, 
                          nSim, method = "lasso")
    # structFDR
    rate.StructFDR <- rateFun(loc.StructFDR_0.01, loc.StructFDR_0.05,
                        loc.StructFDR_0.1, nSim, method = "StructFDR")
    
    # BH
    rate.bh <- rateFun(loc.bh_0.01, loc.bh_0.05, 
                       loc.bh_0.1, nSim, method = "BH")
    
    # HFDR
    rate.HFDR<- rateFun(loc.HFDR_0.01, loc.HFDR_0.05,
                       loc.HFDR_0.1, nSim, method = "HFDR")
    
    
    
    
    rateDat <- rbind(rate.treeclimbR, rate.minP, 
                     rate.ml1, rate.ml2,
                     rate.Lasso, rate.HFDR, rate.bh, rate.StructFDR)
    
    
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
                                      "treeclimbR"),
                           labels = c("lasso", "HFDR",
                                      "StructFDR", "BH", 
                                      "miLineage1", 
                                      "miLineage2", "minP", 
                                      "treeclimbR"))) %>%
    arrange(method) %>%
    filter(!(method == "lasso" & alpha != 0.05)) 

vcolor <- c("treeclimbR" = "#E41A1C", "BH" = "#377EB8",
            "StructFDR" = "#4DAF4A", "HFDR" = "#984EA3",
            "lasso" = "#FF7F00", "minP" = "#A65628",
            "miLineage1" = "#999999", "miLineage2" = "#666666")
vshape <- c("treeclimbR" = 5, "BH" = 16,
            "StructFDR" = 16, "HFDR" = 16,
            "lasso" = 16, "minP" = 16,
            "miLineage1" = 16, "miLineage2" = 16)

vsize <- c("treeclimbR" = 3, "BH" = 2,
           "StructFDR" = 2, "HFDR" = 2,
           "lasso" = 2, "minP" = 2,
           "miLineage1" = 2, "miLineage2" = 2)
prettify <- theme_bw(base_size = 8) + theme(
    aspect.ratio = 1,
    #plot.margin = unit(c(0, 0, 0.2, 0), "cm"),
    panel.grid = element_blank(),
    panel.spacing = unit(0, "lines"),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.key.size= unit(1.5, "mm"),
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
save(p_roc, Dat , file = "summary/roc_L2.RData")


