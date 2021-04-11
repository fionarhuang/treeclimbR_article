
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(RColorBrewer)
    library(treeclimbR)
    library(dplyr)
    library(stringr)
})



# reso: high or low
reso <- "high"


mdir <- "output/RData"
dir <- c("BH.RData", "treeclimbR.RData", "minP.RData",
         "HFDR.RData", "LassoGLM.RData", 
         "miLineage.RData",
         "StructFDR.RData")
(sim <- list.files(file.path(mdir, "Analysis")))
repl <- gsub(pattern = ".*_sp", "", sim)
scenario <- gsub(pattern = "_.*", "", sim)


# remove ancestor
source("summary/rm_ancestor.R")

# calculate TPR/FDR function 
source("summary/rateFun.R")

# load results of lefse
load("output/RData/lefse/out_lefse.RData")

avDat <- vector("list", length(sim))
for (i in seq_along(avDat)) {
    message(sim[i])
    si <- sim[i]
    
    
    load(file.path(mdir, "DataPrep", paste0(si, ".RData")))
    load(file.path(mdir, "parameter", paste0(si, ".RData")))
    lapply(file.path(mdir, "Analysis", si, dir),
           load,.GlobalEnv)
    
    # ==================================================================
    # the number of simulation
    nSim <- length(tse_list)
    
    # tree
    Tree <- rowTree(tse_list[[1]])
    
    # truth
    fc <- lapply(tse_list, FUN = function(x) {metadata(x)$FC})
    truth <- lapply(fc, FUN = function(x) {
        xx <- names(x[x != 1])
        signalNode(tree = Tree, node = xx)
    })
    
    
    
    # lefse
    if (reso == "high") {
        out_lefse <- out_lefse_high
    } else {
        out_lefse <- out_lefse_low
    }
    loc.lefse <- out_lefse[[si]]
    rate.lefse <- rateFun(loc = loc.lefse, nSim = nSim,
                          method = "LEfSe", tree = Tree,
                          truth = truth)
    
    
    
    # minP
    loc.minP <- lapply(c(0.01, 0.05, 0.1), FUN = function(x){
        xl <- lapply(outMin, FUN = function(y) {
            yy <- y$node[y[[paste0("keep_", x)]]]
            if (reso == "high") {
                yy <- rm_ancestor(node = yy, tree = Tree)
            }
            yy
        })
    })
    names(loc.minP) <- c(0.01, 0.05, 0.1)
    rate.minP <- rateFun(loc = loc.minP, nSim = nSim, 
                         method = "minP", tree = Tree, 
                         truth = truth)
    
    # treeclimbR
    loc.treeclimbR <- lapply(list(outsel_0.01, outsel_0.05, outsel_0.1),
                             FUN = function(x) {
                                 xl <- lapply(x, FUN = function(y) {
                                     yy <- y$node[y$signal.node]
                                     if (reso == "high") {
                                         yy <- rm_ancestor(node = yy, tree = Tree)
                                     }
                                     yy
                                 })
                             })
    names(loc.treeclimbR) <- c(0.01, 0.05, 0.1)
    rate.treeclimbR <- rateFun(loc = loc.treeclimbR, nSim = nSim, 
                               method = "treeclimbR", tree = Tree, 
                               truth = truth)
    
    # miLineage part I
    loc1.MLA <- lapply(list(loc1_0.01.MLA, loc1_0.05.MLA,
                            loc1_0.1.MLA),
                       FUN = function(x) {
                           xl <- lapply(x, FUN = function(y) {
                               if (reso == "high") {
                                   y <- rm_ancestor(node = y, tree = Tree)
                               }
                               return(y)
                           })
                       })
    names(loc1.MLA) <- c(0.01, 0.05, 0.1)
    rate.ml1 <- rateFun(loc = loc1.MLA, nSim = nSim,
                        method = "miLineage1", tree = Tree,
                        truth = truth)
    
    
    # miLineage part II
    loc2.MLA <- lapply(list(loc2_0.01.MLA, loc2_0.05.MLA,
                            loc2_0.1.MLA),
                       FUN = function(x) {
                           xl <- lapply(x, FUN = function(y) {
                               if (reso == "high") {
                                   y <- rm_ancestor(node = y, tree = Tree)
                               }
                               return(y)
                           })
                       })
    names(loc2.MLA) <- c(0.01, 0.05, 0.1)
    rate.ml2 <- rateFun(loc = loc2.MLA, nSim = nSim,
                        method = "miLineage2", tree = Tree,
                        truth = truth)
    
    
    # lasso
    loc.Lasso <- lapply(loc.Lasso, FUN = function(x) {
        if (reso == "high") {
            x <- rm_ancestor(node = x, tree = Tree)
        } 
        return(x)
    })
    loc.ll <- list(loc.Lasso, loc.Lasso, loc.Lasso)
    names(loc.ll) <- c(0.01, 0.05, 0.1)
    rate.Lasso <- rateFun(loc = loc.ll, nSim = nSim,
                          method = "lasso", tree = Tree, 
                          truth = truth)
    # StructFDR
    loc.StructFDR <- lapply(list(loc.StructFDR_0.01, loc.StructFDR_0.05, 
                                 loc.StructFDR_0.1), FUN = function(x) {
                                     xl <- lapply(x, FUN = function(y) {
                                         if (reso == "high") {
                                             y <- rm_ancestor(node = y, tree = Tree)
                                         } 
                                         return(y)
                                     })
                                 })
    names(loc.StructFDR) <- c(0.01, 0.05, 0.1)
    rate.StructFDR <- rateFun(loc = loc.StructFDR, nSim = nSim,
                              method = "StructFDR", tree = Tree, 
                              truth = truth)
    
    # BH
    loc.bh <- lapply(list(loc.bh_0.01, loc.bh_0.05, 
                          loc.bh_0.1), FUN = function(x) {
                              xl <- lapply(x, FUN = function(y) {
                                  if (reso == "high") {
                                      y <- rm_ancestor(node = y, tree = Tree)
                                  } 
                                  return(y)
                              })
                          })
    names(loc.bh) <- c(0.01, 0.05, 0.1)
    rate.bh <- rateFun(loc = loc.bh, nSim = nSim,
                       method = "BH", tree = Tree, 
                       truth = truth)
    
    # HFDR
    loc.HFDR <- lapply(list(loc.HFDR_0.01, loc.HFDR_0.05, 
                            loc.HFDR_0.1), FUN = function(x) {
                                xl <- lapply(x, FUN = function(y) {
                                    if (reso == "high") {
                                        y <- rm_ancestor(node = y, tree = Tree)
                                    } 
                                    return(y)
                                })
                            })
    names(loc.HFDR) <- c(0.01, 0.05, 0.1)
    rate.HFDR <- rateFun(loc = loc.HFDR, nSim = nSim,
                         method = "HFDR", tree = Tree, 
                         truth = truth)
    
    
    
    
    rateDat <- rbind(rate.treeclimbR, rate.minP, 
                     rate.ml1, rate.ml2,
                     rate.Lasso, rate.HFDR,
                     rate.lefse, 
                     rate.bh, rate.StructFDR)
    
    
    avDat[[i]] <- rateDat %>%
        group_by(method, alpha) %>% 
        mutate(scenario = scenario[i]) %>% 
        mutate(replicates = repl[i])
    
}

avDat <- lapply(seq_along(avDat), FUN = function(x) {
    avDat[[x]][, "scenario"] <- scenario[x]
    avDat[[x]]
})

rbDat <- do.call(rbind, avDat)

Dat <- rbDat %>%
    group_by(method, alpha, scenario, replicates) %>% 
    summarise(fdr = mean(fdr), tpr = mean(tpr)) %>%
    mutate(pc = factor(scenario, levels = c("BS", "US", "SS"),
                       labels = c("BS", "US", "SS")),
           replicates = as.numeric(replicates)) %>%
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
    filter(!(method == "lasso" & alpha != 0.05))%>%
    ungroup() 

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
    #    facet_wrap(~ pc)+
    facet_grid(rows = vars(pc), cols = vars(replicates)) +
    scale_color_manual("Method", values = vcolor) +
    scale_shape_manual("Method", values = vshape) +
    scale_size_manual("Method", values = vsize) +
    guides(color=guide_legend(ncol = 1, byrow = TRUE,
                              override.aes = list(size = 2.5))) +
    scale_x_continuous(trans = "sqrt", 
                       breaks = c(0.01, 0.05, 0.1, seq(0, 1, 0.2)),
                       limits = c(0, 1), labels = function(u) formatC(u)) +
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    prettify +
    theme(legend.position = "left",
          legend.margin = margin(0, 10, 0, 0))

p_roc

path <- file.path("summary/figure/Supplementary_fdr_tpr_resoH_100.eps")
ggsave(path, units = "in", width = 6, height = 6,
       dpi = 300)

library(tidyr)
head(rbDat)
rbDat <- rbDat %>%
    ungroup() %>%
    mutate(pc = factor(scenario, levels = c("BS", "US", "SS"),
                       labels = c("BS", "US", "SS"))) %>%
    filter(method %in% c("BH", "minP", "StructFDR", "treeclimbR")) %>%
    pivot_longer(cols = c("fdr", "tpr"),
                 names_to = "rate",
                 values_to = "value")  %>%
    mutate(rate = toupper(rate))



ggplot(rbDat %>%
           filter(replicates == "25")) +
    geom_boxplot(aes(x = method, y = value, fill = rate),
                 outlier.size = 1,
                 outlier.stroke = 0.2,
                 notchwidth = 0.2) +
    scale_color_manual("Method", values = vcolor) +
    facet_grid(cols = vars(alpha), rows = vars(pc)) +
    labs(x = "Method", y = "", fill = "",
         subtitle ="25 samples per group") +
    theme(
        aspect.ratio = 1,
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 7),
        legend.position = "right",
        legend.background = element_rect(fill = "NA"),
        strip.background = element_rect(colour = "black", fill = "gray90"),
        strip.text.x = element_text(color = "black", size = 8),
        strip.text.y = element_text(color = "black", size = 8))
path <- file.path("summary/figure/Supplementary_fdr_tpr_sp25_boxplot.eps")
ggsave(path, units = "in", width = 8, height = 8,
       dpi = 300)


ggplot(rbDat %>%
           filter(replicates == "153")) +
    geom_boxplot(aes(x = method, y = value, fill = rate),
                 outlier.size = 1,
                 outlier.stroke = 0.2,
                 notchwidth = 0.2) +
    scale_color_manual("Method", values = vcolor) +
    facet_grid(cols = vars(alpha), rows = vars(pc)) +
    labs(x = "Method", y = "", fill = "",
         subtitle = "153 samples per group") +
    theme(
        aspect.ratio = 1,
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 7),
        legend.position = "right",
        legend.background = element_rect(fill = "NA"),
        strip.background = element_rect(colour = "black", fill = "gray90"),
        strip.text.x = element_text(color = "black", size = 8),
        strip.text.y = element_text(color = "black", size = 8))

path <- file.path("summary/figure/Supplementary_fdr_tpr_sp153_boxplot.eps")
ggsave(path, units = "in", width = 8, height = 8,
       dpi = 300)

# Fold change
df <- lapply(sim, FUN = function(x) {
    message(x)
    load(file.path(mdir, "DataPrep", paste0(x, ".RData")))
    
    fc <- sapply(tse_list, FUN = function(y) {
        FC <- metadata(y)$FC
        unique(abs(log2(FC[FC != 1])))
    })
    
    rep <- lapply(seq_along(fc), FUN = function(x) {
        rep(x, length(fc[[x]]))
    })
    
    data.frame(scene = gsub("_.*", "", x),
               sp = gsub(".*sp", "", x),
               absLFC = unlist(fc),
               repetition = unlist(rep))
})

DF <- do.call(rbind, df) %>%
    mutate(scene = factor(scene, levels = c("BS", "US", "SS")),
           sp = as.numeric(sp))

ggplot(DF) +   
    geom_point(aes(x = repetition, y = absLFC), size = 0.6)+
    labs(y = "abs(logFC)") +
    facet_grid(rows = vars(scene), cols = vars(sp))
head(DF)

path <- file.path("summary/figure/Supplementary_rep100_fc.eps")
ggsave(path, units = "in", width = 8, height = 5,
       dpi = 300)


# branch size
df <- lapply(sim, FUN = function(x) {
    message(x)
    load(file.path(mdir, "DataPrep", paste0(x, ".RData")))
    
    nleaf <- lapply(tse_list, FUN = function(y) {
        br <- c(metadata(y)$br$A, metadata(y)$br$B)
        ds <- findOS(tree = rowTree(y), node = br)
        ds <- lapply(ds, transNode, tree = rowTree(y))
        
        FC <- metadata(y)$FC
        FC <- FC[FC != 1]
        
        rr <- names(FC)
        nn <- lapply(ds, FUN = function(x) {
            length(intersect(rr, x))
        })
        nn <- unlist(nn)
    })
    nleaf <- do.call("rbind", nleaf)
    data.frame(scene = gsub("_.*", "", x),
               sp = gsub(".*sp", "", x),
               nleaf_A = nleaf[, 1],
               nleaf_B = nleaf[, 2],
               repetition = seq_len(nrow(nleaf)))
})


DF <- do.call(rbind, df) %>%
    mutate(scene = factor(scene, levels = c("BS", "US", "SS")),
           sp = as.numeric(sp))


ggplot(DF ) +   
    geom_point(aes(x = repetition, y = nleaf_A), size = 0.6,
               shape = 3)+
    geom_point(aes(x = repetition, y = nleaf_B), size = 0.6)+
    geom_segment(aes(x = repetition, xend = repetition, 
                     y = nleaf_A, yend = nleaf_B), 
                 color = "red") +
    labs(y = "Branch sizes") +
#    facet_wrap(~ scene, nrow = 3) +
    facet_grid(rows = vars(scene), cols = vars(sp)) +
    scale_y_log10(breaks = c(0, 10, 50, 200, 600))

path <- file.path("summary/figure/Supplementary_rep100_br.eps")
ggsave(path, units = "in", width = 8, height = 5,
       dpi = 300)





# Relative abundance
data("throat_v35")
throat_v35 <- parEstimate(obj = throat_v35)
abd <- metadata(throat_v35)$assays.par$pi
df <- lapply(sim, FUN = function(x) {
    message(x)
    load(file.path(mdir, "DataPrep", paste0(x, ".RData")))
    
    ly <- lapply(tse_list, FUN = function(y) {
        FC <- metadata(y)$FC
        namU <- names(FC[FC>1])
        namL <- names(FC[FC<1])
        xc <- c(sum(abd[namU]), sum(abd[namL]))
        sort(xc)
    })
    ly <- do.call("rbind", ly)
    data.frame(scene = gsub("_.*", "", x),
               sp = gsub(".*sp", "", x),
               br1 = ly[, 1],
               br2 = ly[, 2],
               repetition = seq_len(nrow(ly)))
})


DF <- do.call(rbind, df) %>%
    mutate(scene = factor(scene, levels = c("BS", "US", "SS")),
           sp = as.numeric(sp))


ggplot(DF ) +   
    geom_point(aes(x = repetition, y = br1), size = 0.6,
               shape = 3)+
    geom_point(aes(x = repetition, y = br2), size = 0.6)+
    geom_segment(aes(x = repetition, xend = repetition, 
                     y = br1, yend = br2), 
                 color = "red") +
    labs(y = "Relative abundance") +
#    facet_wrap(~ scene, nrow = 3) +
    facet_grid(rows = vars(scene), cols = vars(sp)) +
    scale_y_log10(breaks = c(0, 0.001, 0.005, 0.01, 0.025, 0.1, 0.5)) 

path <- file.path("summary/figure/Supplementary_rep100_abd.eps")
ggsave(path, units = "in", width = 8, height = 5,
       dpi = 300)


