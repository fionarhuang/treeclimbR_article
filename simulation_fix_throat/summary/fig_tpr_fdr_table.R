
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(RColorBrewer)
    library(treeclimbR)
    library(dplyr)
    library(kableExtra)
    library(tidyr)
})

reso <- "high"
sim <- list.files("output/RData/Analysis/")


# reso: high or low
# reso <- "high"
scenario <- gsub(pattern = "_.*", "", sim)
repl <- gsub(pattern = ".*_sp", "", sim)

repl <- as.numeric(repl)
mdir <- "output/RData"
dir <- list.files(file.path(mdir, "Analysis", sim[1]))

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
    dplyr::rename(pc = scenario) %>%
    relocate(alpha, .before = method) %>%
    mutate(pc = factor(pc, levels = c("BS", "US", "SS"),
                       labels = c("BS", "US", "SS"))) 

Dat <- Dat %>%
    mutate(alpha = as.numeric(alpha),
           tpr = round(tpr, digits = 3),
           fdr = round(fdr, digits = 3)) %>%
    arrange(alpha, pc, method, replicates) %>%
    pivot_wider(names_from = pc, values_from = c(tpr, fdr)) %>%
    pivot_wider(names_from = replicates, 
                values_from = c(tpr_BS, fdr_BS, 
                                tpr_US, fdr_US,
                                tpr_SS, fdr_SS))  
head(Dat)
colnames(Dat) <- gsub(".*_", "", colnames(Dat))
kbl(Dat, 
   format = "latex",
    booktabs = T) %>%
    kable_classic(full_width = F) %>%
    kable_styling(latex_options = c("scale_down")) %>%
    column_spec(1, bold = T) %>%
    collapse_rows(columns = 1, valign = "top") %>%
    add_header_above(c(" ", " ", "TPR" = 3, "FDR" = 3, 
                       "TPR" = 3, "FDR" = 3,
                       "TPR" = 3, "FDR" = 3)) %>%
    add_header_above(c(" ", " ", "BS" = 6, "US" = 6, "SS" = 6)) 





