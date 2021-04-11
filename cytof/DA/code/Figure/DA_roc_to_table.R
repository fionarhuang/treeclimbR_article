

library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(treeclimbR)
library(kableExtra)
library(tidyr)

# ---------------------- load data ----------------------
numClus <- rep(c(400, 900, 1600), 3)
pc <- rep(c("0.1pc", "1pc", "5pc"), each = 3)
pcs <- rep(c("weak", "medium", "strong"), each = 3)
seed <- 2020
sim <- paste0("CN_", pc, "_clus", numClus, "_seed", seed)
#sim <- paste0("CN_", pc, "_clus", numClus)
condition <- c("healthy", "CN")

mdir <- "DA/output/RData"
avDat <- vector("list", length(sim))
for (i in seq_along(avDat)) {
    cat(i, "\n")
    cat(numClus[i], "\n")
    si <- sim[i]
    
    load(file.path(mdir, "DataPrep", paste0(si, ".RData")))
    load(file.path(mdir, "parameter", paste0(si, ".RData")))
    load(file.path(mdir, "Analysis", si, "diffcyt.RData"))
    load(file.path(mdir, "Analysis", si, "treeclimbR.RData"))
    load(file.path(mdir, "Analysis", si, "minP.RData"))
    load(file.path(mdir, "Analysis", si, "HFDR.RData"))
    load(file.path(mdir, "Analysis", si, "LassoGLM.RData"))
    load(file.path(mdir, "Analysis", si, "StructFDR.RData"))
    
    # ====================== TPR/FDR function ===========================
    rateFun <- function(loc, method, tree, 
                        condition, df_truth) {
        
        locList <- lapply(loc, FUN = function(x) {
            if (length(x)) {
                if(is.factor(x)) {
                    x <- as.character(x)
                }
                
                xx <- findOS(tree = tree, node = x, 
                             only.leaf = TRUE, self.include = TRUE)
                unlist(unique(xx))  
            } else {
                NULL
            }
            
        })
        
        
        type_column <- ifelse("CN" %in% condition, "n_CN", "n_CBF")
        
        # discovery
        sel <- lapply(locList, FUN = function(x) {
            df_truth$node %in% x
        })
        
        n.pos <- lapply(sel, FUN = function(x) {
            sum(df_truth[x, ][["n_healthy"]], df_truth[x, ][[type_column]])
        })
        
        # truth
        n.truth <- sum(df_truth[[type_column]])
        
        # true discovery
        n.tp <- lapply(sel, FUN = function(x) {
            sum(df_truth[x, ][[type_column]])})
        
        tpr <- mapply(FUN = function(x, y){x/max(c(y, 1))}, n.tp, n.truth)
        
        # false discovery
        n.fp <- mapply(FUN = function(x, y){x - y}, n.pos, n.tp)
        fdr <- mapply(FUN = function(x, y){x /max(c(y, 1))}, n.fp, n.pos)
        
        out <- cbind.data.frame(tpr = tpr, 
                                fdr = fdr,
                                alpha = names(loc),
                                method = method)
        return(out)
    }
    
    # ==================================================================
    # the number of simulation
    nSim <- length(assays(tse))
    nSam <- unique(table(colData(tse)$group))
    Tree <- rowTree(tse)
    
    
    # minP
    loc.minP <- list(outMin$node[outMin$keep_0.01],
                     outMin$node[outMin$keep_0.05],
                     outMin$node[outMin$keep_0.1])
    names(loc.minP) <- c("0.01", "0.05", "0.1")
    rate.minP <- rateFun(loc.minP,  method = "minP", 
                         df_truth = df_truth,
                         tree = Tree, condition = condition)
    
    # treeclimbR
    rate.treeclimbR <- rateFun(loc.treeclimbR,  method = "treeclimbR", 
                               df_truth = df_truth,
                               tree = Tree, condition = condition)
    
    
    # diffcyt
    rate.diffcyt <- rateFun(loc.diffcyt,  method = "diffcyt", 
                            df_truth = df_truth,
                            tree = Tree, condition = condition)
    
    
    
    
    # Lasso has no level (e.g. 0.01, 0.05, 0.1) to specified
    loc.Lasso <- c(loc.Lasso,
                   loc.Lasso,
                   loc.Lasso)
    names(loc.Lasso) <- c("0.01", "0.05", "0.1")
    rate.Lasso <- rateFun(loc.Lasso,  method = "lasso", 
                          df_truth = df_truth,
                          tree = Tree, condition = condition)
    # structFDR
    loc.StructFDR <- c(loc.str_0.01,
                       loc.str_0.05,
                       loc.str_0.1)
    names(loc.StructFDR) <- c("0.01", "0.05", "0.1")
    rate.StructFDR <- rateFun(loc.StructFDR,  method = "StructFDR", 
                              df_truth = df_truth,
                              tree = Tree, condition = condition)
    
    
    # HFDR
    
    rate.HFDR <- rateFun(loc.hc,  method = "HFDR", 
                         df_truth = df_truth,
                         tree = Tree, condition = condition)
    
    
    rateDat <- rbind(rate.treeclimbR, rate.minP, rate.diffcyt,
                     rate.Lasso, rate.HFDR, rate.StructFDR)
    avDat[[i]] <- rateDat %>%
        group_by(method, alpha) %>% 
        summarise(fdr = mean(fdr), tpr = mean(tpr)) %>%
        mutate(pc = pcs[i]) %>% 
        mutate(cluster = numClus[i])
    
}


Dat_DA <- do.call(rbind, avDat)
head(Dat_DA)

Dat_tab <- Dat_DA %>%
    relocate(alpha, .before = method) %>%
    arrange(alpha, pc, method, cluster)%>%
    mutate(tpr = round(tpr, digits = 3),
           fdr = round(fdr, digits = 3)) %>%
    pivot_wider(names_from = pc, values_from = c(tpr, fdr)) %>%
    pivot_wider(names_from = cluster, 
                values_from = c(tpr_weak, fdr_weak, 
                                tpr_medium, fdr_medium,
                                tpr_strong, fdr_strong)) 
colnames(Dat_tab) <- gsub(".*_", "", colnames(Dat_tab))



kbl(Dat_tab, 
    format = "latex", 
    booktabs = T) %>%
    kable_classic(full_width = F) %>%
    kable_styling(latex_options = c("scale_down")) %>%
    column_spec(1, bold = T) %>%
    collapse_rows(columns = 1, valign = "top") %>%
    add_header_above(c(" ", " ", "TPR" = 3, "FDR" = 3, 
                       "TPR" = 3, "FDR" = 3,
                       "TPR" = 3, "FDR" = 3)) %>%
    add_header_above(c(" ", " ", "Weak" = 6, "Medium" = 6, "Strong" = 6)) 
     




