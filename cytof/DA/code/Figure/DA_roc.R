

library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(treeclimbR)

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
Dat_DA <- Dat_DA %>%
    ungroup() %>%
    mutate(pc = factor(pc, levels = c("weak", "medium", "strong"))) %>%
    mutate(method = factor(method,
                           levels = c("lasso", "HFDR", "StructFDR", 
                                      "diffcyt", "minP", "treeclimbR"))) %>%
    arrange(method) %>%
    dplyr::filter(!(method == "lasso" & alpha != 0.05)) 





vshape <- c("lasso" = 16, "minP" = 16, "StructFDR" = 16, "diffcyt" = 16, 
            "treeclimbR" = 5, "HFDR" = 16)

vcolor <- c("treeclimbR" = "#E41A1C", "diffcyt" = "#377EB8",
            "StructFDR" = "#4DAF4A", "HFDR" = "#984EA3",
            "lasso" = "#FF7F00", "minP" = "#A65628")
vsize <- c("treeclimbR" = 3, "diffcyt" = 2,
           "StructFDR" = 2, "HFDR" = 2,
           "lasso" = 2, "minP" = 2)

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



set.seed(2)
rf <- runif(min = -0.02, max = 0.02, n = nrow(Dat_DA))
rf <- 0
p_DA <- ggplot(Dat_DA) +
    geom_vline(xintercept = c(0.01, 0.05, 0.1),
               linetype = "longdash", color = "black",
               size = 0.2) +
    # geom_jitter(aes(x = fdr, y = tpr, color = method),
    #             size = 1, height = 0.025, width = 0) +
    geom_point(aes(x = fdr, y = tpr, color = method, 
                   shape = method, size = method),
               stroke = 0.5) +
    geom_path(aes(x = fdr, y = tpr,
                  color = method, group = method), 
              size = 0.3, show.legend = FALSE) +
    #    facet_wrap(~output) +
    xlab("FDR") +
    ylab("TPR") +
    facet_grid(rows = vars(cluster), cols = vars(pc)) +
    scale_color_manual("Method", values = vcolor) +
    scale_shape_manual("Method", values = vshape) +
    scale_size_manual("Method", values = vsize) +
    guides(color=guide_legend(nrow = 1, byrow = TRUE,
                              override.aes = list(size = 2))) +
    scale_x_continuous(trans = "sqrt", 
                       breaks = c(0.01, 0.05, 0.1, seq(0, 1, 0.2)),
                       limits = c(0, 1), labels = function(u) formatC(u)) +
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    prettify 

p_DA
# outP <- file.path(sprintf("summary/DA_roc_%d.RData", seed))
# save(Dat_DA, file = outP)




# ggsave(figPath, p_DA, units = "cm", width = 9, height = 10,
#        dpi = 300)

