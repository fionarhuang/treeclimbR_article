

library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(treeclimbR)


numClus <- rep(c(100, 400, 900), 3)
pc <- rep(c("less_distinct_75pc", "less_distinct_50pc", "main"), each = 3)
pcs <- rep(c("weak", "medium", "strong"), each = 3)
seed <- 2020
sim <- paste0(pc, "_clus", numClus, "_seed", seed)

sim
mdir <- "DS/output/RData"
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
    load(file.path(mdir, "Analysis", si, "LassoGLM.RData"))
    

# ====================== TPR/FDR function ===========================
rateFun <- function(loc, method, tree, df_truth) {
    
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
    
    locList <- lapply(locList, FUN = function(x) {
        if (is.numeric(x)) {
            x <- transNode(tree = tree, node = x, use.alias = FALSE)
        } else {
            x
        }
        return(x)}
    )
    
    
    
    # discovery
    sel <- lapply(locList, FUN = function(x) {
        df_truth$cluster_id %in% x
    })
    
    n.pos <- lapply(sel, FUN = function(x) {
        sum(df_truth[x, ]$n_cells)
    })
    
    
    prop_column <- "prop_B"
    
    # truth
    n.truth <- sum(df_truth$n_cells * df_truth[[prop_column]])
    
    # true discovery
    n.tp <- lapply(sel, FUN = function(x) {
        sum(df_truth[x, ]$n_cells * df_truth[x, prop_column])})
    
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
Tree <- rowTree(d_medians_all)

loc.treeclimbR <- lapply(loc.treeclimbR, unique)
rate.treeclimbR <- rateFun(loc = loc.treeclimbR,  method = "treeclimbR", 
                           df_truth = df_truth,
                           tree = Tree)


# diffcyt
loc.diffcyt <- lapply(loc.diffcyt, unique)
rate.diffcyt <- rateFun(loc.diffcyt ,  method = "diffcyt", 
                        df_truth = df_truth, tree = Tree)

# Lasso
loc.Lasso <- replicate(3, unique(loc.Lasso), simplify = FALSE)
names(loc.Lasso) <- c("0.01", "0.05", "0.1")

rate.Lasso <- rateFun(loc = loc.Lasso,  method = "lasso", 
                      df_truth = df_truth,
                      tree = Tree)

loc.minP <- lapply(loc.minP, unique)
rate.minP <- rateFun(loc = loc.minP,  method = "minP", 
                      df_truth = df_truth,
                      tree = Tree)

rateDat <- rbind(rate.treeclimbR, rate.diffcyt, rate.Lasso, rate.minP)
avDat[[i]] <- rateDat %>%
    group_by(method, alpha) %>% 
    summarise(fdr = mean(fdr), tpr = mean(tpr)) %>%
    mutate(cluster = numClus[i]) %>%
    mutate(pc = pcs[i]) %>%
    data.frame()

}


Dat_DS <- do.call(rbind, avDat)
Dat_DS <- Dat_DS %>%
    ungroup() %>%
    mutate(pc = factor(pc, levels = c("weak", "medium", "strong"))) %>%
    mutate(method = factor(method,
                           levels = c("lasso", "diffcyt",
                                      "minP", "treeclimbR"))) %>%
    arrange(method) %>%
    dplyr::filter(!(method == "lasso" & alpha != 0.05))  

set.seed(1)
rf <- runif(min = -0.02, max = 0.02, n = nrow(Dat_DS))


vshape <- c("lasso" = 16, "minP" = 16, "diffcyt" = 16, 
            "treeclimbR" = 5)

vcolor <- c("treeclimbR" = "#E41A1C", "diffcyt" = "#377EB8",
             "lasso" = "#FF7F00", "minP" = "#A65628")
vsize <- c("treeclimbR" = 3, "diffcyt" = 2,
           "lasso" = 2, "minP" = 2)



library(ggplot2)
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

#set.seed(1)
p_DS <- ggplot(Dat_DS) +
    geom_vline(xintercept = c(0.01, 0.05, 0.1),
               linetype = "longdash", color = "black",
               size = 0.2) +
    geom_point(aes(x = fdr, y = tpr + rf, color = method, 
                   shape = method, size = method),
               stroke = 0.5) +
    geom_path(aes(x = fdr, y = tpr + rf,
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

p_DS
# outP <- file.path(sprintf("summary/DS_roc_%d.RData", seed))
# save(Dat_DS, file = outP)

# figPath <- file.path(sprintf("~/Documents/ms/cytof/DS_roc_%d.eps", seed))
# ggsave(figPath, p_DS, units = "cm", width = 9, height = 10,
#        dpi = 300)


