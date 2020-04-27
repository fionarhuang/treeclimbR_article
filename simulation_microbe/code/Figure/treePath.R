
# === tool ========
source("/home/fiona/phd/microbes/simulation/Tool/argsR_compile.R")

library(TreeSummarizedExperiment)
library(treeAGG2)
library(dplyr)
#library(RColorBrewer)
library(ggplot2)
library(cowplot)
# =========== source TreeInfer package ===================
R.Version()


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
print(inRDat1)
print(inRDat2)
load(inRDat1)
load(inRDat2)

tree <- rowTree(tse)
dat <- lapply(outAgg_sel, FUN = function(x) {
    
    mm <- cbind(tree$Nnode + 2, x$node)
    # dist <- apply(mm, 1, function(x) {
    #     distNode(tree = tree, node = x) })
    # sc <- x %>%
    #     select(starts_with("treeP_")) 
    
    x <- x %>%
    #    mutate(dist = dist) %>%
        mutate(tp = 1-PValue) %>%
        mutate(score = ifelse(logFC > 0, tp, -tp))
    
})

fc <- metadata(tse)$FC
truth_up <- names(fc[fc > 1])
truth_down <- names(fc[fc < 1])
truth_up <- signalNode(tree = rowTree(tse), node = truth_up)
truth_down <- signalNode(tree = rowTree(tse), node = truth_down)
truth <- c(truth_down, truth_up)

lab <- rep(c("down", "up"), c(length(truth_down), length(truth_up)))
ll <- length(unique(lab))
# clab <- brewer.pal(2, name = "Dark2")
colr <- c("blue", "orange")[seq_len(ll)]
names(colr) <- unique(lab)

figList <- vector("list", length(outAgg_sel))
figLeg <- vector("list", 1)
for (i in seq_along(outAgg_sel)) {
    dat.i <- dat[[i]]
    title.i <- paste("simulation:", i)
    found <- dat.i$node[dat.i$signal.node]
    p <- viewPath(tree = tree, data_frame = dat.i, 
                             node_column = "node", 
                             y_axis_column = "score", 
 #                            x_axis_column = "dist", 
                             node_true = c(truth_down, truth_up), 
                             lab_true = lab, col_true = colr, 
                             node_found = found, x_lab = "N_leaf",
                             log_x = TRUE, log_y = FALSE) +
        ggtitle(title.i) +
        theme(plot.title = element_text(hjust = 0.5,
                                        face = "bold",
                                        size = 12))
    figList[[i]] <- p
    # figLeg[[1]] <- get_legend(p)
    # figList[[i]] <- p + theme(legend.position = "none")
    
}

figAll <- c(figList, figLeg)
pdf(out)
#plot_grid(plotlist = figAll, nrow = 4)
figAll
dev.off()
