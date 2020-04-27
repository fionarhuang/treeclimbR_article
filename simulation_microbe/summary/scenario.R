
# .libPaths("/Users/ruizhu/Library/R/3.6/devel")
suppressPackageStartupMessages({
    library(treeclimbR)
    library(ggtree)
    library(ggplot2)
    library(cowplot)
})
# data("exTree")
set.seed(1)
exTree <- ape::rtree(50)
leaf1 <- findOS(tree = exTree, node = 68, only.leaf = TRUE)
leaf2 <- findOS(tree = exTree, node = 79, only.leaf = TRUE)
gr_leaf <- c(leaf1, leaf2)
names(gr_leaf) <- c("A", "B")
leaf1 <- unlist(leaf1)
leaf2 <- unlist(leaf2)
l1 <- length(leaf1)
l2 <- length(leaf2)
leaf3 <- c(38, 36, 34, 30, 28)


fig1 <- viewBranch(tree = exTree, group_leaf = gr_leaf, 
                   group_color = c("0" = "grey", "A" = "#D69C4E", "B" = "#046C9A"), 
                   zoom_node = c(58, 89, 52), zoom_scale = 1/20) + 
    geom_point2(aes(subset = (node %in% leaf1)), 
                color = "#D69C4E", size =  2) + 
    geom_point2(aes(subset = (node %in% leaf2)), 
                color = "#046C9A", size = 2) + 
   # geom_text(aes(x = 0.1, y = 18, label = "BS"), color = "black") +
    coord_flip() + scale_x_reverse() + theme(legend.position = "none")

set.seed(4)
fig2 <- viewBranch(tree = exTree, group_leaf = gr_leaf, 
                   group_color = c("0" = "grey", "A" = "#D69C4E", "B" = "#046C9A"), 
                   zoom_node = c(58, 89, 52), zoom_scale = 1/20) + 
    geom_point2(aes(subset = (node %in% leaf1)), 
                color = "#D69C4E", size =  runif(l1)*3) + 
    geom_point2(aes(subset = (node %in% leaf2)), 
                color = "#046C9A", size = runif(l2)*3) + 
   # geom_text(aes(x = 0.1, y = 18), label = "US", color = "black") + 
    coord_flip() + scale_x_reverse()+ theme(legend.position = "none")

gr3 <- as.list(leaf3)
names(gr3) <- rep("A", length(gr3))

fig3 <- viewBranch(tree = exTree, group_leaf = c(gr3, "B" = list(leaf1)), 
                   group_color = c("0" = "grey", "A" = "#046C9A", "B" = "#D69C4E"), 
                   zoom_node = c(58, 89, 52), zoom_scale = 1/20, 
                   ladderize = FALSE) + 
    geom_point2(aes(subset = (node %in% leaf1)), 
                color = "#D69C4E", size =  2) + 
    geom_point2(aes(subset = (node %in% leaf3)), 
                color = "#046C9A", size = 2)+ 
   # geom_text(aes(x = 0.1, y = 18, label = "SS"), color = "black") +
    coord_flip() + 
    scale_x_reverse() + 
    theme(legend.position = "none")

save(fig1, fig2, fig3, file = "summary/scenario.RData")

