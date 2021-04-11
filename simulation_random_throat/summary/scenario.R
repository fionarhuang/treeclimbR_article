
suppressPackageStartupMessages({
    library(treeclimbR)
    library(ggtree)
    library(ggplot2)
    library(cowplot)
    library(dplyr)
})


set.seed(1)
exTree <- ape::rtree(50)
br_color <- findOS(tree = exTree, node = c(68, 79), only.leaf = FALSE)
leaf_color <- findOS(tree = exTree, node = c(68, 79), only.leaf = TRUE)
ll <- lapply(leaf_color, length)
df_a <- data.frame(node = showNode(tree = exTree, only.leaf = FALSE)) %>%
    mutate(group = ifelse(node %in% br_color[[1]], "A", 
                         ifelse(node %in% br_color[[2]], "B", "other")))

fig0_a <- ggtree(exTree, branch.length = "none") %<+% df_a +
    aes(color = group) +
    scale_color_manual(values = c("other" = "grey", 
                                  "A" = "#D69C4E", "B" = "#046C9A")) 
fig0_a <- fig0_a %>% 
    scaleClade(node = 58, scale = 1/20) %>%
    scaleClade(node = 89, scale = 1/20) %>%
    scaleClade(node = 52, scale = 1/20) 


df_b <- data.frame(node = showNode(tree = exTree, only.leaf = FALSE)) %>%
    mutate(group = ifelse(node %in% br_color[[1]], "A", 
                          ifelse(node %in% c(38, 36, 34, 30, 28), "B", "other")))
fig0_b <- ggtree(exTree, branch.length = "none") %<+% df_b +
    aes(color = group) +
    scale_color_manual(values = c("other" = "grey", 
                                  "A" = "#D69C4E", "B" = "#046C9A")) 
fig0_b <- fig0_b %>% 
    scaleClade(node = 58, scale = 1/20) %>%
    scaleClade(node = 89, scale = 1/20) %>%
    scaleClade(node = 52, scale = 1/20) 


## ----------------------------BS----------------------------------------------
fig1 <- fig0_a + 
    geom_point2(aes(subset = (node %in% leaf_color[[1]])), 
                color = "#D69C4E", size =  2) + 
    geom_point2(aes(subset = (node %in% leaf_color[[2]])), 
                color = "#046C9A", size = 2) + 
    coord_flip() + scale_x_reverse() + theme(legend.position = "none")



## ----------------------------US----------------------------------------------
set.seed(4)
fig2 <- fig0_a + 
    geom_point2(aes(subset = (node %in% leaf_color[[1]])), 
                color = "#D69C4E", size =  runif(ll[[1]])*3) + 
    geom_point2(aes(subset = (node %in% leaf_color[[2]])), 
                color = "#046C9A", size = runif(ll[[2]])*3) + 
    coord_flip() + scale_x_reverse()+ theme(legend.position = "none")


## ----------------------------SS----------------------------------------------
fig3 <- fig0_b + 
    geom_point2(aes(subset = (node %in% leaf_color[[1]])), 
                color = "#D69C4E", size =  2) + 
    geom_point2(aes(subset = (node %in% c(38, 36, 34, 30, 28))), 
                color = "#046C9A", size = 2)+ 
   # geom_text(aes(x = 0.1, y = 18, label = "SS"), color = "black") +
    coord_flip() + 
    scale_x_reverse() + 
    theme(legend.position = "none")



