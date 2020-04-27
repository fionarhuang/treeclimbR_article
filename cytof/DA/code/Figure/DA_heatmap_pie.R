
suppressPackageStartupMessages({
    library(TreeHeatmap)
    library(treeclimbR)
    library(ggtree)
    library(ggplot2)
    library(dplyr)
    library(TreeSummarizedExperiment)
    library(ggnewscale)
    library(tidyr)
    library(scales)
    library(ggimage)
})

# ---------------------- load data ----------------------
seed <- 2020
nclus <- 900
md <- "DA/output/RData"
sel <- sprintf("CN_1pc_clus%d_seed%d", nclus, seed)
ana <- c("treeclimbR", "diffcyt", "LassoGLM", "minP", "StructFDR", "HFDR")

dataFile <- file.path(md, "DataPrep", paste0(sel, ".RData"))
analysisFile <- paste0(md, "/Analysis/", sel, "/", ana, ".RData")
lapply(c(dataFile, analysisFile), load, .GlobalEnv)

# --------------------------------- result -------------------------------------
loc_da <- list(treeclimbR = loc.treeclimbR$`0.05`,
               StructFDR = unlist(loc.str_0.05),
               HFDR = loc.hc$`0.05`,
               lasso = unlist(loc.Lasso),
               diffcyt = loc.diffcyt$`0.05`,
               minP = outMin$node[outMin$keep_0.05])
loc_da <- lapply(loc_da, FUN = function(x){
    if(length(x)) {x} else{ NULL }})

# ---------------------- visualization level ----------------------
# tree
treeR <- rowTree(tse)

# true CN cell (threshold > 0.1)
ts <- 0.5
cluster_CN <- df_truth$cluster_id[df_truth$prop_CN > ts]
node_CN <- signalNode(tree = treeR, 
                      node = cluster_CN)
leaf_CN <- findOS(tree = treeR, node = node_CN, 
                   only.leaf = TRUE, self.include = TRUE)
names(leaf_CN) <- rep("CN", length(leaf_CN))

# not CN cell
leaf_A <- setdiff(df_truth$cluster_id, cluster_CN)
node_A <- signalNode(tree = treeR, node = leaf_A)
inner_A <- node_A[!isLeaf(tree = treeR, node = node_A)]


# visualize level in non-B branch
vis_A <- c(node_A)
comb_A <- vis_A[!isLeaf(tree = treeR, node =vis_A)]

# visualize level
vis_L <- c(node_A, 
           transNode(treeR, as.character(cluster_CN)))



# --------------------------------- tree ---------------------------------------
treeS <- ggtree(treeR, size = 0.4)
br_CN <- findOS(tree = treeR, node = node_CN, 
                only.leaf = FALSE, self.include = TRUE)
df_CN <- treeS$data %>%
    select(node) %>%
    mutate(Type = ifelse(node %in% unlist(br_CN), "AML", "others"))
treeS <- treeS %<+% df_CN +
    aes(color = Type) + 
    scale_color_manual(values = c("others" = "grey50", "AML" = "mediumvioletred"), 
                       labels = c("NO", "Yes"), guide = FALSE) +
    new_scale_color()


for (i in seq_along(comb_A)) {
    treeS <- treeS %>% collapse(node = comb_A[i])
}
treeS + geom_text2(aes(label = node))

treeS$data <- treeS$data %>%
    mutate(isTip = ifelse(node %in% vis_A, TRUE, isTip)) %>%
    mutate(labs = ifelse(isTip, NA, label)) %>%
    # this is to give more space to present signal branch
    mutate(x = ifelse(node %in% c(924), x*0.7, x),
           x = ifelse(node %in% 1072, x*0.9, x),
           x = ifelse(node %in% 1175, x*1.1, x)) 
treeFig <- treeS + 
    geom_tiplab(align = TRUE, aes(label = labs), show.legend = FALSE)  +
    xlim(0, 5) +
    ylim(c(0, 10))
    
    
treeFig

#---------------------- Truth ------------------------------------------------
# calculate the number of AML cells for each node in visualized level
leaf_L <- findOS(tree = treeR, node = vis_L, 
                 only.leaf = TRUE, self.include = TRUE)
CN_truth <- lapply(seq_along(leaf_L), FUN = function(x) {
    xx <- transNode(tree = treeR, node = leaf_L[[x]], use.alias = FALSE)
    dx <- df_truth %>%
        filter(cluster_id %in% xx) %>%
        mutate(node = vis_L[x],
               CN_n = n_cells*prop_CN,
               nonCN_n = n_cells*(1-prop_CN)) %>%
        summarise(node = unique(node),
                  CN_n = sum(CN_n),
                  nonCN_n = sum(nonCN_n))
    
})
CN_truth <- do.call(rbind, CN_truth)
c <- 2
tree_bar_gap <- 0.6
bar_point_gap <- 0.5
CN_df <- CN_truth %>%
    left_join(treeFig$data) %>%
    select(node, x, y, CN_n, nonCN_n) %>%
    mutate(AML = round(CN_n/(CN_n + nonCN_n), digits = 6),
           others = AML) %>%
    gather(Cell_type, prop, AML:others) %>%
    mutate(x_start = ifelse(Cell_type == "AML", 0, prop*c)+ max(x)+tree_bar_gap,
           x_end = ifelse(Cell_type == "AML", prop*c, 1*c)+ max(x)+tree_bar_gap,
           Total_cell = CN_n+nonCN_n)
fig_0 <- treeFig +
    xlim(c(0, 10))+
    ylim(c(0, 10))+
    geom_segment(data = . %>%
                     mutate(xend = max(x, na.rm = TRUE),
                            x = min(x, na.rm = TRUE),
                            y = max(y, na.rm = TRUE)) %>%
                     select(x, xend, y) %>%
                     distinct(),
                 aes(x = x, xend = xend, y = y + 1.5, yend = y + 1.5),
                 inherit.aes = FALSE, color = "black") +
    annotate("text", x = 0.5*(min(treeFig$data$x, na.rm = TRUE) + 
                                  max(treeFig$data$x, na.rm = TRUE)), 
             y = max(treeFig$data$y, na.rm = TRUE) + 1.8,
             label = "Tree + detected nodes", size = 2)+
    geom_rect(data = CN_df, 
              aes(xmin = x_start, xmax = x_end, 
                  ymin = y - 0.45, ymax = y + 0.45,
                  fill = Cell_type),
              inherit.aes = FALSE) +
    annotate("text", x = 0.5*(min(CN_df$x_start) + max(CN_df$x_end)),  
             y = max(CN_df$y) + 0.9,
             label = "Cell \n type", size = 2) +
    geom_point(data = CN_df %>% filter(Cell_type == "others"),
               aes(x = x_end + bar_point_gap, y = y, size = Total_cell),
               inherit.aes = FALSE) +
    geom_segment(data = CN_df %>%
                     mutate(x_start = min(x_start),
                            x_end = max(x_end) + bar_point_gap,
                            y = max(y)) %>%
                     select(x_start, x_end, y) %>%
                     distinct(), 
                 aes(x = x_start, xend = x_end , y = y + 1.5,
                     yend = y + 1.5), inherit.aes = FALSE) +
    annotate("text", x = 0.5*(min(CN_df$x_start) + max(CN_df$x_end+ bar_point_gap)),  
             y = max(CN_df$y) + 1.8,
             label = "Truth", size = 2)+
    scale_fill_manual(values = c("AML" = "mediumvioletred", "others" = "grey"),
                      guide = guide_legend(order = 1)) +
    scale_size(breaks = c(100, 500, 1E4, 5E4),
               range = c(0.2, 4), guide = guide_legend(order = 2)) +
    new_scale_fill()
fig_0

#---------------------- heatmap ------------------------------------------------
# get the order of tips
ro <- treeFig$data %>%
    dplyr::filter(!is.na(y) & isTip) %>%
    arrange(desc(y)) %>%
    select(label) %>%
    unlist()

# Heatmap
med <- assays(tse)[[1]][rowLinks(tse)$nodeLab %in% ro, ]
med <- apply(med, 2, FUN = function(x) {
    xx <- ifelse(x == 0, x+1, x)
    log2(xx)
})


# column split & labels
ann_c <- gsub(pattern = "_.*", "", colnames(med))
names(ann_c) <- colnames(med)
lab_c <- ifelse(ann_c == "CN", "Diseased", "Healthy")
names(lab_c) <- ann_c

# tree & heatmap
fig_1 <- TreeHeatmap(tree = treeR, tree_fig = fig_0, 
                     hm_data = med, rel_width = 0.8,
                     tree_hm_gap = 3.8, 
                     column_split = ann_c, 
                     column_split_label = lab_c, 
                     column_split_gap = 0.1, 
                     legend_title_hm = "Count (log)",
                     split_label_size = 2,
                     split_label_offset_y = 1,
                     split_label_fontface = "plain" ) +
    scale_fill_gradient(low = "navy", high = "yellow",
                        limits = c(0, 6),
                        breaks=scales::pretty_breaks(n=3),
                        oob = squish)+
    guides(fill = guide_legend(order = 3,
                               override.aes = list(size = 2.5)))  +
    new_scale_fill() 


fig_1 + 
    xlim(c(0, 15))+
    ylim(c(0, 10))
df_1 <- getData(tree_hm = fig_1, type = "heatmap")%>%
    mutate(x1 = min(x)- width/2, y1 = max(y), x2 = max(x)+width/2) %>%
    select(x1, x2, y1) %>%
    distinct()
fig_1 <- fig_1 +
    geom_segment(data = df_1, 
                 aes(x = x1, xend = x2, y = y1 + 1.5,
                     yend = y1 + 1.5), inherit.aes = FALSE) +
    geom_text(data = df_1,
              aes(x = 0.5*(x1+x2), y = y1 + 1.8),
              inherit.aes = FALSE, label = "Observations", size = 2)
fig_1 +
    xlim(c(0, 15))+
    ylim(c(0, 10))

# --------------------------- annotate heatmap rows ---------------------------

df_row <- lapply(loc_da, FUN = function(y){
    if (length(y)) {
        yy <- unlist(findOS(tree = treeR, node = y, 
                            only.leaf = FALSE, 
                            self.include = TRUE))
        ty <- transNode(tree = treeR, node = yy, use.alias = TRUE)
    } else {
        ty <- NULL
    }
    
    ry <- rownames(med) %in% ty
    names(ry) <- rownames(med)
    ry
})

df_row <- do.call(cbind, df_row)    

ann_df <- colnames(df_row)
names(ann_df) <- ann_df
fig_1a <- TreeHeatmap(tree = treeR, tree_fig = fig_1, 
                      hm_data = df_row, rel_width = 0.4,
                      tree_hm_gap = 7.7, 
                      column_split = ann_df, 
                      #column_split_label = ann_df, 
                      column_split_gap = 0.1, 
                      # legend_title_hm = "Results", 
                      show_title = FALSE,
                      split_label_angle = 90,
                      split_label_size = 1.3, 
                      split_label_offset_x = 0.25,
                      split_label_offset_y = 1.7, 
                      split_label_fontface = "plain") +
    scale_fill_manual(values = c('FALSE' = "grey", 'TRUE' = "mediumvioletred"), 
                      guide = FALSE) +
     new_scale_fill() 
fig_1a +
    xlim(0, 15) +
    ylim(c(0, 10))

df_1a <- getData(tree_hm = fig_1a, type = "heatmap") %>%
    mutate(x1 = min(x)- width/2, y1 = max(y), x2 = max(x)+width/2) %>%
    select(x1, x2, y1) %>%
    distinct()
df_1b <- getData(tree_hm = fig_1a, type = "column_split") %>%
    rename(Method = column_split)

fig_1a <- fig_1a +
    geom_segment(data = df_1a, 
                 aes(x = x1, xend = x2, y = y1 + 1.5,
                     yend = y1 + 1.5), inherit.aes = FALSE) +
    geom_text(data = df_1a,
              aes(x = 0.5*(x1+x2), y = y1 + 1.8),
              inherit.aes = FALSE, label = "Results", size = 2)+
    geom_point(data = df_1b, 
               aes(x = x , y = y + 0.8, color = Method), 
               shape = 16, inherit.aes = FALSE) +
    scale_color_manual(values = c("treeclimbR" = "#E41A1C", 
                                 "StructFDR" = "#4DAF4A",
                                 "diffcyt" = "#377EB8",
                                 "HFDR" = "#984EA3",
                                 "lasso" = "#FF7F00", 
                                 "minP" = "#A65628"),
                      guide = guide_legend(order = 4)) +
    xlim(0, 15) +
    ylim(c(0, 10))
fig_1a +
    xlim(0, 14.5) +
    ylim(c(0, 10))

# --------------------------- annotate nodes -----------------------------------
lapply(loc_da, class)

df <- data.frame(node = unique(unlist(loc_da))) %>%
        mutate(treeclimbR = node %in% loc_da$treeclimbR,
               StructFDR = node %in% loc_da$StructFDR,
               diffcyt = node %in% loc_da$diffcyt,
               lasso = node %in% loc_da$lasso,
               minP = node %in% loc_da$minP,
               is_leaf = isLeaf(tree = treeR, node = node))
df[, 2:6] <- t(apply(df[, 2:6], 1, FUN = function(x){x/sum(x)}))



pies <- nodepie(df, 
                cols=2:6,  
                color = c("treeclimbR" = "#E41A1C", 
                          "StructFDR" = "#4DAF4A",
                          "diffcyt" = "#377EB8",
                          "lasso" = "#FF7F00", 
                          "minP" = "#A65628"))

hm_DA <- fig_1a +
    geom_inset(pies, width = 0.4, height = 0.4, vjust = .1) +
    #theme(legend.position = "none")
    theme(legend.box.margin = 
              margin(t = 3, b = 5, 
                     r = -10, l = -25),
          legend.text = element_text(size = 5),
          legend.position = "right",
          legend.spacing.y = unit(0.5, "mm"),
          legend.key.size= unit(2, "mm"),
          legend.title = element_text(size = 6)) +
    xlim(0, 14.5) +
    ylim(c(0, 10))

hm_DA
figPath <- "summary/hm_DA.eps"
ggsave(figPath, hm_DA, units = "in", width = 3, height = 3,
       dpi = 300)

