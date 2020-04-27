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
    library(cowplot)
})


# load data
scene <- "less_distinct_50pc"
load(sprintf("DS/output/RData/DataPrep/%s_clus400_seed2020.RData", scene))
# tree
treeR <- rowTree(d_medians_all)


# load results
path0 <- sprintf("DS/output/RData/Analysis/%s_clus400_seed2020/", scene)
method <- c("diffcyt", "LassoGLM", "minP", "treeclimbR")
path1 <- paste0(path0, method, ".RData")
lapply(path1, load, .GlobalEnv)

find.treeclimbR <- loc.treeclimbR_ab$`0.05`
find.diffcyt <- loc.diffcyt_ab$`0.05`
find.lasso <- loc.Lasso_ab
find.minP <- loc.minP_ab$`0.05`

# main difference: "pS6", "pPlcg2", "pErk"
df <- lapply(c("pS6", "pPlcg2", "pErk", "pNFkB"), FUN = function(x) {
    tcr <- find.treeclimbR[[x]]
    dft <- find.diffcyt[[x]]
    las <- find.lasso[[x]]
    mip <- find.minP[[x]]
    
    xx <- data.frame(node = unique(c(tcr, dft, las, mip))) %>%
        mutate(treeclimbR = node %in% tcr,
               diffcyt = node %in% dft,
               lasso = node %in% las,
               minP = node %in% mip,
               is_leaf = isLeaf(tree = treeR, node = node))
    xx[, 2:5] <- t(apply(xx[, 2:5], 1, FUN = function(x){x/sum(x)}))
    return(xx)
})
names(df) <- c("pS6", "pPlcg2", "pErk", "pNFkB")



# ---------------------------- visualization level -----------------------------


# true B cell (threshold > 0.5)
ts <- 0.5
leaf_B <- df_truth$cluster_id[df_truth$prop_B > ts]
node_B <- signalNode(tree = treeR, node = leaf_B)
branch_B <- findOS(tree = treeR, node = node_B, 
                   only.leaf = TRUE, self.include = TRUE)
names(branch_B) <- rep("B cell", length(branch_B))


# not B cell
leaf_A <- setdiff(df_truth$cluster_id, leaf_B)
node_A <- signalNode(tree = treeR, node = leaf_A)
inner_A <- node_A[!isLeaf(tree = treeR, node = node_A)]


# visualize level in non-B branch
vis_A <- c(node_A)
comb_A <- vis_A[!isLeaf(tree = treeR, node =vis_A)]

# visualize level
vis_L <- c(node_A, 
           transNode(treeR, as.character(leaf_B)))
# --------------------------------- tree ---------------------------------------
treeS <- ggtree(treeR, size = 0.4)
br_B <- findOS(tree = treeR, node = node_B, 
               only.leaf = FALSE, self.include = TRUE)
df_B <- treeS$data %>%
    select(node) %>%
    mutate(Type = ifelse(node %in% unlist(br_B), "B_cell", "others"))
treeS <- treeS %<+% df_B +
    aes(color = Type) + 
    scale_color_manual(values = c("others" = "grey50", "B_cell" = "mediumvioletred"), 
                       labels = c("NO", "Yes"), guide = FALSE) +
    new_scale_color()

for (i in seq_along(comb_A)) {
    treeS <- treeS %>% collapse(node = comb_A[i])
}

treeS$data <- treeS$data %>%
    mutate(isTip = ifelse(node %in% vis_A, TRUE, isTip)) %>%
    mutate(labs = ifelse(isTip, NA, label))  

treeFig <- treeS + 
    geom_tiplab(align = TRUE, aes(label = labs), show.legend = FALSE) +
    ylim(1, 35) +
    xlim(0, 10)
treeFig 
# treeFig +
#     geom_point2(aes(subset = (node %in% 479)))

#---------------------- Truth ------------------------------------------------
# calculate the number of B cells for each node in visualized level
leaf_L <- findOS(tree = treeR, node = vis_L, 
                 only.leaf = TRUE, self.include = TRUE)
B_truth <- lapply(seq_along(leaf_L), FUN = function(x) {
    xx <- transNode(tree = treeR, node = leaf_L[[x]], use.alias = FALSE)
    dx <- df_truth %>%
        filter(cluster_id %in% xx) %>%
        mutate(node = vis_L[x],
               B_n = n_cells*prop_B,
               nonB_n = n_cells*(1-prop_B)) %>%
        summarise(node = unique(node),
                  B_n = sum(B_n),
                  nonB_n = sum(nonB_n))
    
})
B_truth <- do.call(rbind, B_truth)
c <- 0.3 * max(treeFig$data$x, na.rm = TRUE)
tree_bar_gap <- 0.5
bar_point_gap <- 0.5
B_df <- B_truth %>%
    left_join(treeFig$data) %>%
    select(node, x, y, B_n, nonB_n) %>%
    mutate(B_cell = round(B_n/(B_n + nonB_n), digits = 6),
           others = B_cell) %>%
    gather(Type, prop, B_cell:others) %>%
    mutate(x_start = ifelse(Type == "B_cell", 0, prop*c)+ max(x)+tree_bar_gap,
           x_end = ifelse(Type == "B_cell", prop*c, 1*c)+ max(x)+tree_bar_gap,
           Total_cell = B_n+nonB_n)
fig_0 <- treeFig +
    xlim(c(0, 10))+
    ylim(c(0, 40))+
    geom_segment(data = treeFig$data %>%
                     mutate(xstart = min(x, na.rm = TRUE),
                            xend = max(x, na.rm = TRUE),
                            y = max(y, na.rm = TRUE) + 3) %>%
                     select(xstart, xend, y) %>%
                     distinct(),
                 aes(x = xstart, xend = xend,
                     y = y , yend = y) ,
                 inherit.aes = FALSE,
                 color = "black") +
    annotate("text", x = 0.5*(min(treeFig$data$x, na.rm = TRUE) + 
                                  max(treeFig$data$x, na.rm = TRUE)), 
             y = max(treeFig$data$y, na.rm = TRUE) + 3.6,
             label = "Tree + detected nodes", size = 2)+
    # annotate("text", x = 8, y = 34.5,
    #          label = "BCR-XL-sim \n medium, 400", size = 2) +
    geom_rect(data = B_df, 
              aes(xmin = x_start, xmax = x_end, 
                  ymin = y - 0.45, ymax = y + 0.45,
                  fill = Type),
              inherit.aes = FALSE) +
    geom_point(data = B_df %>% filter(Type == "others"),
               aes(x = x_end + bar_point_gap, y = y, size = Total_cell),
               inherit.aes = FALSE) +
    geom_segment(data = B_df %>%
                     mutate(x_start = min(x_start), 
                            x_end = max(x_end) + bar_point_gap,
                            y_end = max(y) + 3) %>%
                     select(x_start, x_end, y_end) %>%
                     distinct(), 
                 aes(x = x_start, xend = x_end, y = y_end,
                     yend = y_end), inherit.aes = FALSE) +
    annotate("text", x = 0.5*(min(B_df$x_start) + max(B_df$x_end+ bar_point_gap)),  
             y = max(B_df$y) + 3.6,
             label = "Truth", size = 2)+
    scale_fill_manual(values = c("B_cell" = "mediumvioletred", "others" = "grey"),
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

med <- lapply(c("pS6", "pPlcg2", "pErk", "pNFkB"), FUN = function(x){
    xx <- assays(d_medians_all)[[x]][rowLinks(d_medians_all)$nodeLab %in% ro, ]
    xx[is.na(xx)] <- 0
    #t(scale(t(xx)))
    return(xx)
})
names(med) <- c("pS6", "pPlcg2", "pErk", "pNFkB")

# column split & labels
ann_c <- gsub(pattern = ".*_", "", colnames(d_medians_all))
names(ann_c) <- colnames(d_medians_all)
lab_c <- ifelse(ann_c == "base", "Control", "Stimulated")
names(lab_c) <- ann_c

# =============================== pS6 ==========================================
# tree & heatmap
fig_1 <- TreeHeatmap(tree = treeR, tree_fig = fig_0, 
                     hm_data = med[[1]], rel_width = 1,
                     tree_hm_gap = 2.5, 
                     column_split = ann_c, 
                     column_split_label = lab_c, 
                     column_split_gap = 0.1, 
                     legend_title_hm = "Expression",
                     split_label_size = 2, 
                     split_label_offset_y = 1.5,
                     split_label_fontface = "plain") +
    scale_fill_gradient(low = "navy", high = "yellow",
                        limits = c(0, 3),
                        oob = squish) +
    guides(fill = guide_legend(order = 3,
                               override.aes = list(size = 2.5))) +
    new_scale_fill() +
    xlim(0, 15) +
    ylim(0, 40)
fig_1
df_1 <- getData(tree_hm = fig_1, type = "heatmap")%>%
    mutate(x1 = min(x)- width/2, y1 = max(y), x2 = max(x)+width/2) %>%
    select(x1, x2, y1) %>%
    distinct()
fig_1 <- fig_1 +
    geom_segment(data = df_1, 
                 aes(x = x1, xend = x2, y = y1 + 3,
                     yend = y1 + 3), inherit.aes = FALSE) +
    geom_text(data = df_1,
              aes(x = 0.5*(x1+x2), y = y1 + 3.6),
              inherit.aes = FALSE, label = "pS6 expression", size = 2)
fig_1

# --------------------------- annotate heatmap rows ---------------------------
ab <- "pS6"
loc <- c(find.treeclimbR[ab], find.diffcyt[ab], find.lasso[ab], find.minP[ab])
names(loc) <- c("treeclimbR", "diffcyt", "lasso", "minP")
df_row <- lapply(loc, FUN = function(y){
    yy <- unlist(findOS(tree = treeR, node = y, 
                        only.leaf = FALSE, 
                        self.include = TRUE))
    ty <- transNode(tree = treeR, node = yy, use.alias = TRUE)
    ry <- rownames(med[[ab]]) %in% ty
    names(ry) <- rownames(med[[ab]])
    ry
})

df_row <- do.call(cbind, df_row)    

ann_df <- colnames(df_row)
names(ann_df) <- ann_df
fig_1a <- TreeHeatmap(tree = treeR, tree_fig = fig_1, 
                      hm_data = df_row, rel_width = 0.2,
                      tree_hm_gap = 6.3, 
                      column_split = ann_df, 
                      # column_split_label = ann_df, 
                      column_split_gap = 0.1, 
                      # legend_title_hm = "Results", 
                      show_title = FALSE,
                      split_label_angle = 45,
                      split_label_size = 1.5, 
                      split_label_offset_x = 0.5,
                      split_label_offset_y = 1.5, 
                      split_label_fontface = "plain") +
    scale_fill_manual(values = c('FALSE' = "grey", 'TRUE' = "mediumvioletred"), 
                      guide = FALSE) +
    new_scale_fill() 
fig_1a
df_1a <- getData(tree_hm = fig_1a, type = "heatmap") %>%
    mutate(x1 = min(x)- width/2, y1 = max(y), x2 = max(x)+width/2) %>%
    select(x1, x2, y1) %>%
    distinct()
df_1b <- getData(tree_hm = fig_1a, type = "column_split") %>%
    rename(Method = column_split)

fig_1a <- fig_1a +
    geom_segment(data = df_1a, 
                 aes(x = x1, xend = x2, y = y1 + 3,
                     yend = y1 + 3), inherit.aes = FALSE) +
    geom_text(data = df_1a,
              aes(x = 0.5*(x1+x2), y = y1 + 3.6),
              inherit.aes = FALSE, label = "Results", size = 2) +
    geom_point(data = df_1b, 
               aes(x = x , y = y + 1.5, color = Method), 
               shape = 16, inherit.aes = FALSE) +
    scale_color_manual(values = c("treeclimbR" = "#E41A1C", "diffcyt" = "#377EB8",
                                  "lasso" = "#FF7F00", "minP" = "#A65628"),
                       guide = guide_legend(order = 4)) +
    xlim(0, 11) +
    ylim(0, 38)
fig_1a +
    xlim(0, 11) +
    ylim(0, 38)
# --------------------------- annotate nodes -----------------------------------

pie_leaf <- nodepie(df[[1]], cols=2:5,  
                    color = c("treeclimbR" = "#E41A1C", "diffcyt" = "#377EB8",
                              "lasso" = "#FF7F00", "minP" = "#A65628"))
hm_DS <- fig_1a +
    #geom_inset(pie_in, width = 0.2, height =  0.2) +
    geom_inset(pie_leaf, width = 0.3, height = 0.3, vjust = 0.1, hjust = 0.05) +
    theme(legend.box.margin = 
              margin(t = 3, b = 5, 
                     r = -10, l = -25),
          legend.text = element_text(size = 5),
          legend.position = "right",
          legend.spacing.y = unit(0.5, "mm"),
          legend.key.size= unit(2, "mm"),
          legend.title = element_text(size = 6))

# ============================= 4 other markers ========================================
#  add text in tree figure  
fig_0 <- fig_0 + 
    annotate("text", x = 0.5 * max(treeFig$data$x, na.rm = TRUE), y = 34.5,
             label = "BCR-XL-sim \n medium, 400", size = 2) 



ab_limit <- list("pPlcg2" = c(1, 2), "pErk" = c(0.5, 2), 
                 "pNFkB" =  c(1.5, 3))
pie_width <- pie_height <- 0.35
fig_list <- vector("list", length(ab_limit))

for (i in seq_along(ab_limit)) {
    ab <- names(ab_limit)[i]
    fig_m <- TreeHeatmap(tree = treeR, tree_fig = fig_0,
                         hm_data = med[[ab]], rel_width = 1,
                         tree_hm_gap = 2.5,
                         column_split = ann_c,
                         column_split_label = lab_c,
                         column_split_gap = 0.1,
                         legend_title_hm = "Expression",
                         split_label_size = 2,
                         split_label_offset_y = 1.5,
                         split_label_fontface = "plain") +
        scale_fill_gradient(low = "navy", high = "yellow",
                            limits = ab_limit[[ab]],
                            oob = squish, 
                            guide = guide_colourbar(order = 3,
                                                    override.aes = list(size = 2.5))) +
        new_scale_fill() +
        xlim(0, 15) +
        ylim(0, 40)
    fig_m
    df_m <- getData(tree_hm = fig_m, type = "heatmap")%>%
        mutate(x1 = min(x)- width/2, y1 = max(y), x2 = max(x)+width/2) %>%
        select(x1, x2, y1) %>%
        distinct()
    fig_m <- fig_m +
        geom_segment(data = df_m,
                     aes(x = x1, xend = x2, y = y1 + 3,
                         yend = y1 + 3), inherit.aes = FALSE) +
        geom_text(data = df_m,
                  aes(x = 0.5*(x1+x2), y = y1 + 3.6),
                  inherit.aes = FALSE, label = paste(ab, "expression", sep = " "), size = 2)
    # --------------------------- annotate heatmap rows ---------------------------
    
    loc <- c(find.treeclimbR[ab], find.diffcyt[ab], find.lasso[ab], find.minP[ab])
    names(loc) <- c("treeclimbR", "diffcyt", "lasso", "minP")
    df_row <- lapply(loc, FUN = function(y){
        if (length(y)) {
            yy <- unlist(findOS(tree = treeR, node = y,
                                only.leaf = FALSE,
                                self.include = TRUE))
            ty <- transNode(tree = treeR, node = yy, use.alias = TRUE)
            
        } else {
            ty <- NULL
        }
        ry <- rownames(med[[ab]]) %in% ty
        names(ry) <- rownames(med[[ab]])
        ry
    })
    
    df_row <- do.call(cbind, df_row)
    
    ann_df <- colnames(df_row)
    names(ann_df) <- ann_df
    fig_ma <- TreeHeatmap(tree = treeR, tree_fig = fig_m,
                          hm_data = df_row, rel_width = 0.2,
                          tree_hm_gap = 6.3,
                          column_split = ann_df,
                          column_split_label = NULL,
                          column_split_gap = 0.1,
                          # legend_title_hm = "Results",
                          show_title = FALSE,
                          split_label_angle = 45,
                          split_label_size = 1,
                          split_label_offset_x = 0.5,
                          split_label_offset_y = 1.5,
                          split_label_fontface = "plain") +
        scale_fill_manual(values = c('FALSE' = "grey", 'TRUE' = "mediumvioletred"),
                          guide = FALSE) +
        #guides(fill = guide_legend(order = 2,
        #                          override.aes = list(size = 2.5))) +
        new_scale_fill() 
    
    df_ma <- getData(tree_hm = fig_ma, type = "heatmap") %>%
        mutate(x1 = min(x)- width/2, y1 = max(y), x2 = max(x)+width/2) %>%
        select(x1, x2, y1) %>%
        distinct()
    df_mb <- getData(tree_hm = fig_ma, type = "column_split") %>%
        rename(Method = column_split)
    
    fig_ma <- fig_ma +
        geom_segment(data = df_ma,
                     aes(x = x1, xend = x2, y = y1 + 3,
                         yend = y1 + 3), inherit.aes = FALSE) +
        geom_text(data = df_ma,
                  aes(x = 0.5*(x1+x2), y = y1 + 3.6),
                  inherit.aes = FALSE, label = "Results", size = 2) +
        geom_point(data = df_mb, 
                   aes(x = x , y = y + 1.5, color = Method), 
                   shape = 16, inherit.aes = FALSE) +
        scale_color_manual(values = c("treeclimbR" = "#E41A1C", "diffcyt" = "#377EB8",
                                      "lasso" = "#FF7F00", "minP" = "#A65628"),
                           guide = guide_legend(order = 4))
    
    
    
    # --------------------------- annotate nodes -----------------------------------
    pie_leaf <- nodepie(df[[ab]], cols=2:5,
                        color = c("treeclimbR" = "#E41A1C", "diffcyt" = "#377EB8",
                                  "lasso" = "#FF7F00", "minP" = "#A65628"))
    Fig_final <- fig_ma +
        geom_inset(pie_leaf, width = pie_width, height =  pie_height, 
                   vjust = 0.1, hjust = 0.05) +
        theme(legend.box.margin = 
                  margin(t = 3, b = 5, 
                         r = -10, l = -25),
              legend.text = element_text(size = 5),
              legend.position = "right",
              legend.spacing.y = unit(0.5, "mm"),
              legend.key.size= unit(2, "mm"),
              legend.title = element_text(size = 6))
    fig_list[[i]] <-  Fig_final +
        labs(title = ab) + 
        theme(
            plot.title = element_text(hjust= 0.5, 
                                      vjust = -5,
                                      face = "plain",
                                      size = 10),
            legend.box.margin = 
                margin(t = 3, b = 5, 
                       r = 5, l = -25),
            legend.text = element_text(size = 5),
            legend.position = "right",
            legend.spacing.y = unit(0.5, "mm"),
            legend.key.size= unit(2, "mm"),
            legend.title = element_text(size = 6)) +
        xlim(0, 11) +
        ylim(0, 38)
} 

fig_more <- plot_grid(plotlist = fig_list, nrow = 2, 
                      rel_heights = c(0.5, 0.5))



figPath <- "summary/cytof_2020_more.eps"
ggsave(figPath, fig_more, units = "in", width = 8, height = 8.5,
       dpi = 300)

