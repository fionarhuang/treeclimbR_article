# tse: the simulated data
# s : which repition was used
# anno_node: nodes annotated in the tree
# loc_da: nodes found by methods 
drawF <- function(tse, s = 5, anno_node, loc_da, scene) {
    
    # ---------------------- visualization level ----------------------
    # tree
    treeR <- rowTree(tse)
    
    # branches: different between groups
    fc <- metadata(tse)$FC
    leaf_diff <- list("increase" = names(fc[fc > 1]),
                      "decrease" = names(fc[fc < 1]))
    node_diff <- lapply(leaf_diff, FUN = function(x) {
        signalNode(tree = treeR, node = x) 
    })
    branch_diff <- findOS(tree = treeR, node = unlist(node_diff), 
                          only.leaf = FALSE, self.include = TRUE)
    names(branch_diff) <- rep(names(leaf_diff), 
                              unlist(lapply(node_diff, length)))
    
    
    # visualize level
    br <- metadata(tse)$br
    br <- findOS(tree = treeR, node = c(br$A, br$B),
                 only.leaf = TRUE, self.include = TRUE)
    leaf <- showNode(tree = treeR, only.leaf = TRUE)
    leaf_same <- setdiff(leaf, unlist(br))
    loc_same <- signalNode(tree = treeR, node = leaf_same)
    inner_A <- loc_same[!isLeaf(tree = treeR, node = loc_same)]
    vis_L <- loc_same
    comb_L <- vis_L[!isLeaf(tree = treeR, node =vis_L)]
    
    
    # --------------------------------- fig: tree -------------------------------
    # each edge connects a parent and a child;
    # edges that has the child node with signal are colored
    node_all <- showNode(tree = treeR, only.leaf = FALSE)
    node_same <- setdiff(node_all, unlist(branch_diff))
    df <- data.frame(node = c(unlist(branch_diff), node_same),
                     Truth  = c(rep(names(branch_diff), 
                                    unlist(lapply(branch_diff, length))),
                                rep("same", length(node_same))))
    
    treeS <- ggtree(treeR, branch.length = "none", 
                    size = 0.4) %<+% df +
        aes(color = Truth) +
        scale_color_manual(values = c("decrease" = "#D69C4E",
                                      "increase" = "#046C9A",
                                      "same" = "grey50"),
                           guide = FALSE)
    
    for (j in seq_along(comb_L)) {
        treeS <- treeS %>% collapse(node = comb_L[j])
    }
    
    treeS$data <- treeS$data %>%
        mutate(isTip = ifelse(node %in% vis_L, TRUE, isTip)) %>%
        mutate(labs = ifelse(isTip, NA, label))
    
    
    treeFig <- treeS + 
    geom_tiplab(align = TRUE, aes(label = labs), 
                show.legend = FALSE) +
    geom_point2(aes(subset = node %in% anno_node), shape =5,
                color = "red", size = 1.8) +
    geom_segment(data = . %>%
                     mutate(xend = max(x, na.rm = TRUE),
                            x = min(x, na.rm = TRUE),
                            y = max(y, na.rm = TRUE)) %>%
                     select(x, xend, y) %>%
                     distinct(),
                 aes(x = x, xend = xend, y = y + 5, yend = y + 5),
                 inherit.aes = FALSE, color = "black") +
    geom_text(data = . %>%
                  mutate(x = mean(range(x, na.rm = TRUE)),
                         y = max(y, na.rm = TRUE) + 6.8) %>%
                  select(x, y) %>%
                  distinct(),
              aes(x, y), inherit.aes = FALSE,
              label =  "Tree + detected nodes", size = 2) +
    geom_text(data = treeS$data %>%
                  mutate(x = 0.1*max(x, na.rm = TRUE),
                         y = max(y, na.rm = TRUE) - 8) %>%
                  select(x, y) %>%
                  distinct(),
              aes(x, y), inherit.aes = FALSE,
              label = scene, size = 3)
    
    treeFig
    # ------------------------------------------------------------------------------
    # get the order of tips & fake tips
    ro <- treeFig$data %>%
        filter(!is.na(y) & isTip) %>%
        arrange(desc(y)) %>%
        select(node) %>%
        unlist() 
    ro_label <- transNode(tree = treeR, node = ro, use.alias = FALSE)
    ro_alias <- transNode(tree = treeR, node = ro, use.alias = TRUE)
    
    
    # Heatmap
    med <- assays(tse)[[s]][ro_alias, ]
    med <- log(med+1)
    med <- data.frame(med)
    
    # column split & labels
    ann_c <- gsub(pattern = "_.*", "", colnames(med))
    names(ann_c) <- colnames(med)
    lab_c <- ifelse(ann_c == "C1", "Treated", "Control")
    names(lab_c) <- ann_c
    
    # tree & heatmap
    fig_1 <- TreeHeatmap(tree = treeR, tree_fig = treeFig, 
                         hm_data = med, rel_width = 0.8,
                         tree_hm_gap = 4, 
                         split_label_offset_y = 2,
                         column_split = ann_c, 
                         column_split_label = lab_c, 
                         column_split_gap = 1.5, 
                         split_label_fontface = "plain",
                         legend_title_hm = "Abundance (log)",
                         split_label_size = 2 ) +
        scale_fill_viridis_c(option = "D", 
                             breaks=scales::pretty_breaks(n=3),
                             limits = c(0, 2), oob = squish,
                             guide = guide_colorbar(order = 3), 
                             direction = 1) +
        # scale_fill_distiller(palette = "PiYG", 
        #                      breaks=scales::pretty_breaks(n=3),
        #                      direction = 1,
        #                      limits = c(0,3), oob = squish)+
        new_scale_fill()
    
    
    df_1 <- getData(tree_hm = fig_1, type = "heatmap")%>%
        mutate(x1 = min(x)- width/2, y1 = max(y), x2 = max(x)+width/2) %>%
        select(x1, x2, y1) %>%
        distinct()
    fig_1 <- fig_1 +
        geom_segment(data = df_1, 
                     aes(x = x1, xend = x2, y = y1 + 5,
                         yend = y1 + 5), inherit.aes = FALSE) +
        geom_text(data = df_1,
                  aes(x = 0.5*(x1+x2), y = y1 + 6.8),
                  inherit.aes = FALSE, label = "Observations", size = 2)
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
                          tree_hm_gap = 45, 
                          column_split = ann_df, 
                          #column_split_label = ann_df, 
                          column_split_gap = 0.2, 
                          legend_title_hm = "Results", 
                          show_title = FALSE,
                          split_label_angle = 90,
                          split_label_size = 1.3, 
                          split_label_offset_x = 0.25,
                          split_label_offset_y = 1.7, 
                          split_label_fontface = "plain") +
        scale_fill_manual(values = c('FALSE' = "grey", 'TRUE' = "mediumvioletred"),
                          labels = c('FALSE' = "same", 'TRUE' = "DA" )) +
        ylim(c(0, 95)) +
        new_scale_fill() +
        new_scale_color()
    fig_1a 
    
    
    # tree & heatmap & column annotations
    df_1a <- getData(tree_hm = fig_1a, type = "heatmap") %>%
        mutate(x1 = min(x)- width/2, y1 = max(y), x2 = max(x)+width/2) %>%
        select(x1, x2, y1) %>%
        distinct()
    df_1b <- getData(tree_hm = fig_1a, type = "column_split") %>%
        dplyr::rename(Method = column_split) %>%
        mutate(Method = factor(Method, 
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
        arrange(Method)
    
    fig_1b <- fig_1a +
        geom_segment(data = df_1a, 
                     aes(x = x1, xend = x2, y = y1 + 5,
                         yend = y1 + 5), inherit.aes = FALSE) +
        geom_text(data = df_1a,
                  aes(x = 0.5*(x1+x2), y = y1 + 6.8),
                  inherit.aes = FALSE, label = "Results", size = 2)+
        geom_point(data = df_1b, 
                   aes(x = x , y = y + 2.5, color = Method, shape = Method), 
                   inherit.aes = FALSE) +
        scale_color_manual(values = c("treeclimbR" = "#E41A1C", 
                                      "BH" = "#377EB8",
                                      "StructFDR" = "#4DAF4A", 
                                      "HFDR" = "#984EA3",
                                      "lasso" = "#FF7F00", 
                                      "minP" = "#A65628",
                                      "miLineage1" = "#999999",
                                      "miLineage2" = "#666666",
                                      "LEfSe" = "#E7298A"),
                           guide = guide_legend(order = 4)) +
        scale_shape_manual(values = c("treeclimbR" = 5, "BH" = 16,
                                      "StructFDR" = 16, "HFDR" = 16,
                                      "lasso" = 16, "minP" = 16,
                                      "miLineage1" = 16, "miLineage2" = 16,
                                      "LEfSe" = 16),
                           guide = guide_legend(order = 4)) +
        theme(legend.box.margin = 
                  margin(t = 3, b = 5, 
                         r = -10, l = -25),
              legend.text = element_text(size = 5),
              legend.position = "right",
              legend.spacing.y = unit(0.5, "mm"),
              legend.key.size= unit(2.5, "mm"),
              legend.title = element_text(size = 6))
    fig_1b
}

