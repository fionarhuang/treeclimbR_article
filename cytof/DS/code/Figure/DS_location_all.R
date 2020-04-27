
library(ggtree)
library(TreeHeatmap)
library(TreeSummarizedExperiment)
library(treeclimbR)
library(ggnewscale)
library(scales)

# ---------------------- load data ----------------------
seed <- 2020
nclus <- 400
md <- "DS/output/RData"
sel <- sprintf("less_distinct_50pc_clus%d_seed%d", nclus, seed)
ana <- c("treeclimbR", "diffcyt", "LassoGLM", "minP")

dataFile <- file.path(md, "DataPrep", paste0(sel, ".RData"))
analysisFile <- paste0(md, "/Analysis/", sel, "/", ana, ".RData")
lapply(c(dataFile, analysisFile), load, .GlobalEnv)


# ---------------------- visualization level ----------------------
# tree
treeR <- rowTree(d_medians_all)

# true B cell (threshold > 0.5)
leaf_B <- df_truth$cluster_id[df_truth$prop_B > 0.1]
node_B <- signalNode(tree = treeR, 
                     node = df_truth$cluster_id[df_truth$prop_B > 0.5])
branch_B <- findOS(tree = treeR, node = node_B, 
                   only.leaf = TRUE, self.include = TRUE)
names(branch_B) <- rep("B cell", length(branch_B))

# not B cell
leaf_A <- setdiff(df_truth$cluster_id, leaf_B)
node_A <- signalNode(tree = treeR, node = leaf_A)
inner_A <- node_A[!isLeaf(tree = treeR, node = node_A)]


# visualize level
vis_L <- c(node_A)
comb_L <- vis_L[!isLeaf(tree = treeR, node =vis_L)]

# --------------------------------- result -------------------------------------
pS6_treecb <- res$`0.05`$output %>%
    filter(marker_id %in% "pS6")
loc_treeclimbR_ds <- pS6_treecb$node[pS6_treecb$signal.node %in% TRUE]

loc_diffcyt_ds <- res_leaf %>%
    data.frame() %>%
    mutate(node = transNode(tree = treeR, node = as.character(cluster_id))) %>%
    filter(p_adj <= 0.05 & marker_id == "pS6") %>%
    select(node) %>% unlist()
loc_minP_ds <- loc.minP_ab$`0.05`$pS6
loc_lasso_ds <- loc.Lasso_ab$pS6

# --------------------------------- tree ---------------------------------------
treeS <- viewBranch(tree = treeR, group_leaf = branch_B, 
                    group_color = c("0" = "grey50", "B cell" = "mediumvioletred"),
                    edge_size = 0.6) +
    scale_color_manual(name = "B_cell (> 50%)", 
                       values = c("0" = "grey50", "B cell" = "mediumvioletred"), 
                       labels = c("NO", "Yes")) +
    guides(color = FALSE)

for (i in seq_along(comb_L)) {
    treeS <- treeS %>% collapse(node = comb_L[i])
}

treeS$data <- treeS$data %>%
    mutate(isTip = ifelse(node %in% vis_L, TRUE, isTip)) %>%
    mutate(labs = ifelse(isTip, NA, label))
treeFig <- treeS + 
    geom_tiplab(align = TRUE, aes(label = labs), show.legend = FALSE) +
    #    geom_point2(aes(subset = node %in% vis_L)) +
    geom_point2(aes(subset = node %in% loc_treeclimbR_ds), shape = 5,
                color = "red", size = 3, stroke = 0.8) +
    geom_point2(aes(subset = node %in% loc_lasso_ds), shape = 16,
                color = "#FF7F00", size = 1.8, stroke = 1) +
    geom_point2(aes(subset = node %in% loc_minP_ds), shape = 16,
                color = "#A65628", size = 2.8) +
    geom_point2(aes(subset = node %in% loc_diffcyt_ds), shape = 16,
                color = "#377EB8", size = 1.2) +
    
    ylim(1, 55) +
    xlim(0, 25)

treeFig



# get the order of tips
ro <- treeFig$data %>%
    dplyr::filter(!is.na(y) & isTip) %>%
    arrange(desc(y)) %>%
    select(label) %>%
    unlist()

ro_leaf <- findOS(tree = treeR, node = ro, 
                  only.leaf = TRUE, self.include = TRUE)
ro_leaf <- lapply(ro_leaf, FUN = function(x){
    transNode(tree = treeR, node = x, use.alias = FALSE)
})
ro_df <- data.frame(group = rep(names(ro_leaf), unlist(lapply(ro_leaf, length))),
                    cluster_id = unlist(ro_leaf))
ro_B <- df_truth %>% left_join(ro_df, by = "cluster_id") %>%
    group_by(group) %>%
    summarise(n_cells = sum(n_cells),
              n_B = sum(prop_B*n_cells),
              prop_B = n_B/n_cells)
pr_B <- ro_B$prop_B
names(pr_B) <- ro_B$group

B_df <- treeFig$data %>%
    dplyr::filter(!is.na(y) & isTip) %>%
    arrange(desc(y)) %>%
    select(node, y) %>%
    mutate("B_pct" = pr_B[ro] *100) 


treeFig <- treeFig + 
    geom_tile(data = B_df, 
              aes(x = 19, y = y, width = 1, fill = B_pct), 
              inherit.aes = FALSE) +
    annotate("text", x = 19, y = 55, label = "A", 
              color = "black", size = 3) +
    scale_fill_gradient(name = "A: B cell (%)",
                        low = "white", high = "mediumvioletred",
                        breaks=c(0, 25, 50, 75, 100),
                        limits = c(0, 100)) +
    # scale_fill_distiller(name = "A: B cell (%)", palette = "PuOr",
    #                      breaks=scales::pretty_breaks(n=5),
    #                      direction = 1) +
    guides(fill = guide_legend(order = 1,
                               override.aes = list(size = 2.5))) +
    new_scale_fill() 


# ----------------- Heatmap -----------------
med <- assays(d_medians_all)[["pS6"]][rowLinks(d_medians_all)$nodeLab %in% ro, ]
med[is.na(med)] <- 0

#med <- sqrt(med)


# column split & labels
ann_c <- gsub(pattern = ".*_", "", colnames(med))
names(ann_c) <- colnames(med)
lab_c <- ifelse(ann_c == "base", "Control", "Stimulated")
names(lab_c) <- ann_c

# tree & heatmap
hm_DS <- TreeHeatmap(tree = treeR, tree_fig = treeFig, 
                    hm_data = med, rel_width = 0.8,
                    tree_hm_gap = 3, 
                    column_split = ann_c, 
                    column_split_label = lab_c, 
                    column_split_gap = 0.5, 
                    legend_title_hm = "Expression",
                    split_label_size = 3, 
                    split_label_offset_y = 2,
                    split_label_fontface = "plain") +
    # scale_fill_distiller(palette = "PRGn",
    #                      breaks=scales::pretty_breaks(n=3),
    #                      direction = 1,
    #                      limits = c(-2,2), oob = squish) +
    scale_fill_gradient(low = "navy", high = "yellow",
                        limits = c(0, 3),
                        oob = squish) +
    guides(fill = guide_legend(order = 2,
                               override.aes = list(size = 2.5))) +
    new_scale_fill() +
    xlim(0, 35) +
    ylim(0, 55)

hm_DS









