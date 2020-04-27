suppressPackageStartupMessages({
    library(TreeSummarizedExperiment)
    library(treeclimbR)
    library(diffcyt)
    library(dplyr)
    library(ape)
    library(tibble)
    library(tidyr)
    library(ggplot2)
})

seedN <- 2020
reso <- c(3, 5, 7, 10, 15, 20, 30, 40, 50)
res_auc <- list(diffcyt = rep(NA, length(reso)),
                treeclimbR = rep(NA, length(reso)))

# load data (directory: cytof/)
d_SE <- readRDS("data/Weber_BCR_XL_sim_less_distinct_less_50pc_SE.rds")

# experiment information
experiment_info <- metadata(d_SE)$experiment_info 

# marker information
marker_info <- as.data.frame(colData(d_SE))
rownames(marker_info) <- NULL  
marker_info

# keep original data
d_SE_original <- d_SE

# check levels of 'sample_id' column in row data are in correct order
rowData(d_SE)
levels(rowData(d_SE)$sample_id)

# split data (note: using 'split.data.frame' to split into list of matrices; see '?split')
d_input <- split.data.frame(assay(d_SE), rowData(d_SE)$sample_id)
d_se <- prepareData(d_input, experiment_info, marker_info)
d_se <- transformData(d_se, cofactor = 5)

# check marker classes
colnames(d_se)[colData(d_se)$marker_class == "type"]
colnames(d_se)[colData(d_se)$marker_class == "state"]

for (i in seq_along(reso)) {
    message(i)
    nx <- ny <- reso[i]
    
    # run clustering (using FlowSOM with 20x20 grid, i.e. 400 clusters)
    d_se <- generateClusters(d_se, xdim = nx, ydim = ny, seed_clustering = seedN)
    rowData(d_se) <- rowData(d_se) %>%
        data.frame() %>%
        mutate(spikein = rowData(d_SE_original)$spikein) %>%
        mutate(B_cell = rowData(d_SE_original)$B_cell) %>%
        mutate(population_id = rowData(d_SE_original)$population_id)
    
    # # -------------------------------------------------------------
    # create the tree
    # -------------------------------------------------------------
    cytof_tree <- buildTree(d_se = d_se)
    
    
    # mapping: cells to nodes
    cell_info <- rowData(d_se) %>%
        data.frame() %>%
        mutate_if(is.factor, as.character) %>%
        rownames_to_column() %>%
        mutate(spikein = rowData(d_SE_original)$spikein,
               B_cell = rowData(d_SE_original)$B_cell,
               population_id = rowData(d_SE_original)$population_id,
               cell_id = paste0("cell_", rowname),
               node = transNode(tree = cytof_tree, node = cluster_id))  %>%
        select(cluster_id, node, cell_id, B_cell) 
    head(cell_info)
    
    
    
    # -------------------------------------------------------------
    # calculate medians, counts for each cluster
    # -------------------------------------------------------------
    # medians
    d_medians_all <- calcTreeMedians(d_se = d_se, tree = cytof_tree)
    
    # counts
    d_counts_all <- calcTreeCounts(d_se = d_se, tree = cytof_tree)
    
    
    
    # ==============================================================================
    ## Analysis
    # ==============================================================================
    # set up design matrix
    # note: include fixed effects for 'patient_id'
    design <- createDesignMatrix(experiment_info, cols_design = 1:2)
    design
    
    # set up contrast matrix
    contrast_vec <- c(0, 1, 0, 0, 0, 0, 0, 0, 0)
    contrast <- createContrast(contrast_vec)
    contrast
    
    ################################ diffcyt ######################################
    # data on the leaf level
    d_counts_leaf <- d_counts_all[rowLinks(d_counts_all)$isLeaf, ]
    rowData(d_counts_leaf) <- droplevels(rowData(d_counts_leaf))
    d_medians_leaf <- d_medians_all[rowLinks(d_medians_all)$isLeaf, ]
    rowData(d_medians_leaf) <- droplevels(rowData(d_medians_leaf))
    
    # run diffcyt
    DS_leaf <- testDS_limma(d_counts = d_counts_leaf, 
                            d_medians = d_medians_leaf, 
                            design, contrast)
    
    
    
    # top DA clusters
    res_leaf <- diffcyt::topTable(DS_leaf, format_vals = TRUE, 
                                  top_n = nrow(DS_leaf))
    
    diffcyt_pS6 <- res_leaf %>%
        data.frame() %>%
        mutate_if(is.factor, as.character) %>%
        filter(marker_id == "pS6" & !is.na(p_val)) %>%
        mutate(node = transNode(tree = cytof_tree, node = cluster_id)) %>%
        select(node, p_val) 
    
    out_diffcyt <- cell_info %>%
        inner_join(diffcyt_pS6)
    
    
    ################################ treeclimbR ######################################
    
    # run limma
    res_tree <- testDS_limma(d_counts_all, d_medians_all, 
                             min_cells = 10, min_samples = 8,
                             design, contrast)
    
    res_DS <- rowData(res_tree) %>%
        data.frame() %>%
        mutate(node = transNode(tree = cytof_tree, 
                                node = as.character(cluster_id))) %>%
        mutate(cluster_id = transNode(tree = cytof_tree, 
                                      node = node)) 
    
    
    # generate candidates
    res_list <- split(res_DS, f = res_DS$marker_id)
    cand <- lapply(res_list, FUN = function(x) {
        xx <- getCand(tree = cytof_tree, 
                      score_data = x, 
                      node_column = "node", p_column = "p_val", 
                      sign_column = "logFC", message = FALSE)
        xx$candidate_list})
    
    # the candidate is selected at FDR 0.05
    best <-  evalCand(tree = cytof_tree, type = "multiple", levels = cand,
                      score_data = res_list, node_column = "node", 
                      p_column = "p_val", sign_column = "logFC",
                      feature_column = "marker_id", limit_rej = 0.05)
    
    treeclimbR_pS6 <- topNodes(object = best, n = Inf, p_value = 1) 
    treeclimbR_pS6 <- best$output %>%
        filter(marker_id == "pS6") %>%
        select(node, p_val) %>%
        mutate(node = findOS(tree = cytof_tree, node = node, 
                             only.leaf = TRUE, self.include = TRUE)) %>%
        unnest(node)
    
    out_treeclimbR <- cell_info %>%
        inner_join(treeclimbR_pS6)
    
    head(out_treeclimbR)
    
    
    ## ============================================================================
    ## calculate AUC
    ## ============================================================================
    auc <- function(p_value, binary_truth, n_rep = 10000){
        p_pos <- 1 - p_value[binary_truth]
        p_neg <- 1 - p_value[!binary_truth]
        mean(sample(p_pos, n_rep, replace=T) > sample(p_neg, n_rep, replace=T))
    }
    res_auc$treeclimbR[i] <-  auc(p_value = out_treeclimbR$p_val, 
                                  binary_truth = out_treeclimbR$B_cell)
    res_auc$diffcyt[i] <- auc(p_value = out_diffcyt$p_val, 
                              binary_truth = out_diffcyt$B_cell)
}

save(res_auc, reso, file = "DS_resolution/res_auc.RData")

df_auc <- data.frame(n_cluster = reso^2, 
                     diffcyt = res_auc$diffcyt, 
                     treeclimbR = res_auc$treeclimbR) %>%
    gather(Method, AUC, -n_cluster)
ggplot(df_auc, aes(x = n_cluster, y = AUC, color = method)) +
    geom_point() +
    geom_line() +
    scale_color_manual(values = c("treeclimbR" = "#E41A1C", "diffcyt" = "#377EB8")) +
    scale_x_sqrt(labels = reso^2, breaks = reso^2) +
    ylim(c(0, 1)) + 
    labs(x = "The number of clusters") +
    theme_bw(base_size = 8) + theme(
        aspect.ratio = 1,
        panel.grid = element_blank(),
        panel.spacing = unit(0, "lines"),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size= unit(2.5, "mm"),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 8),
        legend.position=c(0.5, 0.5),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin=margin(-5,-5,-10,-5),
        strip.background = element_rect(colour = "black", fill = "gray90"),
        strip.text.x = element_text(color = "black", size = 8),
        strip.text.y = element_text(color = "black", size = 8))
ggsave("DS_resolution/Supplementary_AUC_reso.eps", units = "in", width = 4, height = 4,
       dpi = 300)
