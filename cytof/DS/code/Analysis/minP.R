# === tool ========
.libPaths(c(Sys.getenv('R_LIBS_1'), 
            Sys.getenv('R_LIBS_2'),
            Sys.getenv('R_LIBS_3')))
.libPaths()

suppressPackageStartupMessages({
    library(TreeSummarizedExperiment)
    library(treeclimbR)
    library(dplyr)
    library(diffcyt)
})

R.Version()
# ==== arguments from batch R=====================
argsList <- (commandArgs(trailingOnly = TRUE))
argsList
for (i in seq_along(argsList)) {
    if(length(argsList[[i]])>1){
        assign(names(argsList)[i],argsList[[i]])
    }else{eval(parse(text = argsList[[i]]))}
}

## === load data from data preparation step ======
load(inRDat)


## === differential abundance test using edgeR ===
contrast_vec <- c(0, 1, 0, 0, 0, 0, 0, 0, 0)

# set up design matrix
# note: include fixed effects for 'patient_id'
design <- createDesignMatrix(experiment_info, cols_design = 1:2)
design

# set up contrast matrix
contrast <- createContrast(contrast_vec)
contrast

treeR <- rowTree(d_medians_all)
system.time({
    # run limma
    res_tree <- testDS_limma(d_counts_all, d_medians_all, 
                             min_cells = 10, min_samples = 8,
                             design, contrast)
    
    
    # results
    res_DS <- rowData(res_tree) %>%
        data.frame() %>%
        mutate(node = transNode(tree = treeR, 
                                node = as.character(cluster_id))) %>%
        mutate(cluster_id = transNode(tree = treeR, 
                                      node = node)) 
    res_list <- split(res_DS, f = res_DS$marker_id)
    
    # minP
    loc.minP_ab <- lapply(c(0.01, 0.05, 0.1), FUN = function(x){
        xx <- lapply(res_list, FUN = function(y) {
            xy <- getLevel(tree = treeR,
                           score_data = y,
                           drop = p_adj > x,
                           score_column = "p_val",
                           node_column = "node",
                           get_max = FALSE,
                           parent_first = TRUE,
                           message = FALSE)
            xy$node[xy$keep]})
        ll <- unlist(lapply(xx, length)) >0
        xx[ll]
    })
    names(loc.minP_ab) <- c(0.01, 0.05, 0.1)
    
    loc.minP <- lapply(loc.minP_ab, FUN = function(x){
        unique(unlist(x))
    })
})


save(loc.minP, loc.minP_ab, file = outRDat)
