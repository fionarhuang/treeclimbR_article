# === tool ========
.libPaths(c(Sys.getenv('R_LIBS_1'), 
            Sys.getenv('R_LIBS_2'),
            Sys.getenv('R_LIBS_3')))
.libPaths()

suppressPackageStartupMessages({
    library(TreeSummarizedExperiment)
    library(treeclimbR)
    library(dplyr)
    library(edgeR)
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


contrast_vec <- c(0, 1, 0, 0, 0, 0, 0, 0, 0)

# set up design matrix
# note: include fixed effects for 'patient_id'
design <- createDesignMatrix(experiment_info, cols_design = 1:2)
design

# set up contrast matrix
contrast <- createContrast(contrast_vec)
contrast

# run tests
d_counts_leaf <- d_counts_all[rowLinks(d_counts_all)$isLeaf, ]
rowData(d_counts_leaf) <- droplevels(rowData(d_counts_leaf))
d_medians_leaf <- d_medians_all[rowLinks(d_medians_all)$isLeaf, ]
rowData(d_medians_leaf) <- droplevels(rowData(d_medians_leaf))

DS_leaf <- testDS_limma(d_counts = d_counts_leaf, 
                            d_medians = d_medians_leaf, 
                            design, contrast)



# top DA clusters
res_leaf <- diffcyt::topTable(DS_leaf, format_vals = TRUE, 
                              top_n = nrow(DS_leaf))

# number of significant DA clusters at 1% FDR
rej_limit <- c(0.01, 0.05, 0.1)

loc.diffcyt <- lapply(c(0.01, 0.05, 0.1), FUN = function(x) {
    xx <- res_leaf$cluster_id[res_leaf$p_adj <= x]
    xx <- xx[!is.na(xx)]
    unique(as.character(xx))
    
})

treeR <- rowTree(d_medians_all)
loc.diffcyt_ab <- lapply(c(0.01, 0.05, 0.1), FUN = function(x) {
    xx <- res_leaf %>%
        data.frame() %>%
        mutate(node = transNode(tree = treeR, node = as.character(cluster_id))) %>%
        filter(p_adj <= x) %>%
        mutate(marker_id = as.character(marker_id))
    split(xx$node, xx$marker_id)
    
})
names(loc.diffcyt) <- names(loc.diffcyt_ab) <- rej_limit

save(loc.diffcyt, res_leaf, loc.diffcyt_ab, file = outRDat)
