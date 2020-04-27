# === tool ========
.libPaths(c(Sys.getenv('R_LIBS_1'), 
            Sys.getenv('R_LIBS_2'),
            Sys.getenv('R_LIBS_3')))
.libPaths()

suppressPackageStartupMessages({
    library(TreeSummarizedExperiment)
    library(treeclimbR)
    library(diffcyt)
    library(dplyr)
    library(ape)
})



# =========== source TreeInfer package ===================
R.Version()

# ==== arguments from batch R===================
argsList <- (commandArgs(trailingOnly = TRUE))
argsList
for (i in seq_along(argsList)) {
    if(length(argsList[[i]])>1){
        assign(names(argsList)[i],argsList[[i]])
    }else{eval(parse(text = argsList[[i]]))}
}
# ========== read in parameter value ================#
print(inRDat)
load(inRDat)

d_SE <- readRDS(dataFile)



# -------------------------------------------------------------
# count table
# -------------------------------------------------------------

experiment_info <- metadata(d_SE)$experiment_info 
# split expression table into one table per sample
# (note: will update diffcyt 'prepareData' function so this can be done automatically)

# keep original data
d_SE_original <- d_SE

# check levels of 'sample_id' column in row data are in correct order
rowData(d_SE)
levels(rowData(d_SE)$sample_id)

# split data (note: using 'split.data.frame' to split into list of matrices; see '?split')
d_input <- split.data.frame(assay(d_SE), rowData(d_SE)$sample_id)


marker_info <- as.data.frame(colData(d_SE))
rownames(marker_info) <- NULL  # (note: will update diffcyt 'prepareData' so this line is not required)
marker_info

d_se <- prepareData(d_input, experiment_info, marker_info)

# check marker classes
colnames(d_se)[colData(d_se)$marker_class == "type"]
colnames(d_se)[colData(d_se)$marker_class == "state"]

d_se <- transformData(d_se, cofactor = 5)


# run clustering (using FlowSOM with 20x20 grid, i.e. 400 clusters)
d_se <- generateClusters(d_se, xdim = nx, ydim = ny, seed_clustering = seedN)
rowData(d_se) <- rowData(d_se) %>%
    data.frame() %>%
    mutate(spikein = rowData(d_SE_original)$spikein) %>%
    mutate(B_cell = rowData(d_SE_original)$B_cell) %>%
    mutate(population_id = rowData(d_SE_original)$population_id)


# # calculate cluster cell counts
# d_counts <- calcCounts(d_se)
# experiment_info <- experiment_info %>%
#     droplevels()

# 
# # -------------------------------------------------------------
# # meta data
# # -------------------------------------------------------------
# 
# # calculate number and proportion of true spike-in cells per cluster
# # number of cells per sample (including spike-in cells)
# n_cells <- table(rowData(d_SE_original)$sample_id)
# 
# # add spike-in status to rowData of 'd_se' object
# stopifnot(nrow(d_se) == nrow(d_SE_original))
# stopifnot(all(rowData(d_se)$sample_id == rowData(d_SE_original)$sample_id))
# 
# 
df_truth <- rowData(d_se) %>%
    as.data.frame %>%
    group_by(cluster_id) %>%
    summarize(n_cells = n(),
              prop_spikein = mean(spikein), 
              prop_B = mean(B_cell)) %>%
    as.data.frame()
cell_info <- rowData(d_se)
# # -------------------------------------------------------------
# create the tree
# -------------------------------------------------------------
cytof_tree <- buildTree(d_se = d_se)

# -------------------------------------------------------------
# calculate medians, counts for each cluster
# -------------------------------------------------------------
# medians
d_medians_all <- calcTreeMedians(d_se = d_se, tree = cytof_tree)

# counts
d_counts_all <- calcTreeCounts(d_se = d_se, tree = cytof_tree)
    
# -------------------------------------------------------------
# create the TreeSummarizedExperiment
# -------------------------------------------------------------
# contrast_vec <- c(0, 1, 0, 0, 0, 0, 0, 0, 0)
# 
# # set up design matrix
# # note: include fixed effects for 'patient_id'
# design <- createDesignMatrix(experiment_info, cols_design = 1:2)
# design
# 
# # set up contrast matrix
# contrast <- createContrast(contrast_vec)
# contrast
# 
# # run tests
# res <- testDS_limma(d_counts_all, d_medians_all, design, contrast)
# 
# 
# # -------------------------------------------------------------
# # create the TreeSummarizedExperiment
# # -------------------------------------------------------------
# 
# 
# # -------------------------------------------------------------
# # create the TreeSummarizedExperiment
# # -------------------------------------------------------------
# 
# allNode <- unique(sort(as.vector(rowTree(lse)$edge)))
# tse <- aggValue(x = lse, rowLevel = allNode, FUN = sum)
# tse
# if the file doesn't exist or was modified earlier than the R script
if((!file.exists(outRDat)) | (file.mtime(outRDat) < file.mtime(scriptP))){
    save(d_se, d_counts_all, d_medians_all, 
         df_truth, experiment_info, file = outRDat)
}












