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

print(dataFile)
d_SE <- readRDS(dataFile)



# -------------------------------------------------------------
# count table
# -------------------------------------------------------------

experiment_info <- metadata(d_SE)$experiment_info 

# keep original data
d_SE_original <- d_SE

# check levels of 'sample_id' column in row data are in correct order
rowData(d_SE)
levels(rowData(d_SE)$sample_id)

# split data (note: using 'split.data.frame' to split into list of matrices; see '?split')
d_input <- split.data.frame(assay(d_SE), rowData(d_SE)$sample_id)


marker_info <- as.data.frame(colData(d_SE))
rownames(marker_info) <- NULL  
marker_info

d_se <- prepareData(d_input, experiment_info, marker_info)

# check marker classes
colnames(d_se)[colData(d_se)$marker_class == "type"]
colnames(d_se)[colData(d_se)$marker_class == "state"]

d_se <- transformData(d_se, cofactor = 5)


# run clustering (using FlowSOM with 20x20 grid, i.e. 400 clusters)
d_se <- generateClusters(d_se, xdim = nx, ydim = ny, seed_clustering = seedN)
rowData(d_se)$spikein <- rowData(d_SE_original)$spikein

#
# calculate cluster cell counts
d_counts <- calcCounts(d_se[rowData(d_se)$group_id %in% condition, ])
d_counts <- d_counts[,colData(d_counts)$group_id %in% condition]
experiment_info <- experiment_info %>%
    filter(group_id %in% condition) %>%
    droplevels()


# -------------------------------------------------------------
# meta data
# -------------------------------------------------------------

# calculate number and proportion of true spike-in cells per cluster
# number of cells per sample (including spike-in cells)
n_cells <- table(rowData(d_SE_original)$sample_id)

# add spike-in status to rowData of 'd_se' object
stopifnot(nrow(d_se) == nrow(d_SE_original))
stopifnot(all(rowData(d_se)$sample_id == rowData(d_SE_original)$sample_id))

rowData(d_se)$spikein <- rowData(d_SE_original)$spikein
rowData(d_se)



# -------------------------------------------------------------
# create the tree
# -------------------------------------------------------------
cytof_tree <- buildTree(d_se = d_se)


# ------------------------------------------------------------
# information about clusters on the leaf
# ------------------------------------------------------------

df_truth <- rowData(d_se) %>%
    as.data.frame %>%
    mutate(spikein_CN = (group_id %in% "CN") & spikein) %>%
    mutate(spikein_CBF = (group_id %in% "CBF") & spikein) %>%
    group_by(cluster_id) %>%
    summarize(n_cells = n(), 
              n_CN = sum(spikein_CN),
              n_CBF = sum(spikein_CBF),
              n_healthy = sum(!(spikein_CBF|spikein_CN)),
              prop_CN = mean(spikein_CN), 
              prop_CBF = mean(spikein_CBF)) %>%
    mutate(node = transNode(tree = cytof_tree, node = as.character(cluster_id)))



# -------------------------------------------------------------
# create the TreeSummarizedExperiment
# -------------------------------------------------------------

# counts
tse <- calcTreeCounts(d_se = d_se, tree = cytof_tree)
tse <- tse[, colData(tse)$group_id %in% condition]

colData(tse) <- colData(tse) %>%
    as.data.frame() %>%
    select(patient_id, group_id) %>%
    droplevels() %>%
    DataFrame()
tse
# if the file doesn't exist or was modified earlier than the R script
if((!file.exists(outRDat)) | (file.mtime(outRDat) < file.mtime(scriptP))){
    save(tse, df_truth, d_counts, experiment_info, file = outRDat)
}












