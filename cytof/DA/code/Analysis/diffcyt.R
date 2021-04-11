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
    library(StructFDR)
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


# set up design matrix
# note: include fixed effects for 'patient_id'
design <- createDesignMatrix(experiment_info, cols_design = 1:2)
design

# set up contrast (entries corresponding to columns of design matrix)
contrast <- createContrast(c(0, 1, 0, 0, 0, 0))
contrast


# run tests
# using default filtering
out_DA <- testDA_edgeR(d_counts, design, contrast)

# show results
rowData(out_DA)

# top DA clusters
top_DA <- diffcyt::topTable(out_DA, format_vals = TRUE, top_n = nrow(d_counts))

# number of significant DA clusters at 10% FDR
loc.diffcyt <- lapply(c(0.01, 0.05, 0.1), FUN = function(x) {
    xx <- top_DA$cluster_id[top_DA$p_adj <= x]
    xx <- xx[!is.na(xx)]
    xx <- as.character(xx)
    transNode(tree = rowTree(tse), node = xx)
})
names(loc.diffcyt) <- c(0.01, 0.05, 0.1)

save(loc.diffcyt, file = outRDat)


sessionInfo()
