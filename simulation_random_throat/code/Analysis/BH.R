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
(nSIM <- length(tse_list))

system.time({
    lse_list <- lapply(tse_list, FUN = function(x) {
        x[rowLinks(x)$isLeaf, ]})
    
    out <- lapply(seq_len(nSIM), 
                  FUN = function(x) {
                      cat(x, "out of ", nSIM, "has been done", "\n")
                      res <- runDA(TSE = lse_list[[x]], 
                                   feature_on_row = TRUE, 
                                   assay = 1, option = "glm", 
                                   filter_min_count = 1,
                                   normalize = TRUE, 
                                   group_column = "group",
                                   design_terms = "group")
                      out <- nodeResult(res, n = Inf) 
                      })
    loc.bh_0.05 <- lapply(seq_len(nSIM), FUN = function(x) {
        out.i <- out[[x]]
        out.i$node[out.i$FDR <= 0.05]
    })
    
    loc.bh_0.01 <- lapply(seq_len(nSIM), FUN = function(x) {
        out.i <- out[[x]]
        out.i$node[out.i$FDR <= 0.01]
    })
    
    loc.bh_0.1 <- lapply(seq_len(nSIM), FUN = function(x) {
        out.i <- out[[x]]
        out.i$node[out.i$FDR <= 0.1]
    })
})


save(loc.bh_0.01, loc.bh_0.05, loc.bh_0.1, file = outRDat)

sessionInfo()
