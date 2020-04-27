# === tool ========
source("/home/fiona/phd/microbes/simulation/Tool/argsR_compile.R")

.libPaths(c(Sys.getenv('R_LIBS_1'), 
            Sys.getenv('R_LIBS_2'),
            .libPaths()))
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
args <- (commandArgs(trailingOnly = TRUE))
args
argsList <- argRun(args, grp.pattern = "dirGRP")
argsList
for (i in seq_along(argsList)) {
    if(length(argsList[[i]])>1){
        assign(names(argsList)[i],argsList[[i]])
    }else{eval(parse(text = argsList[[i]]))}
}

## === load data from data preparation step ======
load(inRDat)

## === differential abundance test using edgeR ===
(nSIM <- length(assays(tse)))

system.time({
    lse <- tse[rowLinks(tse)$isLeaf, ]
    out <- lapply(seq_len(nSIM), 
                  FUN = function(x) {
                      res <- runDA(TSE = lse, 
                                   feature_on_row = TRUE, 
                                   assay = x, option = "glm", 
                                   filter_min_count = 0,
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
