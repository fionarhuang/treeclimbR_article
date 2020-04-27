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
    outMin <- lapply(seq_len(nSIM), FUN = function(x) {
        res <- runDA(TSE = tse, feature_on_row = TRUE, 
                     assay = x, option = "glm", 
                     filter_min_count = 1, 
                     normalize = TRUE, group_column = "group",
                     design_terms = "group")
        out <- nodeResult(res, n = Inf)
        
    out$keep_0.05 <- getLevel(tree = rowTree(tse),
                    score_data = out,
                    drop = FDR > 0.05,
                    score_column = "PValue",
                    node_column = "node",
                    get_max = FALSE,
                    parent_first = TRUE,
                    message = FALSE)$keep
    out$keep_0.01 <- getLevel(tree = rowTree(tse),
                              score_data = out,
                              drop = FDR > 0.01,
                              score_column = "PValue",
                              node_column = "node",
                              get_max = FALSE,
                              parent_first = TRUE,
                              message = FALSE)$keep
    out$keep_0.1 <- getLevel(tree = rowTree(tse),
                              score_data = out,
                              drop = FDR > 0.1,
                              score_column = "PValue",
                              node_column = "node",
                              get_max = FALSE,
                              parent_first = TRUE,
                              message = FALSE)$keep
    return(out)
})
})


save(outMin, file = outRDat)

sessionInfo()