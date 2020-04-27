# === tool ========
.libPaths(c(Sys.getenv('R_LIBS_1'), 
            Sys.getenv('R_LIBS_2'),
            Sys.getenv('R_LIBS_3')))
.libPaths()

suppressPackageStartupMessages({
    library(TreeSummarizedExperiment)
    library(treeclimbR)
    library(dplyr)
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
system.time({
    res_agg <- runDA(TSE = tse, 
                     feature_on_row = TRUE, assay = 1,
                     option = "glm", group_column = "group_id",
                     design_terms = c("patient_id", "group_id"))
    out_tree <-  nodeResult(object = res_agg, n = Inf)
    
    
    
    # get candidate levels
    cand <- getCand(tree = rowTree(tse), 
                    score_data = out_tree,
                    node_column = "node",
                    sign_column = "logFC", 
                    p_column = "PValue",
                    message = FALSE)
     
    
    # find and test on the best level
    limit_rej <- c(0.01, 0.05, 0.1)
    best <- vector("list", length(limit_rej))
    names(best) <- limit_rej
    for (i in seq_along(limit_rej)) {
        best[[i]] <- evalCand(tree = rowTree(tse), 
                              levels = cand$candidate_list,
                              score_data = out_tree, 
                              node_column = "node", 
                              sign_column = "logFC",
                              p_column = "PValue", 
                              method = "BH", 
                              limit_rej = limit_rej[i],
                              use_pseudo_leaf = FALSE)}
    
    
    loc.treeclimbR <- lapply(best, FUN = function(x) {
        xx <- x$output
        xx$node[xx$signal.node %in% TRUE]
    })
    
})

save(out_tree, cand, best, loc.treeclimbR,  file = outRDat)


sessionInfo()