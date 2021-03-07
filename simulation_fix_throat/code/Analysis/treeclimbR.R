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
(nSIM <- length(assays(tse)))
system.time({
    # data analysis: get p value
    message("Calculating p values ...")
    out <- lapply(seq_len(nSIM), FUN = function(x) {
        cat(x, "out of ", nSIM, "has been done", "\n")
        res <- runDA(TSE = tse, feature_on_row = TRUE, 
                     assay = x, option = "glm", 
                     filter_min_count = 1, 
                     normalize = TRUE, group_column = "group",
                     design_terms = "group")
        out <- nodeResult(object = res, n = Inf)
        return(out)
    })
    
    # get candidate levels
    cand <- lapply(seq_along(out), FUN = function(x) {
        cat("searching for candidates: ", x,
            " out of ", nSIM, " simulations ...\n")
        lev <- getCand(tree = rowTree(tse), 
                       score_data = out[[x]],
                       node_column = "node",
                       sign_column = "logFC", 
                       p_column = "PValue",
                       message = FALSE)
        return(lev)
    })
    
    # find and test on the best level
    limit_rej <- c(0.01, 0.05, 0.1)
    best <- vector("list", length(limit_rej))
    names(best) <- limit_rej
    for (i in seq_along(limit_rej)) {
        best[[i]] <- lapply(seq_along(out),
                            FUN = function(x) {
            cat("get the best: ", x, " out of ", 
                nSIM, " simulations ...\n")
            xx <- evalCand(tree = rowTree(tse), 
                          levels = cand[[x]]$candidate_list,
                          score_data = out[[x]], 
                          node_column = "node", 
                          sign_column = "logFC",
                          p_column = "PValue", 
                          method = "BH", 
                          limit_rej = limit_rej[i], 
                          use_pseudo_leaf = FALSE)
            return(xx)
        })
    }
    
    
   
outsel_0.01 <- lapply(best$`0.01`, FUN = function(x) {
    x$output })

outsel_0.05 <- lapply(best$`0.05`, FUN = function(x) {
    x$output })

outsel_0.1 <- lapply(best$`0.1`, FUN = function(x) {
    x$output })
})

save(out, cand, best, outsel_0.01, outsel_0.05, outsel_0.1, 
     file = outRDat)

sessionInfo()