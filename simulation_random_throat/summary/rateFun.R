rateFun <- function(loc, nSim, method, tree, truth) {
    
    res <- vector("list", nSim)
    for (i in seq_len(nSim)) {
        res[[i]] <- lapply(loc, FUN = function(x){x[[i]]})
    }
    
    rate <- lapply(seq_along(res), FUN = function(i) {
        message("Calculating FDR-TPR: ", i, " out of ",
                length(res), " is done")
        res.i <- res[[i]]
        truth.i <- truth[[i]]
        
        rate.i <- lapply(res.i, FUN = function(y) {
            fdr.y <- fdr(tree = tree, truth = truth.i,
                         found = y, only.leaf = TRUE)
            tpr.y <- tpr(tree = tree, truth = truth.i, 
                         found = y, only.leaf = TRUE)
            c(fdr.y, tpr.y)
        })
        do.call(rbind, rate.i)
    })
    
    rate <- lapply(seq_along(rate), FUN = function(x) {
        cbind.data.frame(rate[[x]], alpha = names(loc), 
                         sim = x, method = method)
    })
    rate <- do.call(rbind, rate)
    rownames(rate) <- NULL
    return(rate)
}
