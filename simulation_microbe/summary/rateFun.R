rateFun <- function(loc, nSim, method) {
    res <- vector("list", nSim)
    for (i in seq_len(nSim)) {
        res[[i]] <- lapply(loc, FUN = function(x){x[[i]]})
    }
    rate <- lapply(res, FUN = function(x) {
        xx <- lapply(x, function(y) {
            fdr.y <- fdr(tree = rowTree(tse), truth = truth, 
                         found = y, only.leaf = TRUE)
            tpr.y <- tpr(tree = rowTree(tse), truth = truth, 
                         found = y, only.leaf = TRUE)
            c(fdr.y, tpr.y)
        })
        do.call(rbind, xx)
    })
    
    rate <- lapply(seq_along(rate), FUN = function(x) {
        cbind.data.frame(rate[[x]], alpha = names(loc), 
                         sim = x, method = method)
    })
    rate <- do.call(rbind, rate)
    rownames(rate) <- NULL
    return(rate)
}
