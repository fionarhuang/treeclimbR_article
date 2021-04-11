# === tool ========
.libPaths(c(Sys.getenv('R_LIBS_1'), 
            Sys.getenv('R_LIBS_2'),
            Sys.getenv('R_LIBS_3')))
.libPaths()

suppressPackageStartupMessages({
    library(TreeSummarizedExperiment)
    library(miLineage)
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




## === load data from data preparation step ===
print(inRDat)
load(inRDat)

# ------------------------------- data ---------------------------------
# tree
treeX <- lapply(lse, rowTree)
taxTree <- lapply(seq_along(treeX), FUN = function(i) {
    message(i, " out of ", length(treeX), " is done")
    
    y <- treeX[[i]]
    xx <- matTree(tree = y)
    rownames(xx) <- transNode(tree = y, node = xx[, 1],
                              use.alias = TRUE)
    rmTree <- apply(xx, 1, FUN = function(x){
        lx <- length(x)
        lv <- rev(x[!is.na(x)])
        la <- transNode(tree = y, node = lv, use.alias = TRUE,
                        message = FALSE)
        fv <- c(la, rep("Unclassified", lx - length(lv)))
        fv
    })
    taxTree <- t(rmTree)
    colnames(taxTree) <- paste("Rank", seq_len(ncol(taxTree)), sep = "")
    return(taxTree)
    })


# count table
taxCount <- lapply(seq_along(lse), FUN = function(i) {
    message(i, " out of ", length(lse), " is done")
    x <- lse[[i]]
    xx <- assays(x)[[1]]
    rownames(xx) <- rowLinks(x)$nodeLab
    
    xt <-t(xx)
    cn <- transNode(tree = rowTree(x), node = colnames(xt))
    colnames(xt) <- transNode(tree = rowTree(x), 
                              node = cn,
                              use.alias = TRUE)
    xt <- xt[, rownames(taxTree[[i]])]
    return(xt)
})


# case
case <- lapply(seq_along(lse), FUN = function(i) {
    message(i, " out of ", length(lse), " is done")
    x <- lse[[i]]
    gr <- as.factor(colData(x)$group)
    matrix(as.numeric(gr), ncol = 1,
           dimnames = list(rownames(taxCount[[i]]), "case"))
})


# check OTU order
f1 <- function(x, y) {all(rownames(x) == colnames(y))}
f2 <- function(x, y) {all(rownames(x) == rownames(y))}
all(mapply(FUN = f1, taxTree, taxCount))

# check sample order
all(mapply(FUN = f2, case, taxCount))


# analysis 
# one-part analysis
miLineage1 <- vector("list", length(lse))
names(miLineage1) <- names(lse)
miLineage2 <- miLineage1

for (i in seq_along(taxCount)) {
    message(i, "out of ", length(taxCount), "has been done", "\n")
    miLineage1[[i]] <- system.time({
        res <- QCAT(taxCount[[i]], case[[i]], 1, taxTree[[i]],
                    fdr.alpha = 0.05)
        vres <- res$sig.lineage
        ures <- unique(unlist(vres))
    })
}

# two - part analysis

for (i in seq_along(taxCount)) {
    message(i, "out of ", length(taxCount), "has been done", "\n")
    miLineage2[[i]] <- system.time({
        res <- QCAT_GEE(taxCount[[i]], case[[i]], 1, case[[i]], 1, taxTree[[i]],
                        fdr.alpha = 0.05)
        vres <- res$sig.lineage
        ures <- unique(unlist(vres))
    })
}

# print out the run time
miLineage1
miLineage2

save(miLineage1, miLineage2, file=outRDat)

sessionInfo()