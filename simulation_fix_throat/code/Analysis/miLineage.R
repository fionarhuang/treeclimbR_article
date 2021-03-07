# === tool ========
.libPaths(c(Sys.getenv('R_LIBS_1'), 
            Sys.getenv('R_LIBS_2'),
            Sys.getenv('R_LIBS_3')))
.libPaths()

suppressPackageStartupMessages({
    library(edgeR)
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


print(inRDat)

## === load data from data preparation step ===
load(inRDat)
limFDR <- 0.05
# ================= prepare data ====================
# tree
Tree <- rowTree(tse)
mTree <- matTree(tree = Tree)
rownames(mTree) <- transNode(tree = Tree, node = mTree[, 1], use.alias = TRUE)
rmTree <- apply(mTree, 1, FUN = function(x){
    lx <- length(x)
    lv <- rev(x[!is.na(x)])
    la <- transNode(tree = Tree, node = lv, use.alias = TRUE,
                    message = FALSE)
    fv <- c(la, rep("Unclassified", lx - length(lv)))
    fv
})
taxTree <- t(rmTree)
colnames(taxTree) <- paste("Rank", seq_len(ncol(taxTree)), sep = "")
# count table
tipDAT <- tse[rowLinks(tse)$isLeaf, ]
tipCount <- assays(tipDAT, withDimnames = TRUE)
tipCount <- lapply(tipCount, FUN = function(x) {
    rownames(x) <- rowLinks(tipDAT)$nodeLab
    return(x)
})
taxCount <- lapply(tipCount, 
                   FUN = function(x){
                       xt <-t(x)
                       cn <- transNode(tree = rowTree(tse), node = colnames(xt))
                       colnames(xt) <- transNode(tree = rowTree(tse), 
                                                 node = cn,
                                                 use.alias = TRUE)
                       xt <- xt[, rownames(taxTree)]
                       return(xt)
                   })

# case
cv <- gsub(pattern = "C", "", colData(tse)$group)
cv <- factor(cv)
case <- matrix(as.numeric(cv), ncol = 1,
               dimnames = list(rownames(taxCount[[1]]), "case"))

# check OTU order
all(rownames(taxTree) == colnames(taxCount[[1]]))
# check sample order
all(rownames(case) == rownames(taxCount[[1]]))

# analysis (resample 100) find nothing
# one-part analysis
# system.time({
# loc1_rs100_0.05.MLA <- lapply(seq_along(taxCount), FUN = function(x){
#     res <- QCAT(taxCount[[x]], case, 1, taxTree,
#                 fdr.alpha = limFDR, n.resample = 100)
#     vres <- res$sig.lineage
#     cat(x, "out of ", length(taxCount), "has been done", "\n")
#     return(vres)
# })
# })
# # two - part analysis
# system.time({
# loc2_rs100_0.05.MLA <- lapply(seq_along(taxCount), FUN = function(x){
# 
#     res <- QCAT_GEE(taxCount[[x]], case, 1, case, 1, taxTree,
#                     fdr.alpha = limFDR, n.resample = 100)
#     vres <- res$sig.lineage
#     ures <- unique(unlist(vres))
#     cat(x, "out of ", length(taxCount), "has been done", "\n")
#     return(ures)
# })
# })

# analysis 
# one-part analysis
system.time({
    loc1_0.05.MLA <- lapply(seq_along(taxCount), FUN = function(x){
    res <- QCAT(taxCount[[x]], case, 1, taxTree,
                fdr.alpha = limFDR)
    vres <- res$sig.lineage
    ures <- unique(unlist(vres))
    ures <- transNode(tree = Tree, node = ures)
    cat(x, "out of ", length(taxCount), "has been done", "\n")
    return(ures)
})
})

# two - part analysis
system.time({
    loc2_0.05.MLA <- lapply(seq_along(taxCount), FUN = function(x){
    res <- QCAT_GEE(taxCount[[x]], case, 1, case, 1, taxTree,
                    fdr.alpha = limFDR)
    vres <- res$sig.lineage
    ures <- unique(unlist(vres))
    ures <- transNode(tree = Tree, node = ures)
    cat(x, "out of ", length(taxCount), "has been done", "\n")
    return(ures)
})
})

# ------------------0.01 -------------------------------------------
system.time({
    loc1_0.01.MLA <- lapply(seq_along(taxCount), FUN = function(x){
        res <- QCAT(taxCount[[x]], case, 1, taxTree,
                    fdr.alpha = 0.01)
        vres <- res$sig.lineage
        ures <- unique(unlist(vres))
        ures <- transNode(tree = Tree, node = ures)
        cat(x, "out of ", length(taxCount), "has been done", "\n")
        return(ures)
    })
})

# two - part analysis
system.time({
    loc2_0.01.MLA <- lapply(seq_along(taxCount), FUN = function(x){
        res <- QCAT_GEE(taxCount[[x]], case, 1, case, 1, taxTree,
                        fdr.alpha = 0.01)
        vres <- res$sig.lineage
        ures <- unique(unlist(vres))
        ures <- transNode(tree = Tree, node = ures)
        cat(x, "out of ", length(taxCount), "has been done", "\n")
        return(ures)
    })
})

# ------------------0.1 -------------------------------------------
system.time({
    loc1_0.1.MLA <- lapply(seq_along(taxCount), FUN = function(x){
        res <- QCAT(taxCount[[x]], case, 1, taxTree,
                    fdr.alpha = 0.1)
        vres <- res$sig.lineage
        ures <- unique(unlist(vres))
        ures <- transNode(tree = Tree, node = ures)
        cat(x, "out of ", length(taxCount), "has been done", "\n")
        return(ures)
    })
})

# two - part analysis
system.time({
    loc2_0.1.MLA <- lapply(seq_along(taxCount), FUN = function(x){
        res <- QCAT_GEE(taxCount[[x]], case, 1, case, 1, taxTree,
                        fdr.alpha = 0.1)
        vres <- res$sig.lineage
        ures <- unique(unlist(vres))
        ures <- transNode(tree = Tree, node = ures)
        cat(x, "out of ", length(taxCount), "has been done", "\n")
        return(ures)
    })
})



ls.mla <- ls()
# save result
ind <- grepl(".MLA$", ls.mla)
obj.mla <- ls.mla[ind]
save(list=obj.mla,file=outRDat)

sessionInfo()