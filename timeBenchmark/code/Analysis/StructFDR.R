# === tool ========
.libPaths(c(Sys.getenv('R_LIBS_1'), 
            Sys.getenv('R_LIBS_2'),
            Sys.getenv('R_LIBS_3')))
.libPaths()

suppressPackageStartupMessages({
    library(TreeSummarizedExperiment)
    library(dplyr)
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
print(inRDat)
load(inRDat)

# ------------------------------- functions ---------------------------------
# wilcox.test
test.func <- function(X, Y) {  
    Y <- as.numeric(factor(Y))
    obj <- apply(X, 1, function(x) {                
        p.value <- suppressWarnings(wilcox.test(x ~ Y)$p.value)
        e.sign <- sign(mean(x[Y == 2]) - mean(x[Y == 1]))
        c(p.value, e.sign)          
    })
    return(list(p.value=obj[1, ], e.sign=obj[2, ])) 
}

# Define permutation function, simple group label permutation
perm.func <- function (X, Y) {  
    return(list(X=X, Y=sample(Y))) 
}

# ------------------------------- data ---------------------------------
# structFDR takes too long to run on tree with 1E4 leaves, so we would not include
# it 
nam <- names(lse)
lse <- lse[!grepl("leaf10000", x = nam)]

X <- lapply(lse, FUN = function(x) {
    assays(x)[[1]]
})
Y <- lapply(lse, FUN = function(x) {
    colData(x)$group
})
treeX <- lapply(lse, rowTree)


# Call TreeFDR
structFDR <- vector("list", length(lse))
names(structFDR) <- names(lse)
for (i in seq_along(lse)) {
    message(i, " out of ", length(lse), " has been done.")
    structFDR[[i]] <- system.time({
        tree.fdr.obj <- TreeFDR(X[[i]], Y[[i]], treeX[[i]],
                                test.func, perm.func)
    })
}

# print out the running time
structFDR

save(structFDR, file = outRDat)

sessionInfo()
