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

# ------------------------------- data ---------------------------------
X <- lapply(lse, FUN = function(x) {
    assays(x)[[1]]
})
Y <- lapply(lse, FUN = function(x) {
    colData(x)$group
})
treeX <- lapply(lse, rowTree)

# ------------------------------- run time ---------------------------------
BH <- vector("list", length(lse))
names(BH) <- names(lse)

for (i in seq_along(lse)) {
    message(i, " out of ", length(lse), " is done.")
    BH[[i]] <- system.time({
        resW <- test.func(X[[i]],Y[[i]])
        adjP <- p.adjust(p = resW$p.value, method = "BH")
        names(adjP)[adjP <= 0.05]
    })
}

# print out the run time
BH

save(BH, file = outRDat)


sessionInfo()
