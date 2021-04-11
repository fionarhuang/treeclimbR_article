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
X <- lapply(tse, FUN = function(x) {
    assays(x)[[1]]
})
Y <- lapply(tse, FUN = function(x) {
    colData(x)$group
})
treeX <- lapply(tse, rowTree)

# ------------------------------- run time ---------------------------------
treeclimbR <- vector("list", length(tse))
names(treeclimbR) <- names(tse)
for (i in seq_along(tse)) {
    message(i, " out of ", length(tse), " is done.")
    treeclimbR[[i]] <- system.time({
        resW <- test.func(X[[i]],Y[[i]])
        outW <- data.frame(node = rowLinks(tse[[i]])$nodeNum,
                           pvalue = resW$p.value,
                           sign = resW$e.sign)
        cand <- getCand(tree = treeX[[i]], score_data = outW, 
                        node_column = "node", p_column = "pvalue",
                        threshold = 0.05,
                        sign_column = "sign", message = TRUE)
        best <- evalCand(tree = treeX[[i]], levels = cand$candidate_list, 
                         score_data = outW, node_column = "node",
                         p_column = "pvalue", sign_column = "sign",
                         message = FALSE)
        outB <- topNodes(object = best, n = Inf, p_value = 0.05)
    })
}


treeclimbR_all <- mapply(FUN = function(x, y) {
    x + y
}, Time_agg, treeclimbR, SIMPLIFY = FALSE)

# print out the run time
treeclimbR
treeclimbR_all

save(treeclimbR, treeclimbR_all, file = outRDat)


sessionInfo()