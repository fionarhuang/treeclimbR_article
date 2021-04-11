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
minP <- vector("list", length(tse))
names(minP) <- names(tse)
for (i in seq_along(tse)) {
    message(i, " out of ", length(tse), " is done.")
    minP[[i]] <- system.time({
        resW <- test.func(X[[i]],Y[[i]])
        outW <- data.frame(node = rowLinks(tse[[i]])$nodeNum,
                           pvalue = resW$p.value,
                           sign = resW$e.sign) %>%
            mutate(adjP = p.adjust(pvalue, "BH"))
        
        outW$keep_0.05 <- getLevel(tree = treeX[[i]],
                                  score_data = outW,
                                  drop = adjP > 0.05,
                                  node_column = "node",
                                  score_column = "adjP",
                                  get_max = FALSE,
                                  parent_first = TRUE,
                                  message = FALSE)$keep
    })
    
}


minP_all <- mapply(FUN = function(x, y) {
    x + y
}, Time_agg, minP, SIMPLIFY = FALSE)

# print out the run time
minP_all
minP
save(minP, minP_all, file = outRDat)

sessionInfo()