# === tool ========
.libPaths(c(Sys.getenv('R_LIBS_1'), 
            Sys.getenv('R_LIBS_2'),
            Sys.getenv('R_LIBS_3')))
.libPaths()

suppressPackageStartupMessages({
    library(TreeSummarizedExperiment)
    library(treeclimbR)
    library(dplyr)
    library(structSSI)
    library(ape)
    library(igraph)
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

# ------------------------------- data ---------------------------------
X <- lapply(lse, FUN = function(x) {
    assays(x)[[1]]
})
Y <- lapply(lse, FUN = function(x) {
    colData(x)$group
})
treeX <- lapply(lse, rowTree)


HFDR <- vector("list", length(lse))
names(HFDR) <- names(lse)

for (i in seq_along(lse)) {
    message(i, " out of ", length(lse), " has been done")
    HFDR[[i]] <- system.time({
        chl.tree <- get.edgelist(as.igraph(treeX[[i]]))
        chl.pval <- treePValues(chl.tree, X[[i]], Y[[i]])
        # the root has NA p value because all samples have the same library size
        # to run next step successfully, p values is set to 0.5 because 
        # there are no difference between two groups
        chl.pval[is.na(chl.pval)] <- 0.5
        chl.hfdr <- hFDR.adjust(chl.pval, chl.tree, alpha = 0.05)
        
        # Discovery on the leaf level
        find.hFDR <- rownames(chl.hfdr@p.vals)[chl.hfdr@p.vals$adjp <= 0.05]
        find.hFDR <- find.hFDR[!is.na(find.hFDR)]
        leaf.hFDR <- find.hFDR[startsWith(find.hFDR, "t")]  
    })
    
}

# print out the run time
HFDR

save(HFDR, file = outRDat)


sessionInfo()