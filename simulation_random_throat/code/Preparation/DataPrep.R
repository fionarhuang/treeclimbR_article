# === tool ========
.libPaths(c(Sys.getenv('R_LIBS_1'), 
            Sys.getenv('R_LIBS_2'),
            Sys.getenv('R_LIBS_3')))
.libPaths()

suppressPackageStartupMessages({
    library(TreeSummarizedExperiment)
    library(treeclimbR)
})

# =========== source TreeInfer package ===================
R.Version()

# ==== arguments from batch R===================
argsList <- (commandArgs(trailingOnly = TRUE))
argsList
for (i in seq_along(argsList)) {
    if(length(argsList[[i]])>1){
        assign(names(argsList)[i],argsList[[i]])
    }else{eval(parse(text = argsList[[i]]))}
}

# ========== read in parameter value ================#
print(inRDat)
load(inRDat)

set.seed(seedNum)

data("throat_v35")
throat_v35 <- parEstimate(obj = throat_v35)
lib <- apply(assays(throat_v35)[[1]], 2, sum)


# ---------------------------------------------------------------------------
# As the simulation is based on real data, it might not exist for two branches that fullfill the restricted requirements (e.g., branch size, fold changes) to swap; or it might exist too many possibilities under the specific requirement
# Also, a tree structure naturally has more nodes with small descendant leaves than nodes with large descendant leaves. In this simulation, we would like to have signal branches at different sizes. 
# Therefore, we use 'selNode' before 'simData' to speed up the simulation and to specificy the number of simulations at different size ranges of signal branches.
# ---------------------------------------------------------------------------

lev <- cbind(seq(from = 1, by = 5, length.out = 10),
             seq(from = 5, by = 5, length.out = 10))
nr <- nrow(lev)
# the number of nodes in the tree has the corresponding branch sizes 
sapply(seq_len(nr), FUN = function(x) {
    xx <- selNode(obj = throat_v35, minTip = lev[x, 1], 
                  maxTip = lev[x, 2],  all = TRUE)
    nrow(xx)
})

# To speed up: randomly select 100 nodes
selA <- sapply(seq_len(nr), FUN = function(x) {
    xx <- selNode(obj = throat_v35, minTip = lev[x, 1], 
                  maxTip = lev[x, 2],  all = TRUE)
    sample(xx$nodeNum, 10, replace = TRUE)
})


lse_list <- vector("list", nSIM)

i = 1
while (i <= nSIM) {
    message(i)
    
    # select a ratio randomly from the range below
    fcRange <- seq(from = 2, to = 8, by = 0.01)
    fc <- sample(fcRange, 1)
    
    # randomly select a branch from 'selA', and then randomly select another
    # branches to swap. The swap should cause a fold change close to 'fc' that
    # is randomly selected from 'fcRange'.
    # Normally the true ratio is not exactly equal to the specified 'fc'. 
    lse <- simData(obj = throat_v35,
                   from.A = sample(selA, 1),
                   ratio = fc,
                   mu = lib,
                   nSam = nSam,
                   n = 1,
                   scenario = scene,
                   pct = pct)
    
    # the true ratio (for BS & SS)
    true_ratio <- metadata(lse)$branch$ratio
    
    # for US:
    # 1) requires the abundance ratio of two swapped branches is above 1.5 to 
    #   have reasonable difference between two groups
    if (scene %in% "US") {
        FC <- metadata(lse)$FC
        if (true_ratio > 1.5) {
            lse_list[[i]] <- lse
            i <- i + 1
        }
    }
    
    # for SS:
    # 1) requires the fold change of two branches to be above 1.5 to have 
    #    reasonable difference between two groups
    if (scene %in% c("BS", "SS")) {
        FC <- metadata(lse)$FC
        uFC <- unique(FC[FC>1])
        if (uFC > 1.5) {
            lse_list[[i]] <- lse
            i <- i + 1
        }
    }
}

allNode <- showNode(tree = rowTree(throat_v35), only.leaf = FALSE)
tse_list <- lapply(lse_list, FUN = function(x) {
    aggValue(x = x, rowLevel = allNode, FUN = sum)
})

# if the file doesn't exist or was modified earlier than the R script
if((!file.exists(outRDat)) | (file.mtime(outRDat) < file.mtime(scriptP))){
    save(tse_list, file = outRDat)
}



sessionInfo()








