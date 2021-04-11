# === tool ========
.libPaths(c(Sys.getenv('R_LIBS_1'), 
            Sys.getenv('R_LIBS_2'),
            Sys.getenv('R_LIBS_3')))
.libPaths()

suppressPackageStartupMessages({
    library(TreeSummarizedExperiment)
    library(ape)
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


ls.before <- ls()
# ------------- Data on the leaf level ---------------------------
# the number of leaves
sz_tree <- c("leaf100" = 100, "leaf500" = 500, 
             "leaf1000"= 1000, "leaf3000" = 3000,
             "leaf6000" = 6000, "leaf10000" = 10000)

# the number of samples per group
sz_samp <- c("sample50" = 50, "sample250" = 250)

# the library size
sz_lib <- sz_tree*20

# proportions of entities
set.seed(2020)
pr <- lapply(sz_tree, FUN = function(x) {
    # proportions are sum to 1
    rx <- runif(n = x)
    rx <- rx/sum(rx)
    return(rx)
})

# trees
trees <- lapply(sz_tree, rtree)


# count matrices
lse <- vector("list", length(sz_samp))
names(lse) <- names(sz_samp)

for (i in seq_along(sz_samp)) {
    samp_i <- sz_samp[i]
    
    # generate matrices
    lse_i <- lapply(seq_along(sz_tree), FUN = function(j) {
        lib_j <- sz_lib[j]
        xx <- rmultinom(n = samp_i*2, size = lib_j, prob = pr[[j]])
        
        # select 30% entities, multiply their counts by 8 in the first 50% samples
        # To test the time, this step is not mandatory 
        # but we have it to generate some differences between groups
        swp <- seq_len(0.3*sz_tree[j])
        xx[swp, seq_len(samp_i)] <- xx[swp, seq_len(samp_i)] * 8
        rownames(xx) <- trees[[j]]$tip.label
        colnames(xx) <- paste0("sample_", seq_len(ncol(xx)))
        
        # generate information of samples
        gr <- data.frame(group = rep(LETTERS[1:2], each = samp_i),
                         row.names = colnames(xx))
        
        dse <- TreeSummarizedExperiment(assays = list(xx),
                                        colData = gr,
                                        rowTree = trees[[j]])
        return(dse)
    })
    names(lse_i) <- names(sz_tree)
    lse[[i]] <- lse_i
}

lse <- unlist(lse, recursive = FALSE)
names(lse)
names(lse) <- gsub(pattern = "[.]", "_", names(lse))
names(lse)

Time_agg <- tse <- setNames(vector("list", length(lse)), names(lse))
for (i in seq_along(lse)) {
    message(i, " out of ", length(lse), " are done")
    lse_i <- lse[[i]]
    Time_agg[[i]] <- system.time({
        nodes <- showNode(tree = rowTree(lse_i), only.leaf = FALSE)
        tse[[i]] <- aggValue(x = lse_i, rowLevel = nodes, FUN = sum)
    })
}
    


ls.after <- ls()
ls.save <- setdiff(ls.after, c(ls.before, "ls.before"))
# if the file doesn't exist or was modified earlier than the R script
if((!file.exists(outRDat)) | (file.mtime(outRDat) < file.mtime(scriptP))){
    save(list = ls.save, file = outRDat)
}



sessionInfo()








