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
load(inRDat)

## === differential abundance test using edgeR ===
(nSIM <- length(tse_list))
rejLimit <- 0.05

treeR <- rowTree(tse_list[[1]])
ed <- treeR$edge
system.time({
    
    # data analysis: get p value
    out <- lapply(seq_len(nSIM), FUN = function(x) {
        cat(x, "out of ", nSIM, "has been done", "\n")
        res <- runDA(TSE = tse_list[[x]], feature_on_row = TRUE, 
                     assay = 1, option = "glm", 
                     filter_min_count = 0,  # turn off filtering; otherwise HFDR gets error due to NA value
                     normalize = TRUE, group_column = "group",
                     design_terms = "group")
        out <- nodeResult(res, n = Inf)
        return(out)
    })
    
    loc.HFDR_0.05 <- loc.HFDR_0.01 <- loc.HFDR_0.1 <- vector("list", nSIM)
    
    # edges
    chl.tree <- apply(ed, 2, FUN = function(x) {
        transNode(tree = treeR, node = x, use.alias = TRUE)})
    
    for (i in seq_len(nSIM)) {
        cat(i, "\n")
        chl.pval <- out[[i]]$PValue
        names(chl.pval) <- transNode(tree = treeR, node = out[[i]]$node, 
                                use.alias = TRUE)
        chl.hfdr_0.05 <- hFDR.adjust(chl.pval, chl.tree, alpha = 0.05)
        chl.hfdr_0.01 <- hFDR.adjust(chl.pval, chl.tree, alpha = 0.01)
        chl.hfdr_0.1 <- hFDR.adjust(chl.pval, chl.tree, alpha = 0.1)
        
        find_0.05 <- rownames(chl.hfdr_0.05@p.vals)[chl.hfdr_0.05@p.vals$adjp <= 0.05]
        find_0.05 <- find_0.05[!is.na(find_0.05)]
        find_0.05 <- transNode(tree = treeR, node = find_0.05)
        loc.HFDR_0.05[[i]] <- find_0.05[isLeaf(tree = treeR, node = find_0.05)]
        
        find_0.01 <- rownames(chl.hfdr_0.01@p.vals)[chl.hfdr_0.01@p.vals$adjp <= 0.01]
        find_0.01 <- find_0.01[!is.na(find_0.01)]
        find_0.01 <- transNode(tree = treeR, node = find_0.01)
        loc.HFDR_0.01[[i]] <- find_0.01[isLeaf(tree = treeR, node = find_0.01)]
        
        find_0.1 <- rownames(chl.hfdr_0.1@p.vals)[chl.hfdr_0.1@p.vals$adjp <= 0.1]
        find_0.1 <- find_0.1[!is.na(find_0.1)]
        find_0.1 <- transNode(tree = treeR, node = find_0.1)
        loc.HFDR_0.1[[i]] <- find_0.1[isLeaf(tree = treeR, node = find_0.1)]
        
        
    }
    
})

save(loc.HFDR_0.05, loc.HFDR_0.01, loc.HFDR_0.1, file = outRDat)


sessionInfo()