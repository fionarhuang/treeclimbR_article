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
rejLimit <- 0.05

system.time({
    # calculate p values
    res_agg <- runDA(TSE = tse, 
                    feature_on_row = TRUE, assay = 1,
                    option = "glm", group_column = "group_id",
                    design_terms = c("patient_id", "group_id"))
    out_agg <-  nodeResult(object = res_agg, n = Inf)
    
    
    # p values
    chl.pval <- out_agg$PValue
    names(chl.pval) <- out_agg$nodeLab
    
    # edges
    treeR <- rowTree(tse)
    ed <- treeR$edge
    chl.tree <- apply(ed, 2, FUN = function(x) {
        transNode(tree = treeR, node = x, use.alias = FALSE)
    })
    
    # leaf nodes
    leaf <- setdiff(chl.tree[, 2], chl.tree[, 1])
    
    # found
    rej_limit <- c(0.01, 0.05, 0.1)
    chl.hfdr <- loc.hc <- vector("list", 3)
    names(chl.hfdr) <- names(loc.hc) <- rej_limit
    
    for (i in seq_along(rej_limit)) {
        hi <- hFDR.adjust(chl.pval, chl.tree, 
                          alpha = rej_limit[i])
        chl.hfdr[[i]] <- hi
        found <- rownames(hi@p.vals)[hi@p.vals$adjp <= rej_limit[i]]
        found <- found[!is.na(found)]
        
        # leaf level
        loc.hc[[i]] <- intersect(found, leaf)
    }
   
    
})

save(loc.hc, file = outRDat)
sessionInfo()