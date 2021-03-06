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
load(inRDat)

## === differential abundance test using edgeR ===
(nSIM <- length(tse_list))
test.func <- function(X, Y) {
    runP <- function(X, Y,
                     normalize = TRUE,
                     adjust.method = "BH", 
                     prior.count = 0.125) {
        
        y <- DGEList(X, remove.zeros = FALSE)
        design <- model.matrix(~Y)
        
        
        # do normalisation
        if (normalize) {
            y <- calcNormFactors(y, method = "TMM")
        }
        
        # estimate dispersion
        y <- estimateDisp(y, design = design)
        
        # build model
        fit <- glmFit(y, design = design, prior.count = prior.count)
        
        # extract results
        lrt <- glmLRT(fit, contrast = NULL)
        tt1 <- topTags(lrt, n = Inf, adjust.method = adjust.method,
                       sort.by = "none")$table
        tt2 <- tt1[rownames(X), ]
        
        
        return(tt2)
    }
    out <- runP(X, Y)
    ro <- out %>% 
        data.frame %>%
        select(PValue, logFC) %>%
        mutate(sign = ifelse(logFC < 0, -1, 1))
    pv <- ro$PValue
    es <- ro$sign
    names(pv) <- names(es) <- rownames(ro)
    list(p.value = pv, e.sign = es)
}


perm.func <- function (X, Y) {  
    return(list(X=X, Y=sample(Y))) 
}


system.time({
    lse_list <- lapply(tse_list, FUN = function(x) {
        x[rowLinks(x)$isLeaf, ]})
    
    loc.StructFDR_0.05 <- loc.StructFDR_0.01 <- loc.StructFDR_0.1 <- vector("list", nSIM)
    for (i in seq_len(nSIM)) {
        cat(i, "\n")
        X <- assays(lse_list[[i]])[[1]]
        rownames(X) <- rowLinks(lse_list[[i]])$nodeLab
        Y <- as.numeric(colData(lse_list[[i]])$group) - 1
        
        tree.fdr.obj <- TreeFDR(X, Y, rowTree(lse_list[[i]]), test.func, perm.func) 
        
        # 0.05
        sel <- tree.fdr.obj$p.adj <= 0.05
        sx <- rownames(X)[sel]
        loc.StructFDR_0.05[[i]] <- transNode(tree = rowTree(lse_list[[i]]), 
                                             node = as.character(sx))
        
        # 0.01
        sel <- tree.fdr.obj$p.adj <= 0.01
        sx <- rownames(X)[sel]
        loc.StructFDR_0.01[[i]] <- transNode(tree = rowTree(lse_list[[i]]), 
                                       node = as.character(sx))
        
        # 0.1
        sel <- tree.fdr.obj$p.adj <= 0.1
        sx <- rownames(X)[sel]
        loc.StructFDR_0.1[[i]] <- transNode(tree = rowTree(lse_list[[i]]), 
                                       node = as.character(sx))
    }
    
})


save(loc.StructFDR_0.01, loc.StructFDR_0.05, loc.StructFDR_0.1, file = outRDat)

sessionInfo()
