# === tool ========
.libPaths(c(Sys.getenv('R_LIBS_1'), 
            Sys.getenv('R_LIBS_2'),
            Sys.getenv('R_LIBS_3')))
.libPaths()

suppressPackageStartupMessages({
    library(TreeSummarizedExperiment)
    library(glmnet)
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

## === load data from data preparation step ===
load(inRDat)

## ================ functions =================

glmnet <- function(countTab, y, 
                   cv.fit = TRUE, ...){

    # define y response, base as 0; spike as 1.
    # y <- rep(c(0,1), nSam)
    if(!class(countTab) %in% c("data.frame", "matrix")){
        stop("countTab must be a data.frame or matrix")
    }else{
        countTab <- as.matrix(countTab)
        }

    mat <- t(countTab)

    # lasso regularisation
    if(cv.fit){
        cv.mod <- glmnet::cv.glmnet(x = mat, y = y,
                                    family="binomial",...)
        return(cv.mod)
    }else{
        mod <- glmnet::glmnet(x = mat, y = y, family="binomial",
                              ...)
        return(mod)
    }
}

outLasso <- function(cvfit, fit){

    lambda <- cvfit$lambda.min
    coefs <- glmnet::coef.glmnet(fit, s = lambda)

    vSig <- rownames(coefs)[coefs[,1] != 0]
    vSig1 <- setdiff(vSig, "(Intercept)")
    return(vSig1)

}

# -------------- Data analysis -----------------
is_state <- metadata(d_medians_all)$id_state_markers
assayDAT <- assays(d_medians_all, withDimnames = FALSE)[is_state]
fullDAT <- lapply(seq_along(assayDAT), FUN = function(x){
    xx <- assayDAT[[x]]
    rownames(xx) <- paste0(rownames(xx), ".", names(assayDAT)[x])
    return(xx)
})
fullDAT <- do.call(rbind, fullDAT)

y <- as.numeric(colData(d_medians_all)$group_id)-1
notNa <- apply(fullDAT, 1, FUN = function(x) {
        all(!is.na(x))})
    
cvmod.Lasso <- glmnet(countTab = fullDAT[notNa, ], y = y,
                     cv.fit = TRUE)
    
mod.Lasso  <- glmnet(countTab = fullDAT[notNa, ], y = y,
                  cv.fit = FALSE)
   

node_find <-  outLasso(cvfit = cvmod.Lasso, fit = mod.Lasso)
loc.Lasso_ab <- split(x = sub(pattern = "\\..*", replacement = "", node_find), 
                       f = gsub(pattern = ".*\\.", replacement = "", node_find))
loc.Lasso_ab <- lapply(loc.Lasso_ab, FUN = function(x){
    transNode(tree = rowTree(d_medians_all), 
                       node = x, use.alias = FALSE)})
loc.Lasso <- unique(unlist(loc.Lasso_ab))


# save result
save(loc.Lasso_ab, loc.Lasso, file = outRDat)
