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
print(inRDat)
load(inRDat)

# ------------------------------- data ---------------------------------
X <- lapply(tse, FUN = function(x) {
    assays(x)[[1]]
})
Y <- lapply(tse, FUN = function(x) {
    gr <- as.factor(colData(x)$group)
    as.numeric(gr) - 1 
})
treeX <- lapply(tse, rowTree)

## ================ functions =================

glmnet <- function(countTab, y, cv.fit = TRUE, ...){

    # define y response, site 1 as 0; site 2 as 1.
    countTab <- as.matrix(countTab)
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

Lasso <- vector("list", length(tse))
names(Lasso) <- names(tse)


for (i in seq_along(tse)) {
    Lasso[[i]] <- system.time({
        cv.mod <- glmnet(countTab = X[[i]], y = Y[[i]], cv.fit = TRUE)
        mod <- glmnet(countTab = X[[i]], y = Y[[i]], cv.fit = FALSE)
        vv <- outLasso(cvfit = cv.mod, fit = mod)
        transNode(tree = treeX[[i]], node = vv)
    })
    
}


Lasso_all <- mapply(FUN = function(x, y) {
    x + y
}, Time_agg, Lasso, SIMPLIFY = FALSE)

# print out the run time
Lasso

Lasso_all

save(Lasso, Lasso_all, file=outRDat)

sessionInfo()
