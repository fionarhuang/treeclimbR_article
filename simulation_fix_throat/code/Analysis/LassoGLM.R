# === tool ========
.libPaths(c(Sys.getenv('R_LIBS_1'), 
            Sys.getenv('R_LIBS_2'),
            Sys.getenv('R_LIBS_3')))
.libPaths()

suppressPackageStartupMessages({
    library(edgeR)
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

glmnet <- function(countTab, y, isTip,
#                   normalize = FALSE,
#                   method = "TMM",
                   cv.fit = TRUE, ...){

    # define y response, site 1 as 0; site 2 as 1.
    # y <- rep(c(0,1), nSam)
    if(!class(countTab) %in% c("data.frame", "matrix")){
        stop("countTab must be a data.frame or matrix")
    }else{
        countTab <- as.matrix(countTab)
        }

    libSize <- rowsum(countTab[isTip, ], group = rep(1, sum(isTip)))[1, ]


    # create DEGList
    dy <- DGEList(countTab, remove.zeros = TRUE)
    dy$samples$lib.size <- libSize

    cd <- edgeR::cpm(dy, log = FALSE)
    mat <- t(cd)


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
fullDAT <- assays(tse, withDimnames = FALSE)
isLeaf <- rowLinks(tse)$isLeaf
y <- as.numeric(colData(tse)$group)-1


cvmod.Lasso <- lapply(seq_along(fullDAT), FUN = function(x){
    datx <- fullDAT[[x]]
    cv.mod <- glmnet(countTab = datx, y = y, isTip = isLeaf,
 #                    method = "TMM",
                     cv.fit = TRUE)
    cat(x, "out of ", length(fullDAT), "has been done", "\n")
    return(cv.mod)
})


mod.Lasso <- lapply(seq_along(fullDAT), FUN = function(x){
    datx <- fullDAT[[x]]
    mod <- glmnet(countTab = datx, y = y, isTip = isLeaf,
                  cv.fit = FALSE)
    cat(x, "out of ", length(fullDAT), "has been done", "\n")
    return(mod)
})

loc.Lasso <- lapply(seq_along(mod.Lasso), FUN = function(x){
    vv <- outLasso(cvfit = cvmod.Lasso[[x]], fit = mod.Lasso[[x]])
    cat(x, "out of ", length(mod.Lasso), "has been done", "\n")
    return(vv)
})

loc.Lasso <- lapply(loc.Lasso, FUN = function(x) {
    xx <- as.numeric(x)
    rowLinks(tse)$nodeNum[xx]
})

LS.2B <- ls()
# save result
ind <- grepl(".Lasso$",LS.2B)
obj.2B <- LS.2B[ind]
save(list=obj.2B,file=outRDat)

sessionInfo()
