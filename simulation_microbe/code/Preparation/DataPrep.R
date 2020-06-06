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

data("throat_v35")
lib <- apply(assays(throat_v35)[[1]], 2, sum)
## ========== Select branches to swap ================
(sel1 <- selNode(obj = throat_v35, minTip = numTip1[1], 
                 maxTip = numTip1[2],
                 minPr = pr[1], maxPr = pr[2],  all = TRUE))

selR <- vector("list", nrow(sel1))
for (i in seq_len(nrow(sel1))) {
    cat(i, "\n")
    pr2 <- ratio*sel1[i, "proportion"] + c(-erp, erp)
    
    res <- try({
        sel2 <- selNode(obj = throat_v35,
                        minTip = numTip2[1], maxTip = numTip2[2],
                        minPr = pr2[1], maxPr = pr2[2],
                        all = TRUE, skip = sel1[i, "nodeNum"])},
        silent = TRUE)
    if (class(res) != "try-error") {
        selR[[i]] <- rbind.data.frame(sel1[i, ], sel2[1, ])
    } else {
        selR[[i]] <- NA
    }
}
selR
selD <- unlist(lapply(selR, FUN = function(x){!is.na(x) && all(!is.na(x$nodeNum))}))
if (sum(selD)) {
    (selR <- selR[selD][[1]])
}else{stop("No results available")}

(sp <- selR$nodeNum)


## ========== simulate data ================
if (scene %in% c("S1", "S2")) {
    pct <- NULL
}

set.seed(2020)
lse <- simData(obj = throat_v35,
               from.A = sp[1],
               from.B = sp[2],
               mu = lib,
               nSam = nSam,
               n = nSIM,
               scenario = scene,
               pct = pct)

allNode <- unique(sort(as.vector(rowTree(throat_v35)$edge)))
tse <- aggValue(x = lse, rowLevel = allNode, FUN = sum)
tse

# if the file doesn't exist or was modified earlier than the R script
if((!file.exists(outRDat)) | (file.mtime(outRDat) < file.mtime(scriptP))){
    save(tse, file = outRDat)
}



sessionInfo()








