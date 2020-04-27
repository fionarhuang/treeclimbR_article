# === tool ========
.libPaths(c(Sys.getenv('R_LIBS_1'), 
            Sys.getenv('R_LIBS_2'),
            Sys.getenv('R_LIBS_3')))
.libPaths()

suppressPackageStartupMessages({
    library(TreeSummarizedExperiment)
    library(treeclimbR)
    library(dplyr)
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
(nSIM <- length(assays(tse)))
system.time({
   
    resMin <- runDA(TSE = tse, 
                     feature_on_row = TRUE, assay = 1,
                     option = "glm", group_column = "group_id",
                     design_terms = c("patient_id", "group_id"))
    outMin <-  nodeResult(object = resMin, n = Inf)
    
    keep_df <- lapply(c(0.01, 0.05, 0.1), FUN = function(x){
        getLevel(tree = rowTree(tse),
                 score_data = outMin,
                 drop = FDR > x,
                 score_column = "PValue",
                 node_column = "node",
                 get_max = FALSE,
                 parent_first = TRUE,
                 message = FALSE)$keep
    })
    outMin <- outMin %>%
        mutate(keep_0.01 = keep_df[[1]],
               keep_0.05 = keep_df[[2]],
               keep_0.1 = keep_df[[3]])

})


save(outMin, file = outRDat)

sessionInfo()