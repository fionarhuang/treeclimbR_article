# === tool ========
.libPaths(c(Sys.getenv('R_LIBS_1'), 
            Sys.getenv('R_LIBS_2'),
            Sys.getenv('R_LIBS_3')))
.libPaths()

suppressPackageStartupMessages({
    library(ape)
    library(dplyr)
    library(TreeSummarizedExperiment)
    library(stringr)
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

## ===============================================
print(dataFile)
print(txtFile)
print(txtPath)
print(outTree)


if (!file.exists(txtPath)) {
    dir.create(txtPath)
}

ss <- gsub(pattern = ".*/", "", dataFile)
scene <- gsub("[.]RData", "", ss)


outFile <- file.path(txtPath, paste0(scene, ".txt"))
 

load(dataFile)
# ------------------------------- data format ---------------------------------
# tree is the same in 100 repetitions
tree <- rowTree(tse)
save(tree, file = outTree)

path <- matTree(tree = tree)
pathL <- lapply(seq_len(nrow(path)), FUN = function(x) {
    xx <- path[x, ]
    xx <- xx[!is.na(xx)]
    tx <- transNode(tree = tree, node = xx, use.alias = TRUE)
    paste(rev(tx), collapse = "|")
})

# tabular data required by lefse
oo <- transNode(tree = tree, node = path[, 1], use.alias = TRUE)
lse <- tse[rowLinks(tse)$isLeaf, ]


# repetitions
count <- assays(lse)
class <- colData(lse)$group
for (j in seq_along(count)) {
    count_j <- count[[j]][oo, ]
    rownames(count_j) <- unlist(pathL)
    
    df <- rbind.data.frame(
        class,
        colnames(count_j),
        count_j)
    rownames(df) <- c("class", "sample_id", rownames(count_j))
    
    # output
    file_j <- gsub(pattern = "\\.txt", paste0("_rep", j, ".txt"), outFile)
    message(file_j)
    write.table(x = df, file = file_j, sep = "\t",
                row.names = TRUE, col.names = FALSE, quote = FALSE)
}


sink(txtFile)
cat("successfully export txt files... ")
sink()

sessionInfo()