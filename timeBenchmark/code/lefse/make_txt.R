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
print(dataPath)
print(txtPath)
print(outDat)

if (!dir.exists(txtPath)) {
    dir.create(txtPath)  
}


# load data
load(dataPath)

# trees
treeX <- lapply(lse, rowTree)


# the time to create text files
txtTime <- vector("list", length(lse))
names(txtTime) <- names(lse)

for (i in seq_along(lse)) {
    txtTime[[i]] <- system.time({
        # hierarchy
        tree <- rowTree(lse[[i]])
        path <- matTree(tree = tree)
        pathL <- lapply(seq_len(nrow(path)), FUN = function(x) {
            xx <- path[x, ]
            xx <- xx[!is.na(xx)]
            tx <- transNode(tree = tree, node = xx, use.alias = TRUE)
            paste(rev(tx),collapse = "|")
        })
        
        # tabular data required by lefse
        oo <- transNode(tree = tree, node = path[, 1], use.alias = FALSE)
        count <- assays(lse[[i]])[[1]]
        class <- colData(lse[[i]])$group
        
        count_i <- count[oo, ]
        rownames(count_i) <- unlist(pathL)
        
        df <- rbind.data.frame(
            class,
            colnames(count_i),
            count_i)
        rownames(df) <- c("class", "sample_id", rownames(count_i))
        
        # output
        file_i <- file.path(txtPath, paste0(names(lse)[i], ".txt"))
        message(file_i)
        write.table(x = df, file = file_i, sep = "\t",
                    row.names = TRUE, col.names = FALSE, quote = FALSE)
    })
}


save(treeX, txtTime, file = outDat)

