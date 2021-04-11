# === tool ========
# .libPaths(c(Sys.getenv('R_LIBS_1'), 
#             Sys.getenv('R_LIBS_2'),
#             Sys.getenv('R_LIBS_3')))
# .libPaths()
.libPaths()

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
print(inPath)
print(txtPath)
print(resPath)
print(outFile)
print(dataPath)

load(dataPath)

nam <- names(lse)
if (!dir.exists(inPath)) {
    dir.create(inPath)
}
if (!dir.exists(resPath)) {
    dir.create(resPath)
}


# text files
txts <- list.files(txtPath)

# parameters
anova.alpha <- 0.05
lda.cutoff <- 2

lefse <- vector("list", length(lse))
names(lefse) <- nam

for (i in seq_along(nam)) {
    message(i, " out of ", length(lse), " have been done.")
    lefse[[i]] <- system.time({
        fi <- txts[grep(pattern = nam[i], x = txts)]
        
        txtFile <- file.path(txtPath, fi)
        formatFile <- file.path(inPath, 
                                gsub(pattern = ".txt", ".in", fi))
        
        # step 1: format input
        format_command <- paste("format_input.py", txtFile, 
                                formatFile,
                                "-c 1 -u 2 -o 1000000")
        system(format_command)
        
        # step 2: run lefse analysis
        resFile <- gsub(pattern = "\\.txt", ".res", fi)
        resFile <- file.path(resPath, resFile)
        lefse_command <- paste("run_lefse.py", formatFile, resFile, 
                               "-a", anova.alpha) 
        system(lefse_command)
    })
}

save(lefse, file = outFile)


sessionInfo()