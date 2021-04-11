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

dir.create(inPath)
dir.create(resPath)

# text files
txts <- list.files(txtPath)

# parameters
anova.alpha <- c(0.01, 0.05, 0.1)
lda.cutoff <- 2



for (i in seq_along(txts)) {
    fi <- txts[i]
    
    txtFile <- file.path(txtPath, fi)
    formatFile <- file.path(inPath, 
                            gsub(pattern = ".txt", ".in", fi))
    
    # step 1: format input
    format_command <- paste("format_input.py", txtFile, 
                            formatFile,
                            "-c 1 -u 2 -o 1000000")
    system(format_command)
   
    # step 2: run lefse analysis
     for (j in seq_along(anova.alpha)) {
         message(i, "-", j)
        resFile <- gsub(pattern = "\\.txt",
                        paste0("@", anova.alpha[j], ".res"), fi)
        resFile <- file.path(resPath, resFile)
        lefse_command <- paste("run_lefse.py", formatFile, resFile, 
                               "-a", anova.alpha[j]) 
        system(lefse_command)
        
    }
}


sink(outFile)
cat("successfully run lefse... ")
sink()