txts <- list.files("lefse/output/txt/")

txtPath <- "lefse/output/txt"
inPath <- "lefse/output/in"
resPath <- "lefse/output/res"

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
                        paste0(":", anova.alpha[j], ".res"), fi)
        resFile <- file.path(resPath, resFile)
        lefse_command <- paste("run_lefse.py", formatFile, resFile, 
                               "-a", anova.alpha[j]) 
        system(lefse_command)
        
    }
}


# txtFile <- "BS_sp50_sim5_r4_L2\\&5.txt"
# inFile <- "BS_sp50_sim5_r4_L2\\&5.in"
# resFile <- "BS_sp50_sim5_r4_L2&5.res"
# anova.alpha <- 0.05
# lda.cutoff <- 2
# 
# print(txtFile)
# print(inFile)
# 
# format_command <- paste("format_input.py", txtFile, inFile,
#                         "-c 1 -u 2 -o 1000000")
# system(format_command)
# 
# lefse_command <- paste("run_lefse.py", inFile, resFile, 
#                        "-a", anova.alpha, 
#                        "-l", lda.cutoff, 
#                        "-y", 1)
# system(lefse_command)
