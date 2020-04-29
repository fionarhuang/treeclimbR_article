# === tool ========
R.Version()

# =================== arguments from batch R=============
argsList <- (commandArgs(trailingOnly = TRUE))
argsList
for (i in seq_along(argsList)) {
    if(length(argsList[[i]])>1){
        assign(names(argsList)[i],argsList[[i]])
    }else{eval(parse(text = argsList[[i]]))}
}

# ======================== parameters ===================
ls.before <- ls()


nSam <- c(25, 25)
nSIM <- 5
limFDR <- 0.05
ratio <- 4
scene <- "BS"
pr <- c(0.01, 0.05)
numTip1 <- c(11, 15)
erp <- 0.01
numTip2 <- c(0, 40)

ls.after <- ls()
ls.save <- setdiff(ls.after, c(ls.before, "ls.before"))
# if the file doesn't exist or was modified earlier than the R script
if((!file.exists(outRDat)) | (file.mtime(outRDat) < file.mtime(scriptP))){
    save(list = ls.save, file = outRDat)
}