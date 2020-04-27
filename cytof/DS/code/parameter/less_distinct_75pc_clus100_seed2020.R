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

# the data path
dataFile <- "../data/Weber_BCR_XL_sim_less_distinct_less_75pc_SE.rds"
    
# the number of cluster
nx <- 10
ny <- 10

nclus <- 100

# the seed to do cluster
seedN <- 2020



ls.after <- ls()
ls.save <- setdiff(ls.after, c(ls.before, "ls.before"))

# if the file doesn't exist or was modified earlier than the R script

if((!file.exists(outRDat)) | (file.mtime(outRDat) < file.mtime(scriptP))){
    save(list = ls.save, file = outRDat)
}