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
dataFile <- "../data/Weber_AML_sim_main_0.1pc_SE.rds"
    
# the number of cluster
nx <- 30
ny <- 30

# the disease type
condition <- c("healthy", "CN")

# the seed to do cluster
type <- "DA"
seedN <- 2020

# define clusters with above 50% spike-in cells as true spike-in cluster
thr_spikein <- 0.5

scene <- "CN, 0.1%"

ls.after <- ls()
ls.save <- setdiff(ls.after, c(ls.before, "ls.before"))

# if the file doesn't exist or was modified earlier than the R script

if((!file.exists(outRDat)) | (file.mtime(outRDat) < file.mtime(scriptP))){
    save(list = ls.save, file = outRDat)
}