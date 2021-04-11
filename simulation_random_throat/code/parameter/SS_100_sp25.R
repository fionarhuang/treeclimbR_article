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
nSIM <- 100
scene <- "SS"
# seed: we use the same R file (code/Preparation/DataPrep.R) to simulate data for three scenarios: BS, US, SS
# if randomly selected branches needs to be different among three scenarios,
# then different seedNum should be set.
seedNum <- 1234

if (scene %in% c("BS", "SS")) {
    pct = 1
} else {
    pct = 0.6
}

ls.after <- ls()
ls.save <- setdiff(ls.after, c(ls.before, "ls.before"))
# if the file doesn't exist or was modified earlier than the R script
if((!file.exists(outRDat)) | (file.mtime(outRDat) < file.mtime(scriptP))){
  save(list = ls.save, file = outRDat)
}
