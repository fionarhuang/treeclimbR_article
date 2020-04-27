# === tool ========
source("/home/fiona/phd/microbes/simulation/Tool/argsR_compile.R")


library(treeAGG2)
library(ggtree)
library(ggplot2)
library(cowplot)
R.Version()

# ==== arguments from batch R===================
args <- (commandArgs(trailingOnly = TRUE))
args
argsList <- argRun(args, grp.pattern = "dirGRP")
argsList
for (i in seq_along(argsList)) {
    if(length(argsList[[i]])>1){
        assign(names(argsList)[i],argsList[[i]])
    }else{eval(parse(text = argsList[[i]]))}
}

# ==========read in parameter value================#
load(inRDat)

nn <- length(assays(tse))
n <- ncol(tse)
nNode <- 2*rowTree(tse)$Nnode + 1

for (i in seq_len(nn)){
vl <- as.vector(assays(tse)[[i]])
dat <- data.frame(node = rep(seq_len(nNode), n),
                  sample = rep(seq_len(n), each = nNode),
                  group = rep(1:2, each = nNode*n/2),
                  value = vl)
otf <- gsub(pattern = "Fig.txt", 
            replacement = paste0("Fig_", i,".png"), out)

png(otf, height = 600, width = 600, bg = NA)
fig <- ggplot(dat[dat$node %in% seq_len(nNode), ]) +
    geom_boxplot(aes(x = factor(group), y = value, group = group)) +
    geom_jitter(aes(x = group, y = value), size = 0.2,
                alpha = 0.4, color = "blue", width = 0.1, height = 0) +
    facet_wrap(~ node, scales = "free") +
    theme(strip.text = element_text(size = 18),
          panel.background = element_rect(fill = "transparent", color = NA))
print(fig)
dev.off()

capture.output(writeLines("success"), file = out)
}