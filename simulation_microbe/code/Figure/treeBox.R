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
    # file
    otf <- gsub(pattern = "Box.txt", 
                replacement = paste0("Box_", i,".png"), out)
    
    vl <- as.vector(assays(tse)[[i]])
    dat <- data.frame(node = rep(seq_len(nNode), n),
                      sample = rep(seq_len(n), each = nNode),
                      group = rep(1:2, each = nNode*n*0.5),
                      value = vl)
    px <- vector("list", nNode)
    for (i in seq_len(nNode)) {
        px[[i]] <- ggplot(dat[dat$node %in% i, ]) +
            geom_boxplot(aes(x = as.factor(group), y = value, group = group)) +
            geom_jitter(aes(x = group, y = value), size = 0.2,
                        alpha = 0.8, color = "blue", width = 0.15, height = 0) +
      #      facet_wrap(~ node, scales = "free") +
            facet_wrap(~ node) +
            theme(strip.text = element_text(size = 15),
                  axis.text.x = element_blank(),
                  panel.background = element_rect(fill = "white")) +
            xlab("") +
            ylab("")
    }
    
    
    png(otf, height = 650, width = 800, bg = NA)
    treeP <- ggtree(rowTree(tse), branch.length = "none",
                    color = "grey20")+ xlim(c(0, 8))+
        geom_hilight(node = 18, fill = "blue", alpha = 0.4) +
        geom_hilight(node = 13, fill = "orange", alpha = 0.4)
    
    names(px) <- seq_len(nNode)
    p <- inset(treeP, px[1:10], width = 1.2, height = 1.6, vjust = 0.2)
    pp <- inset(p, px[11:19], width = 1.4, height = 2.2, vjust = 0.2)
    print(pp)
    dev.off()
    
    capture.output(writeLines("success"), file = out)
}