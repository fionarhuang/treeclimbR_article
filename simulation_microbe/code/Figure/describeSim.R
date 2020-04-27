
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
(pd <- data.frame(value = c(pr[, "p1"], pr[, "p2"]), 
                  Group = rep(LETTERS[1:2], each = 10),
                  Entity = rep(1:10, 2),
                  Different = rep(rep(c(TRUE, FALSE), each = 5), 2)))

png(outP, height = 300, width = 600, bg = NA)
f1 <- ggtree(tinyTree, size = 2)+
    geom_tiplab(aes(label = node), color = "darkblue",
                hjust = -0.5, vjust = 0.7, size = 6) +
    geom_hilight(node = 18, fill = "blue", alpha = 0.3) +
    geom_hilight(node = 13, fill = "orange", alpha = 0.3) +
    xlim(c(0, 3))
f2 <- ggplot(pd) +
    geom_point(aes(x = Entity, y = value,
                   color = Group, shape = Group),
               size = 4) +
    scale_x_discrete(name = "Entities", limits = as.character(1:10)) +
    scale_y_continuous(name = "Proportion",
                       limits = c(0, 0.20)) +
    theme(legend.position = "top",
          axis.text = element_text(size = 18),
          axis.title = element_text(size = 20, face = "bold"))
plot_grid(f1, f2, rel_widths = c(1,1.5))
dev.off()


