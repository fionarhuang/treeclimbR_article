

suppressPackageStartupMessages({
  library(TreeHeatmap)
  library(ggnewscale)
  library(ggplot2)
  library(ggnewscale)
  library(dplyr)
  library(ggtree)
  library(cowplot)
  library(scales)
})


# ---------------------- load data ----------------------
### Edit this!!!
# resolution: low or high
reso <- "high"

# schematic scenarios: BS, US, SS
source("summary/scenario.R")

# TPR vs FDR 
source(sprintf("summary/fig_tpr_%s_reso.R", reso))

# heatmap: bottom (panel c)
s <- 5 # 5th repetition
cbPath <- file.path(sprintf("summary/figure/microbe_%s_simu.eps", reso))
source("summary/DA_heatmap.R")

fig_scene <- plot_grid(fig1 +
                         theme(plot.margin = unit(c(2, -5, 0, 5),"mm")),
                       fig2 +
                         theme(plot.margin = unit(c(2, -5, 0, 5),"mm")), 
                       fig3 +
                         theme(plot.margin = unit(c(2, -5, 0, 5),"mm")),
                       labels = c("BS", "US", "SS"), nrow = 3, 
                       hjust = -8, vjust = 3,
                       label_size = 10)
fig_up <- plot_grid(fig_scene, 
                    p_roc +
                      theme(plot.margin = unit(c(2, -10, 0, -20),"mm")), 
                    nrow = 1, 
                    labels = c("a.", "b."),
                    rel_widths = c(0.4, 1), rel_heights = c(1, 1))
fig_cb <- plot_grid(fig_up, fig_bottom, nrow = 2, 
                    labels = c(" ", "c."),
                    vjust = 0.5,
                    rel_heights = c(1.6, 1))


ggsave(cbPath, fig_cb, units = "in", width = 8, height = 8,
       dpi = 300)

