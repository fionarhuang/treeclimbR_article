suppressPackageStartupMessages({
    library(cowplot)
    library(ggplot2)
    library(TreeHeatmap)
    library(ggnewscale)
    library(treeclimbR)
    library(scales)
    library(dplyr)
})


#---------------------------------------------------------------------------
# TPR vs FDR
#---------------------------------------------------------------------------
# work directory is in "cytof/
source("DA/code/Figure/DA_roc.R")
source("DS/code/Figure/DS_roc.R")

## ROC
vshape <- c("lasso" = 16, "minP" = 16, "StructFDR" = 16, "diffcyt" = 16, 
            "treeclimbR" = 5, "HFDR" = 16)

vcolor <- c("treeclimbR" = "#E41A1C", "diffcyt" = "#377EB8",
            "StructFDR" = "#4DAF4A", "HFDR" = "#984EA3",
            "lasso" = "#FF7F00", "minP" = "#A65628")
vsize <- c("treeclimbR" = 2.5, "diffcyt" = 2,
           "StructFDR" = 2, "HFDR" = 2,
           "lasso" = 2, "minP" = 2)
prettify <- theme_bw(base_size = 8) + theme(
    aspect.ratio = 1,
    #plot.margin = unit(c(0, 0, 0.2, 0), "cm"),
    panel.grid = element_blank(),
    panel.spacing = unit(0, "lines"),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.key.size= unit(1.5, "mm"),
    plot.title = element_text(hjust = 0.5),
    legend.text = element_text(size = 7),
    legend.position="bottom",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin=margin(-5,-5,-10,-5),
    strip.background = element_rect(colour = "black", fill = "gray90"),
    strip.text.x = element_text(color = "black", size = 8),
    strip.text.y = element_text(color = "black", size = 8))

set.seed(2)
rS <- runif(min = -0.02, max = 0.02, n = nrow(Dat_DS))
p_DS <- ggplot(Dat_DS) +
    geom_vline(xintercept = c(0.01, 0.05, 0.1),
               linetype = "longdash", color = "black",
               size = 0.2) +
    geom_point(aes(x = fdr, y = tpr + rS, color = method, 
                   shape = method, size = method),
               stroke = 0.8) +
    geom_path(aes(x = fdr, y = tpr + rS,
                  color = method, group = method), 
              size = 0.3, show.legend = FALSE) +
    #    facet_wrap(~output) +
    xlab("FDR") +
    ylab("TPR") +
    facet_grid(rows = vars(cluster), cols = vars(pc)) +
    scale_color_manual("Method", values = vcolor) +
    scale_shape_manual("Method", values = vshape) +
    scale_size_manual("Method", values = vsize) +
    guides(color=guide_legend(nrow = 1, byrow = TRUE,
                              override.aes = list(size = 2))) +
    scale_x_continuous(trans = "sqrt", 
                       breaks = c(0.01, 0.05, 0.1, seq(0, 1, 0.2)),
                       limits = c(0, 1), labels = function(u) formatC(u)) +
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    prettify 
p_DS
set.seed(2)
rA <- runif(min = -0.02, max = 0.02, n = nrow(Dat_DA))

p_DA <- ggplot(Dat_DA) +
    geom_vline(xintercept = c(0.01, 0.05, 0.1),
               linetype = "longdash", color = "black",
               size = 0.2) +
    geom_point(aes(x = fdr, y = tpr + rA, color = method, 
                   shape = method, size = method),
               stroke = 0.8) +
    geom_path(aes(x = fdr, y = tpr + rA,
                  color = method, group = method), 
              size = 0.3, show.legend = FALSE) +
    xlab("FDR") +
    ylab("TPR") +
    facet_grid(rows = vars(cluster), cols = vars(pc)) +
    scale_color_manual("Method", values = vcolor) +
    scale_shape_manual("Method", values = vshape) +
    scale_size_manual("Method", values = vsize) +
    guides(color=guide_legend(nrow = 2, byrow = TRUE,
                              override.aes = list(size = 2))) +
    scale_x_continuous(trans = "sqrt", 
                       breaks = c(0.01, 0.05, 0.1, seq(0, 1, 0.2)),
                       limits = c(0, 1), labels = function(u) formatC(u)) +
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    prettify 
p_DA
seed <- 2020

roc_both <- plot_grid(
    p_DA + theme(
        plot.margin = unit(c(3.5, -5, 3.5, -5),"mm"),
        legend.position = "none") + 
        labs(title = "AML-sim"),
    p_DS + theme(plot.margin = unit(c(3.5, -5, 3.5, -5),"mm"),
                 legend.position = "none") + 
        labs(title = "BCR-XL-sim"),
    labels = c("a", "c"),
    hjust = c(-2, -0.3),
    ncol = 2
)
legend <- get_legend(p_DA + 
                         theme(legend.key.size= unit(3, "mm"),
                               legend.text = element_text(size = 8),
                               legend.box.margin = margin(t = 0, 0, 20, 0)))





source("DA/code/Figure/DA_heatmap_pie.R")
hm_DA

space <- ggdraw()
hm_DA_legend <- plot_grid(legend, 
                          hm_DA +
                              labs(title = "AML-sim, medium, 900") + 
                              theme(plot.margin = unit(c(-4, 7, 0, 0),"mm"),
                                    plot.title = element_text(hjust= 0.5, 
                                                              vjust = -2,
                                                              face = "plain",
                                                              size = 10),
                                    legend.box.margin = 
                                        margin(t = 3, b = 5, 
                                               r = -5, l = 2),
                                    legend.text = element_text(size = 5),
                                    legend.position = c(1, 0.55),
                                    legend.spacing.y = unit(0.5, "mm"),
                                    legend.key.size= unit(2, "mm"),
                                    legend.title = element_text(size = 6)),
                          space, 
                          ncol = 1,
                          rel_heights = c(0.08, 0.3, 0.04),
                          labels = c("", "b"),
                          hjust = c(0, -2.2),
                          vjust = c(1.5, 0))
hm_DA_legend 
source("DS/code/Figure/DS_heatmap_pie.R")
hm_DS
hm_both <- plot_grid(hm_DA_legend,
                     hm_DS +
                         labs(title = "BCR-XL-sim, medium, 400") + 
                         theme(plot.margin = unit(c(-4, 3.5, 0, 0),"mm"),
                               plot.title = element_text(hjust= 0.5, 
                                                         vjust = -5,
                                                         face = "plain",
                                                         size = 10),
                               legend.box.margin = 
                                   margin(t = 3, b = 5, 
                                          r = -10, l = -25),
                               legend.text = element_text(size = 5),
                               legend.position = "right",
                               legend.spacing.y = unit(0.5, "mm"),
                               legend.key.size= unit(2, "mm"),
                               legend.title = element_text(size = 6)), 
                     labels = c("", "d"), 
                     rel_widths = c(0.5, 0.5),
                     hjust = c(-2, -2.3),
                     vjust = c(2, 0.6),
                     ncol = 2)
hm_both

cb_both <- plot_grid(roc_both, hm_both, nrow = 2, rel_heights = c(0.55, 0.45))
figPath <- file.path(sprintf("summary/cytof_%d.eps", seed))
ggsave(figPath, cb_both, units = "in", width = 8, height = 8,
       dpi = 300)


# run below code to get fig_list
source("DS/code/Figure/DS_heatmap_pie.R")

load("DS_resolution/AUC_reso.RData")
df_auc <- data.frame(n_cluster = reso^2, 
                     diffcyt = res_auc$diffcyt, 
                     treeclimbR = res_auc$treeclimbR) %>%
    gather(Method, AUC, -n_cluster)

fig_auc <- ggplot(df_auc, aes(x = n_cluster, y = AUC, color = Method)) +
    geom_point() +
    geom_line() +
    scale_color_manual(values = c("treeclimbR" = "#E41A1C", "diffcyt" = "#377EB8")) +
    scale_x_sqrt(labels = reso^2, breaks = reso^2) +
    ylim(c(0, 1)) + 
    labs(x = "The number of clusters", title = "BCR-XL-sim, medium") +
    theme_bw(base_size = 8) + theme(
        aspect.ratio = 1,
        panel.grid = element_blank(),
        panel.spacing = unit(0, "lines"),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size= unit(2.5, "mm"),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 8),
        legend.position=c(0.5, 0.5),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin=margin(-5,-5,-10,-5),
        strip.background = element_rect(colour = "black", fill = "gray90"),
        strip.text.x = element_text(color = "black", size = 8),
        strip.text.y = element_text(color = "black", size = 8))

fig_list <- c(fig_list, list(fig_auc))
fig_more <- plot_grid(plotlist = fig_list, nrow = 2, 
                      rel_heights = c(0.5, 0.5))
figPath <- "summary/Supplementary_cytof_2020_more.eps"
ggsave(figPath, fig_more, units = "in", width = 8, height = 8.5,
       dpi = 300)
