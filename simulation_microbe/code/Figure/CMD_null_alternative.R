library(ape)
library(TreeSummarizedExperiment)
library(treeAGG2)
library(dplyr)

# example from S2, 50 replicates
mdir <- "/Volumes/fiona/phd/microbes/simulation/wflow_microbes/output/RData"
si <- "S2_lib1E5_sp50_sim5_r4_L2"
load(file.path(mdir, "Analysis", si, "edgeR_aggP.RData"))
load(file.path(mdir, "DataPrep", paste0(si, ".RData")))

# tree
tr <- rowTree(tse)

# signal branch
br <- metadata(tse)$br

# branches without signal
fc <- metadata(tse)$FC
no_leaf <- names(fc)[fc == 1]
no_leaf <- transNode(tree = tr, node = no_leaf, use.alias = FALSE)
no_br <- signalNode(tree = tr, node = no_leaf)
no_node <- findOS(tree = tr, node = no_br, 
                  only.leaf = FALSE, self.include = TRUE)

# branches with signal
# yes_leaf <- names(fc)[fc != 1]
# yes_br <- signalNode(tree = tr, node = yes_leaf)
# yes_node <- findOS(tree = tr, node = yes_br, 
#                 only.leaf = FALSE, self.include = TRUE)
yes_leaf <- findOS(tree = tr, node = c(br$A, br$B), only.leaf = TRUE)
yes_br <- c(br$A, br$B)
yes_node <- findOS(tree = tr, node = c(br$A, br$B), only.leaf = FALSE,
                   self.include = TRUE)
true_leaf <- findOS(tree = tr, node = c(br$A, br$B), only.leaf = TRUE)
true_br <- c(br$A, br$B)
true_node <- findOS(tree = tr, node = c(br$A, br$B), only.leaf = FALSE,
                    self.include = TRUE)

# the first simulation
ff <- function(x) {
    xx <- findOS(tree = tr, node = x, only.leaf = TRUE, self.include = TRUE)
    unlist(lapply(xx, length))
}
out <- outsel_0.05[[1]]
out <- out %>%
    mutate(len = ff(node))

out %>%
    mutate(PV = round(PValue, 4)) %>%
    filter(node %in% unlist(yes_node)) %>%
    select(node, PV, PValue, logFC)

lev <- getCand(tree = rowTree(tse), 
               score_data = out,
               node_column = "node", 
               p_column = "PValue",
               sign_column = "logFC",
               message = TRUE)

outF <- getBest(tree = rowTree(tse), 
               levels = lev$candidate_list, 
               score_data = out, 
               p_column = "PValue",
               sign_column = "logFC",
               method = "BH", 
               node_column = "node", 
               limit_rej = 0.05)



names(true_leaf) <- LETTERS[1:2]
# viewBranch(tree = tr,  
#            hlight_node = c(br$A, br$B), 
#            hlight_fill = c("orange", "blue"), 
#            point_node = sel_i, point_color = "darkviolet", point_size = 1,
#            zoom_node = br$A, zoom_scale = 40, 
#            group_leaf = c(true_leaf, list("0" = no_leaf)), 
#            group_color = c("A" = "orange", "B" = "blue", "0" = "grey"), 
#            edge_size = 0.5)

# level i
tt <- outF$level_info$T
df <- vector("list", length(tt))
for (j in seq_along(tt)) {
    sel_i <- lev$candidate_list[[j]]
    #yes_i <- intersect(unlist(yes_node), sel_i)
    #yes_i <- setdiff(sel_i, unlist(no_node))
    yes_i <- intersect(unlist(yes_node), sel_i)
    no_i <- intersect(unlist(no_node), sel_i)
    
    
    
    df_no <- out %>%
        data.frame() %>%
        dplyr::filter(node %in% sel_i) %>%
        mutate(adj.p = p.adjust(PValue, method = "BH")) %>%
        dplyr::filter(node %in% no_i) %>%
        dplyr::mutate(t = tt[j]) %>%
        dplyr::mutate(truth = "NO") %>%
        dplyr::mutate(tP = ifelse(logFC > 0, PValue, -PValue)) %>%
        dplyr::select(PValue, adj.p, t, truth, node, len, tP) 
    
    df_yes <- out %>%
        data.frame() %>%
        dplyr::filter(node %in% sel_i) %>%
        mutate(adj.p = p.adjust(PValue, method = "BH")) %>%
        dplyr::filter(node %in% yes_i) %>%
        dplyr::mutate(t = tt[j]) %>%
        dplyr::mutate(truth = "YES") %>%
        dplyr::mutate(tP = ifelse(logFC > 0, PValue, -PValue)) %>%
        dplyr::select(PValue, adj.p, t, truth, node, len, tP) 
    
    df[[j]] <- rbind(df_no, df_yes)
    
    # hist(df_no, freq = FALSE, main = outF$threshold$T[j], breaks = 40)
    
}

# leaf level
df_leaf <- out %>%
    data.frame() %>%
    mutate(isLeaf = isLeaf(tree = tr, node = node)) %>%
    dplyr::filter(isLeaf) %>%
    mutate(adj.p = p.adjust(PValue, method = "BH")) %>%
    mutate(len = ff(node)) %>%
    dplyr::filter(node %in% c(unlist(yes_leaf), no_leaf)) %>%
    dplyr::mutate(t = 0) %>%
    dplyr::mutate(truth = ifelse(node %in% unlist(yes_leaf), "YES", "NO")) %>%
    dplyr::mutate(tP = ifelse(logFC > 0, PValue, -PValue)) %>%
    dplyr::select(PValue, adj.p, t, truth, node, len, tP) 


df <- c(list(df_leaf), df)



# t : distinct result
df_yes <- lapply(df, FUN = function(x) {
    x  %>% 
        filter(truth == "YES") %>% select(-t)
})
sel_yes <- lapply(df_yes, FUN = function(x) {
    xl <- lapply(df_yes, FUN = function(y) {
    all(all.equal(y, x) == TRUE)
    })
    unlist(xl)
    })
sel_yes <- lapply(unique(sel_yes), FUN = function(x) {which(x)[1]})
sel_yes <- unlist(sel_yes)

# H1:
df_yes <- lapply(df[sel_yes], FUN = function(x) {
    x  %>% filter(truth == "YES")
})
df_yes <- do.call(rbind, df_yes)

df_no <- lapply(df, FUN = function(x) {
    x %>% filter(truth == "NO")
})
df_no <- do.call(rbind, df_no)



library(viridis)
# ggplot() +
#     stat_ecdf(data = df[df$truth == "NO", ],
#               aes(PValue, group = t, color = t), 
#               geom = "step") +
#     geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
#     scale_color_manual(values = brewer.pal(21, "RdBu")[-c(6:15)])
# #    scale_color_viridis(discrete = TRUE)
# mcol <- colorRampPalette(c("orange", "red", "violet", "blue"))( 21 )
# f <- function(x, y){
#     ff <- ecdf(x)
#     ff(y)
# }
df_yes<- df_yes %>%
    slice(rep(1:n(), len))
 # %>% group_by(t) %>%
 #    mutate(label = f(adj.p, t*0.05))


prettify <- theme_bw(base_size = 8) + theme(
    rect = element_rect(fill = "transparent"),
    aspect.ratio = 1,
    #plot.margin = unit(c(0, 0, 0.2, 0), "cm"),
    panel.grid = element_blank(),
    panel.spacing = unit(0, "lines"),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.key.size= unit(2, "mm"),
    plot.title = element_text(hjust = 0.5),
    legend.text = element_text(size = 6),
    legend.position="bottom",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin=margin(-10, -10,-5,-5),
    strip.background = element_rect(colour = "black", fill = "gray90"),
    strip.text.x = element_text(color = "black", size = 8),
    strip.text.y = element_text(color = "black", size = 8, angle = 0)) 

p0 <- ggplot() +
    stat_ecdf(data = df_no,
              aes(PValue, group = t), 
              geom = "step", alpha = 1, size = 0.1,
              color = "black") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    labs(x = "p value", y = "") + 
    prettify  +
    theme(axis.title = element_text(size = 5),
          axis.text.x = element_text(size = 5, angle = 0),
          axis.text.y = element_text(size = 5))
p0 + prettify


ggplot() +
    stat_ecdf(data = df_no,
              aes(tP, group = t), 
              geom = "step", alpha = 1, size = 0.1,
              color = "black") +
    geom_abline(intercept = 0.5, slope = 0.5, linetype = "dashed") +
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    labs(x = "p value", y = "") + 
    prettify  +
    theme(axis.title = element_text(size = 5),
          axis.text.x = element_text(size = 5, angle = 0),
          axis.text.y = element_text(size = 5))
p0 + prettify
library(viridis)
library(ggrepel)
library(cowplot)
p1 <- ggplot() +
    stat_ecdf(data = df_yes, aes(adj.p, group = t, color = factor(t)),
              geom = "step") +
    geom_vline(xintercept = 0.05, linetype = "dashed") +
    # scale_color_viridis(discrete = TRUE, direction = -1, option ="C",
    #                     labels = c("leaf", "0.05", "0.1, 0.15", "0.2", 
    #                                "0.25", ">= 0.3")) +
    scale_color_manual(values = plasma(7), 
                        labels = c("leaf", "0.05", "0.1, 0.15", "0.2", 
                                   "0.25", ">= 0.3")) +
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    scale_x_continuous(trans = "sqrt", breaks = c(0.05, seq(0, 1, by = 0.2)),
                       labels = function(u) formatC(u)) +
    labs(x = "adjusted p value", y = "probablility",
         color = "levels") + prettify +
    theme(legend.position = "right",
          axis.text.x = element_text(angle = 0))
p1
gp <- ggdraw(p1) + draw_plot(p0, .28, .15, .55, .55) 
ggsave("~/Documents/ms/microbes/microbes_H01.pdf", gp, 
       units = "cm", width = 8, height = 6,
       dpi = 300, useDingbats = FALSE)
p0_1 <- ggplot() +
    stat_ecdf(data = df_no,
              aes(PValue, group = t, color = as.factor(t)), 
              geom = "step", alpha = 1, size = 0.3) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    guides(color=guide_legend(NULL, ncol= 2, byrow = FALSE,
                              override.aes = list(size = 1)))+
    labs(x = "p value", y = "", color = "Levels",
         title = "non-signal") + 
    prettify +
    scale_color_viridis_d(direction = -1, 
                          labels = c("leaf", as.character(seq(0.05, 1, by = 0.05)))) + 
    theme(legend.position = "right",
          plot.title = element_text(size = 8),
          axis.text.x = element_text(angle = 0),
          legend.key.size= unit(2, "mm"),
          plot.margin = margin(l = -3, b = -20, t = -10))

p1_1 <- ggplot() +
    stat_ecdf(data = df_yes, aes(adj.p, group = t, color = factor(t)),
              geom = "step") +
    geom_vline(xintercept = 0.05, linetype = "dashed") +
    # scale_color_viridis(discrete = TRUE, direction = -1, option ="C",
    #                     labels = c("leaf", "0.05", "0.1, 0.15", "0.2", 
    #                                "0.25", ">= 0.3")) +
    scale_color_viridis_d(direction = -1, 
                       labels = c("leaf", as.character(seq(0.05, 1, by = 0.05)))) +
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    scale_x_continuous(trans = "log10", breaks = c(0.05),
                       labels = function(u) formatC(u)) +
    labs(x = "log10(adj.p)", y = "",
         color = "Levels", title = "signal") + 
    prettify +
    theme(legend.position = "right",
          plot.title = element_text(size = 8),
          axis.text.x = element_text(angle = 0),
          legend.key.size= unit(2, "mm"),
          plot.margin = margin(l = -3, b = -20, t = -10))
p01_1 <- plot_grid(p0_1 + theme(legend.position = "none"), 
                   p1_1 + theme(legend.position = "none"),
                   labels = LETTERS[1:2])
leg <- get_legend(p1_1 + theme(legend.box.margin = margin(0, 0, 0, 0)))
sepF <- plot_grid(p01_1, leg, rel_widths = c(2.4, .5))
ggsave("~/Documents/ms/microbes/microbes_H01_1.pdf", sepF, 
       units = "cm", width = 10, height = 5,
       dpi = 300, useDingbats = FALSE)


# =============================================================
head(out)
pv0 <- out %>%
    filter(node %in% no_leaf) %>%
    mutate(s = ifelse(logFC < 0, -PValue, PValue)) %>%
    select(s) %>%
    unlist()
