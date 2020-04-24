# function to get node at each of (4 levels)
fun_inner_label <- function(tree, tree_df, column_leaf) {
    trD <- ggtree(tr = tree)$data
    nodeList <- split(trD$node, trD$x)
    nodeL <- lapply(nodeList, unique)[as.character(1:4)]
    
    mT <- data.frame(matTree(tree))
    tree_df[["L1"]] <- transNode(tree = tree, node = as.character(tree_df[[column_leaf]]))
    tree_df <- tree_df %>%
        left_join(mT)
    list_df <- list("1" = tree_df %>% select(genomic_cluster, L4) %>% 
                        rename(node = L4),
                    "2" = tree_df %>% select(primary_transcript, L3) %>% 
                        rename(node = L3),
                    "3" = tree_df %>% select(miRNA, L2) %>% rename(node = L2),
                    "4" = tree_df %>% select(miRNA, L1) %>% rename(node = L1))
    out <- mapply(FUN = function(x, y){
        xy <- x %>% 
            filter(node %in% y) %>%
            distinct()
        vx <- xy$node
        names(vx) <- xy %>% select(-node) %>% unlist() %>% as.character()
        return(vx)
    }, x = list_df, y = nodeL)
    return(out)
}