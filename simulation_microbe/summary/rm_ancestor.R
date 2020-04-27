rm_ancestor <- function(node, tree) {
    if (!length(node)) {
        node
    } else {
        loc <- node
        # each row is a path connecting a leaf to the root
        mat <- matTree(tree = tree)
        # which paths are detected nodes on?
        lat <- lapply(loc, FUN = function(x){
            which(mat == x, arr.ind = TRUE)
        })
        lat1 <- do.call(rbind.data.frame, lat)
        lat1 <- lat1 %>% 
            group_by(row) %>% 
            arrange(col) %>%
            as.matrix()
        
        # nodes in paths where exist only one detected nodes
        lat2 <- lat1[!duplicated(lat1[, 1]), ]
        k2 <- unique(mat[lat2])
        # nodes in paths where exist more one detected nodes
        # (node that is closest to the leaf level not included)
        lat3 <- lat1[duplicated(lat1[, 1]), ]
        k3 <- unique(mat[lat3])
        
        # remove nodes 
        kk <- intersect(k2, k3)
        setdiff(loc, kk)  
    }
}