
require(dplyr)
require(grid)
# utility functions used for plotting

drop_genecol <- function(df){
  df <- as.data.frame(df)
  stopifnot("gene" %in% colnames(df))
  genes <- df$gene
  df <- df %>% select(-gene)
  rownames(df) <- genes
  return(df)
}

# clustering based on correlation as implemented by MaSigPro
clusterfun <- function(clusterdata, k){
  dcorrel <- matrix(rep(1, nrow(clusterdata)^2), 
                    nrow(clusterdata), nrow(clusterdata)) - 
    cor(t(clusterdata), use = "pairwise.complete.obs")
  clust <- hclust(as.dist(dcorrel), method = "ward.D")
  
  cut <- cutree(clust, k = k)
  
  OUTPUT <- list(cut, clust$order)
  names(OUTPUT) <- c("cut", "order")
  OUTPUT
}