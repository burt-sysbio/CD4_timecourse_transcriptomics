
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



make_heatmap <- function(df, out, nclust, filename){
  
  # reorder df based on correlation clustering
  df <- df[out$order, ]
  
  timepoints <- c(0,3,6,12,24,35,48,73,96,120)
  timepoints <- as.character(timepoints)
  
  # induce some gaps where clusters change
  test <- out$cut[out$order]
  idx <- which(diff(test)!=0)
  
  pheatmap(df,
           scale = "none",
           breaks = seq(-2,2, length.out = 101),
           color = heatmap_colors,
           cluster_rows = F,
           cluster_cols = F,
           show_rownames = F,
           cellheight = 0.03,
           cellwidth = 12,
           labels_col = timepoints,
           gaps_row = idx,
           annotation_legend = F,
           filename = filename)
}