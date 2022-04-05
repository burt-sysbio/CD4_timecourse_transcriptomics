
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



get_volcano_sigs <- function(comparison, sigs, sigs2){
  if(comparison == "Th1_Th2"){
    sigs <- sigs$Th2vsTh1
    myfilter <- "cor_th1_th2"
  } else if (comparison == "Th1_ThMix"){
    sigs <- sigs$ThMixvsTh1
    myfilter <- "cor_th1_thmix"
  } else if (comparison == "Th1_Th0"){
    sigs <- sigs$Th0vsTh1
    myfilter <- "cor_th1_th0"
  } else if (comparison == "Th2_Th0"){
    sigs <- sigs2$Th0vsTh2
    myfilter <- "cor_th2_th0"
  } else if (comparison == "Th2_ThMix"){
    sigs <- sigs2$ThMixvsTh2
    myfilter <- "cor_th2_thmix"
  }
  return(list(sigs, myfilter))
}

# functions for correlation volcano plot

prep_volcano_data <- function(comparison, pvals, df_cor, sigs, sigs2, candidates, cor_thres){
  
  # helper function to get comparisons right
  comp_info <- get_volcano_sigs(comparison, sigs, sigs2)
  myfilter <- comp_info[[2]]
  sigs <- comp_info[[1]]
  
  # filter data
  df_cor <- df_cor %>% filter(correlation == myfilter)
  df_pvals <- prep_volcano_pvals(pvals)
  
  # combine with pvals
  df_join <- dplyr::inner_join(df_cor, df_pvals)
  
  # add some output annotations
  out <- proc_volcano_data(df_join, sigs, candidates, cor_thres)
  return(out)
}

# process pvalues into data frame
prep_volcano_pvals <- function(pvals){
  p_adj <- pvals$p.vector
  gene_names <- rownames(pvals$p.vector)
  df_pvals <- data.frame(gene_names, p_adj)
  colnames(df_pvals) <- c("gene", "padj")
  return(df_pvals)
}

# add annotation dfs and additional readouts for coloring volcano plot
proc_volcano_data <- function(df_join, sigs, candidates, cor_thres){
  
  df_join <- df_join %>% mutate(cor_idx = 1 - value)
  
  # assign signif genes based on pval and correlation threshold
  df_join <- df_join %>% mutate(sig = ifelse((cor_idx >= cor_thres) & (padj < 0.05), "sig", "ns"))
  
  # add color based on DEG information from masigpro full run
  df_join$deg <- "ns"
  # if found in sig array, annotate
  df_join$deg[df_join$gene %in% sigs] <- "quant. DEG"
  df_join$deg[(df_join$gene %in% sigs) & (df_join$sig == "sig")] <- "qual. DEG"
  
  # convert p values to -log10
  df_join <- df_join %>% mutate(padj_log10 = -log10(padj))
  
  return(df_join)
}
