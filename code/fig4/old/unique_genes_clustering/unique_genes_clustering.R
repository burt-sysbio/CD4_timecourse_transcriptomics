require(pheatmap)
require(dplyr)
require(RColorBrewer)
require(stringr)
require(stringi)
require(cowplot)
source("code/paper_theme.R")
# define unqiue and hybrid genes for thm and th0
# approach: take th1 th2 genes as background
# approach2: filter genes with high cor across all replicates
# then do hclust and assign cluster names manually

# get data and degs th1 vs th2
df_cor <- readRDS("data/correlation_df.RDS")

mincor_df <- read.csv("data/minimal_correlation.csv")
mincor_df <- mincor_df %>% filter(sig == "sig")


# keep only genes that I used for minimal correlation volcano plot
df_filt <- df_cor[df_cor$Symbol %in% mincor_df$gene,]



# split df for thm and th0
df_thm <- df_filt[c("Symbol", "cor_th1_th2", "cor_th1_thm", "cor_th2_thm")]
df_th0 <- df_filt[c("Symbol", "cor_th1_th2", "cor_th1_th0", "cor_th2_th0")]


# take df, cluster, then reorder clusters and patch together with annotation
myfun <- function(df, 
                  contrast,
                  manual_clust = NULL,
                  nclust = 10){

  # run initial heatmap clust to get ordering, then annotate in second heatmap
  filename <- paste0("figures/fig3/clustering_degs/heatmaps/hm_", contrast, "_3.png")
  
  if(contrast=="Th0"){
    labels_col <- c("Th1vsTh2", "Th1vsTh0", "Th2vsTh0")
  } else{
    labels_col <- c("Th1vsTh2", "Th1vsTh1/2", "Th2vsTh1/2")
  }
  
  
  p <- pheatmap(df[, 2:ncol(df)], 
                scale = "row",
                cellwidth = 40,
                width = 7,
                color = heatmap_colors,
                height = 12,
                cluster_cols = F, 
                treeheight_row = 0,
                fontsize = 22,
                show_rownames = F,
                cutree_rows = nclust,
                labels_col = labels_col,
                filename = filename)
  
  clust <- p$tree_row$order
  out <- cutree(p$tree_row, k=nclust)
  df$cluster <- out
  
  df <- df[clust,]
  return(df)
}

nclust <- 10
out_th0 <- myfun(df_th0, contrast = "Th0", nclust = nclust)
# run again with manual clustering
write.csv(out_th0, "code/fig3/th0_hclust.csv")

out_thm <- myfun(df_thm, contrast = "Thm", nclust = nclust)
write.csv(out_thm, "code/fig3/thm_hclust.csv")


