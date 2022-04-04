# load data and add minimal of both pvalues for sorting

require(dplyr)
require(pheatmap)
source("code/paper_theme.R")
baseline <- "cor0"
x1 <- read.csv(paste0("genesets_output/quantDEGS/degs_Th1_ThMix_", baseline, ".csv"))
x2 <- read.csv(paste0("genesets_output/quantDEGS/degs_Th2_ThMix_", baseline, ".csv"))

mygenes <- intersect(x1$gene, x2$gene)

plot_heatmap <- function(baseline, mygenes){
 
  df <- read_csv("data/peine_data_log2fc_annot.csv")
  # get dfs for individual celltypes
  df_th0 <- df[, 1:10]
  df_th1 <- df[, 11:20]
  df_th2 <- df[, 21:30]
  df_thm <- df[, 31:40]
  df_list <- list(df_th0, df_th1, df_thm, df_th2)
  
  genes <- df$gene
  df <- do.call("cbind", df_list)
  rownames(df) <- genes
  df <- as.data.frame(df)
  
  breaks <- seq(-1,1, length.out = 101)
  timepoints <- c(0,3,6,12,24,35,48,76,96,120)
  labels_col <- as.character(rep(timepoints, 4))
  
  
  sname <-paste0("figures/fig3/heatmap_venn_categories_", baseline)
  df2 <- df[rownames(df) %in% mygenes,]
  p <- pheatmap(df2,
                breaks = breaks,
                color = heatmap_colors,
                cluster_cols = F,
                cluster_rows = T,
                show_rownames = F,
                show_colnames = F,
                labels_col = labels_col,
                gaps_col = c(10,20,30),
                fontsize = 8,
                fontsize_row = 8,
                scale = "row",
                legend = T,
                filename =  paste0(sname, ".png")
  )
}
plot_heatmap(baseline, mygenes)