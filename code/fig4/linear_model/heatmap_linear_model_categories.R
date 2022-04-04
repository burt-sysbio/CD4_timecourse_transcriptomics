require(pheatmap)
require(dplyr)
require(tidyr)
require(readr)
source("code/paper_theme.R")

plot_heatmap_categories <- function(baseline){
  mygenes <- read.csv(paste0("data/linear_model/linear_model_processed_", baseline, ".csv"))
  mygenes <- mygenes %>% filter(model == "Th1/2") 
  categories <- unique(mygenes$category)
  
  
  df <- read_csv("data/peine_data_log2fc_annot.csv")
  df <- df %>% separate(gene, into = c("gene", "iso"), sep = "\\.") %>% distinct(gene, .keep_all = T)
  # get dfs for individual celltypes
  df_th0 <- df[, 1:10]
  df_th1 <- df[, 11:20]
  df_th2 <- df[, 21:30]
  df_thm <- df[, 31:40]
  df_list <- list(df_th0, df_th1, df_thm, df_th2)
  
  
  genes <- df$gene
  df <- do.call("cbind", df_list)
  
  rownames(df) <- genes
  
  breaks <- seq(-1,1, length.out = 101)
  timepoints <- c(0,3,6,12,24,35,48,76,96,120)
  labels_col <- as.character(rep(timepoints, 4))
  
  
  
  df <- as.data.frame(df)
  
  
  for(cat in categories){
    sname <-paste0("figures/fig4/linear_model/heatmaps/hm_linmo_", baseline, "_", cat)
    cat_genes <- mygenes$gene[mygenes$category == cat]
    df2 <- df[rownames(df) %in% cat_genes,]
    
    # reorder df based on car gene sort
    if(nrow(df2)<100){
      show_rownames <- T
      cellwidth = 8
      cellheight = 8
    } else {
      show_rownames <- F
      cellwidth = NA
      cellheight = NA
    }
    
    
    
    p <- pheatmap(df2,
                  breaks = breaks,
                  color = heatmap_colors,
                  cluster_cols = F,
                  cluster_rows = T,
                  cellwidth = cellwidth,
                  cellheight = cellheight,
                  show_rownames = show_rownames,
                  treeheight_row = 0,
                  show_colnames = T,
                  labels_col = labels_col,
                  gaps_col = c(10,20,30),
                  fontsize = 8,
                  fontsize_row = 8,
                  scale = "row",
                  legend = T,
                  filename =  paste0(sname, ".png")
    )
    
  }
}
plot_heatmap_categories("degs_Th1_Th2_cor0.3")