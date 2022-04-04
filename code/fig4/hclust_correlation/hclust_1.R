require(pheatmap)
require(dplyr)
require(RColorBrewer)
require(stringr)
require(stringi)
require(tidyr)
require(cowplot)
source("code/paper_theme.R")
source("code/utils.R")
require(grid)


# take df, cluster, then reorder clusters and patch together with annotation
myfun <- function(df, 
                  contrast,
                  nclust = 10,
                  labels = NULL,
                  cluster_labels = NULL){
  
  # assign labels
  if(contrast=="Th0"){
    labels_col <- c("Th1vsTh0", "Th2vsTh0")
  } else{
    labels_col <- c("Th1vsTh1/2", "Th2vsTh1/2")
  }
  
  mycolors <- heatmap_colors
  mybreaks <- seq(-1,1,length.out = 101)
  
  # annotate rownames
  rownames(df) <- df$gene
  labels_row <- rownames(df)
  if(!is.null(labels)){
    labels_row[!(labels_row %in% labels)] = ""
    show_rownames <- T
  } else {
    show_rownames <- F
  }  
  
  # filenames
  f0 <- paste0("figures/fig4/clustering_correlation/hclust_corr_", contrast)
  filename1 <- paste0(f0, ".png")
  filename2 <- paste0(f0, ".svg")
  
  p <- pheatmap(df[, 3:ncol(df)], 
                scale = "none",
                cellwidth = 12,
                cellheight = 0.5,
                breaks = mybreaks,
                color = mycolors,
                cluster_cols = F, 
                treeheight_row = 0,
                labels_row = labels_row,
                show_rownames = show_rownames,
                cutree_rows = nclust,
                labels_col = labels_col,
                fontsize = 8,
                fontsize_row = 8,
                filename = filename1,
                height = 4.5,
                width = 1.2)
  
  save_pheatmap(p, filename2, width =2, height = 5)
  
  
  df <- sort_clusters(df, p, nclust, cluster_labels)

  
  return(df)
}


sort_clusters <- function(df, p, nclust, cluster_labels){
  # take data frame and pheatmap object and reassign cluster labels
  
  # get cluster info from pheatmap object and reorder data frame accordingly
  clust <- p$tree_row$order
  out <- cutree(p$tree_row, k=nclust)
  df$cluster <- out
  
  df <- df[clust,]
  
  # create dummy df to rename the clusters from top to bottom in asc. order
  clust_order <- unique(df$cluster)
  clust_neworder <- seq(1:length(clust_order))
  
  dummy <- data.frame(cluster = clust_order, new_cluster = clust_neworder)
  
  if(!is.null(cluster_labels)){
    print("manual cluster labels should be checked for each newly provided cluster gene list!")
    stopifnot(length(cluster_labels) == length(clust_order))
  
    dummy$category <- cluster_labels
  }
  print(colnames(df))
  print(colnames(dummy))
  df <- dplyr::left_join(df, dummy, by = "cluster")
  return(df)
}

# remove specific cluster and plot again
myfun2 <- function(df, contrast, clust_id, df_candidates, nclust = 10, labels = NULL, cluster_labels = NULL){
  df <- df[!(df$new_cluster %in% clust_id),]
  
  # get annotation from df candidates
  df_annot <- add_target_annot(df, df_candidates)
  rownames(df) <- df$gene
  # further process df
  #df <- drop_genecol(df)
  #df <- df[c(2,3)]
  
  mycolors <- heatmap_colors
  mybreaks <- seq(-1,1,length.out = 101)
  
  f0 <- paste0("figures/fig4/clustering_correlation/hclust_corr_without_superpos_", contrast)
  filename1 <- paste0(f0, ".png")
  filename2 <- paste0(f0, ".svg")
  
  
  labels_row <- rownames(df)
  if(!is.null(labels)){
    labels_row[!(labels_row %in% labels)] = ""
    show_rownames <- T
  } else {
    show_rownames <- F
  }  
  
  mycols <- c("0" = "white", "1" = "darkslategrey")
  ann_colors = list(
    Gata_Th1 = mycols,
    Gata_Th2 = mycols,
    Tbet_pos = mycols,
    Tbet_neg = mycols,
    Tbet_Gata_target = mycols
  )
  
  p <- pheatmap(df[c(3,4)], 
                scale = "none",
                breaks = mybreaks,
                color = mycolors,
                cluster_cols = F, 
                cluster_rows = T,
                cutree_rows = nclust,
                cellwidth = 16,
                annotation_colors = ann_colors,
                annotation_row = df_annot,
                cellheight = 3,
                treeheight_row = 0,
                labels_row = labels_row,
                show_rownames = show_rownames,
                fontsize = 8,
                fontsize_row = 6,
                filename = filename1)
  
  save_pheatmap(p, filename2, width =3, height = 4)
  df <- sort_clusters(df, p, nclust, cluster_labels)
  return(df)
}

add_target_annot <- function(df, df_candidates){
  
  # input should be tidy df_candidate frame wit columns gs_name and gene
  # process candidate dataframe and make wide
  # returns annotated data frame matching gene names from df with annotation for the candidates
  df_deg <- df_candidates %>% filter(gene %in% df$gene)
  df_deg$isdeg <- 1
  df_deg <- df_deg %>% pivot_wider(names_from = gs_name, values_from = isdeg, values_fill = 0)
  
  # get genes that are candidates but not in df_clust (and therefore not kinetic)
  kin_genes <- df$gene[!(df$gene %in% df_deg$gene)]
  
  # make matrix filled with zeros for kinetic genes
  df_kin <- data.frame(matrix(0, ncol(df_deg)-1, nrow = length(kin_genes)))
  df_kin <- df_kin %>% mutate(gene = kin_genes) %>% select(gene, everything())
  colnames(df_kin) <- colnames(df_deg)

  df_annot <- bind_rows(df_deg, df_kin)
  df_annot <- drop_genecol(df_annot)
  return(df_annot)
}
# define unqiue and hybrid genes for thm and th0
# approach: take th1 th2 DEGS as background (these genes have already been filtered for correlation)

# load data
df_cor <- read.csv("data/correlation/correlation_df.csv")
cor_thres <- "0.3"
degs <- read.csv(paste0("data/classswitch_quantDEGS/quantDEGS_summary_cor", cor_thres, ".csv"))

# keep only th1 th2 degs
degs_th1th2 <- degs %>% filter(deg == "th1_th2")

# keep only genes that I used for minimal correlation volcano plot
df_filt <- df_cor[df_cor$gene %in% degs_th1th2$gene,]

# split df for thm and th0
df_thm <- df_filt[c("gene", "cor_th1_th2", "cor_th1_thmix", "cor_th2_thmix")]
df_th0 <- df_filt[c("gene", "cor_th1_th2", "cor_th1_th0", "cor_th2_th0")]



genes_caro <- read.csv("genesets_literature/candidate_genes_sorted_short.csv", header = T)
candidates <- genes_caro$gene

#nclust <- 20
nclust <- 10
# you can provide cluster labels but these should be added after running it without additional cluster labels

use_clust_annot <- F
if(use_clust_annot){
  df_clust_annot <- read.csv("code/fig4/hclust_correlation/hclust_annot.csv")
  clabels_th0 <- df_clust_annot[,1]
  clabels_thm <- df_clust_annot[,2]
} else{
  clabels_th0 <- NULL
  clabels_thm <- NULL
}


out_th0 <- myfun(df_th0, contrast = "Th0", nclust = nclust, labels = NULL, cluster_labels = NULL)
out_thm <- myfun(df_thm, contrast = "ThMix", nclust = nclust, labels = NULL, cluster_labels = NULL)

#sdir <- "data/correlation/hclust_output/"

#s_th0 <- paste0(sdir, "hclust_corr_th0.csv")
#s_thm <- paste0(sdir, "hclust_corr_thm.csv")
#write.csv(out_th0, s_th0, row.names = F)
#write.csv(out_thm, s_thm, row.names = F)
