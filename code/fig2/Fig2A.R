
# make clustered heatmaps using masigpro for each cells kinetic genes
require(readr)
require(tidyr)
require(dplyr)
require(ggplot2)
require(pheatmap)
source("code/paper_theme.R")
source("code/utils.R")


# load data
df <- read.csv("data/peine_data_log2fc_annot.csv")
df <- drop_genecol(df)

# filter kinetic genes
df_kin <- read.csv("data/kinetic_genes/kinetic_genes_combined.csv")
# if a gene is kinetic in any cell type, keep it
kinetic_genes <- df_kin$gene[apply(df_kin[,2:5], 1, any)]

# reorder
df <- df[rownames(df) %in% kinetic_genes,]

df_th0 <- df[,1:10]
df_th1 <- df[,11:20]
df_th2 <- df[,21:30]
df_thm <- df[31:40]

nclust <- 3


out_th0 <- clusterfun(df_th0, nclust)
out_th1 <- clusterfun(df_th1, nclust)
out_th2 <- clusterfun(df_th2, nclust)
out_thm <- clusterfun(df_thm, nclust)

# plot the heatmaps
make_heatmap(df_th0, out_th0, nclust, "figures/fig2A_Th0.png")
make_heatmap(df_th1, out_th1, nclust, "figures/fig2A_Th1.png")
make_heatmap(df_th2, out_th2, nclust, "figures/fig2A_Th2.png")
make_heatmap(df_thm, out_thm, nclust, "figures/fig2A_Th12.png")
