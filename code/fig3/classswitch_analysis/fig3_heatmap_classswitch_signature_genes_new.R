
# signature genes
require(readr)
require(dplyr)
require(pheatmap)
library(cowplot)
require(tidyr)

source("code/utils.R")
source("code/paper_theme.R")

df_clust <- read.csv("data/clustered_genes/kinetic_genes_20210714_3clust.csv")

df_degs <- read.csv("data/classswitch_quantDEGS/quantDEGS_summary.csv")
df <- read_csv("data/peine_data_log2fc_annot.csv")
caro_genes <- read_csv("genesets/genesets_literature/genes_caro.csv", col_names = "gene")

# from caro genes only keep those that had been found kinetic in at least one cluster
df_clust <- df_clust %>% filter(gene %in% caro_genes$gene)
df_clust <- df_clust %>% pivot_wider(names_from = celltype, values_from = annot)

df_annot <- data.frame(cat = "ns", gene = df_clust$gene)
df_annot$cat[df_annot$gene %in% df_degs$gene] <- "quantDEG"


df_clust <- drop_genecol(df_clust)
df_annot <- drop_genecol(df_annot)


mycolors <- colorRampPalette(brewer.pal(n = 9, name = "Greys"))(4)
mycolors <- mycolors[2:4]
p <- pheatmap(df_clust,
             cluster_cols = F,
            cellwidth = 10,
         cellheight = 10,
         border_color = "black",
         color = mycolors,
         treeheight_row =  0,
         treeheight_col = 0,
         legend = FALSE)

save_pheatmap(p, "figures/supp/fig_S4_classswitch_quantDEG.svg")