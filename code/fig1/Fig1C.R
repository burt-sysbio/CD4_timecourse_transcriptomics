
# requirements
require(readr)
require(dplyr)
require(pheatmap)
require(tidyr)
source("code/paper_theme.R")
source("code/utils.R")


df <- read_csv("data/peine_data_log2fc_annot.csv")
# split isoforms
df <- df %>% separate(gene, into = c("candidate", "iso"), sep = "\\.", remove = F)


candidate_genes <- read.csv("genesets_literature/candidate_genes.csv", header = TRUE) %>% select(gene)

# filter dfs for candidate genes
df <- df %>% filter(gene %in% candidate_genes$gene)

# reorder columns
df_th0 <- df[, 1:10]
df_th1 <- df[, 11:20]
df_th2 <- df[, 21:30]
df_thm <- df[, 31:40]
df_list <- list(df_th0, df_th1, df_thm, df_th2)

# add gene names to reordered df
genes <- df$gene
df <- do.call("cbind", df_list)
rownames(df) <- genes

# sort df acc. to order in candidate list
df <- as.data.frame(df)
df <- df[match(candidate_genes$gene, rownames(df)),]
df <- na.omit(df)


# plotting parameters
width = 8
height = 6
breaks <- seq(-1.5,1.5, length.out = 101)
cellw <- 9
cellh <- 9
timepoints <- c(0,3,6,12,24,35,48,76,96,120)
labels_col <- as.character(rep(timepoints, 4))


p <- pheatmap(df,
         breaks = breaks,
         color = heatmap_colors,
         cluster_cols = F,
         cluster_rows = F,
         show_rownames = T,
         show_colnames = T,
         gaps_col = c(10,20,30),
         gaps_row = c(9,16,27),
         scale = "row",
         cellwidth = cellw,
         cellheight = cellh,
         labels_col = labels_col,
         width = width,
         height = height,
         legend = T,
         filename =  "figures/fig1C.png"
         )
