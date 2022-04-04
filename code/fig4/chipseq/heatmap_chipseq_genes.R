require(dplyr)
require(tidyr)
require(readr)
require(pheatmap)
source("code/paper_theme.R")

# read data
df <- read_csv("data/peine_data_log2fc_annot.csv")
df_targets <- read.csv("genesets_literature/references/chipseq_genes_all.csv")

# get dfs for individual celltypes and rearrange columns
df_th0 <- df[, 1:10]
df_th1 <- df[, 11:20]
df_th2 <- df[, 21:30]
df_thm <- df[, 31:40]
df_list <- list(df_th0, df_th1, df_thm, df_th2)


genes <- df$gene
df <- do.call("cbind", df_list)

rownames(df) <- genes

# parameters for heatmap
breaks <- seq(-1.5,1.5, length.out = 101)
timepoints <- c(0,3,6,12,24,35,48,76,96,120)
labels_col <- as.character(rep(timepoints, 4))


df <- as.data.frame(df)


# get unique gene sets from chipseq data
targets <- unique(df_targets$gs_name)
for(target in targets){
  mygenes <-df_targets$gene_symbol[df_targets$gs_name == target]
  mydf <- df[rownames(df) %in% mygenes,]
  
  sname <-paste0("figures/fig4/chipseq/heatmaps/heatmap_chipseqtargets_", target)
  
  p <- pheatmap(mydf,
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

