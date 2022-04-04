
# signature genes with annotation kinetic / not kinetic
require(readr)
require(dplyr)
require(pheatmap)
library(cowplot)
require(tidyr)
require(readxl)

source("code/paper_theme.R")
source("code/utils.R")


df <- read_csv("data/peine_data_log2fc_annot.csv")

# make a table and print for each DEG type the number of DEGs and class switches
df_switch <- read.csv("data/classswitch_quantDEGS/class_switch_summary_3clust.csv")
df_switch <- df_switch %>% separate(gene, into = c("gene", "iso"), sep = "\\.") %>% distinct(gene, .keep_all = TRUE)

# remove the switch for th0 vs thm
df_switch <- df_switch %>% select(-s_th0_thm)


contrast <- "qualDEG"
stopifnot(contrast %in% c("quantDEG", "qualDEG"))

if(contrast == "quantDEG"){
  df_deg <- read_xlsx("figures/supp/table_S4_QuantDEG_switch.xlsx")
  cellh <- 1
  sname <-"figures/fig3/heatmap_classswitches_quantDEG"
  
  width = 5
  height = 4
  
} else{
  df_deg <- read_xlsx("figures/supp/table_S3_QualDEG_switch.xlsx")
  cellh <- 1
  sname <-"figures/fig3/heatmap_classswitches_qualDEG"
  
  width = 4
  height = 3.5
}


genes_switch <- unlist(df_deg, use.names=FALSE)

# filter dfs for caros genes but also remove some because of space constraints
df <- df %>% filter(gene %in% genes_switch)
df_switch <- df_switch %>% filter(gene %in% genes_switch)
df_switch <- df_switch %>% filter(gene %in% df$gene)

# only focus on switch as index
df_switch <- drop_genecol(df_switch)
df_switch <- df_switch[, 6:10]

# convert boolean to integer
df_switch <- df_switch*1

# get dfs for individual celltypes
df_th0 <- df[, 1:10]
df_th1 <- df[, 11:20]
df_th2 <- df[, 21:30]
df_thm <- df[, 31:40]
df_list <- list(df_th0, df_th1, df_thm, df_th2)

# get hclust ordering for th1
genes <- df$gene
df <- do.call("cbind", df_list)
rownames(df) <- genes

breaks <- seq(-1.5,1.5, length.out = 101)
cellw <- 4

timepoints <- c(0,3,6,12,24,35,48,76,96,120)
labels_col <- as.character(rep(timepoints, 4))

df <- as.data.frame(df)

# color annotation for clusters
mycolors <- colorRampPalette(brewer.pal(n = 9, name = "Greys"))(4)
mycols <- c("0" = "white", "1" = "black")
ann_colors = list(
  s_th1_th0 = mycols,
  s_th1_th2 = mycols,
  s_th1_thm = mycols,
  s_th2_th0 = mycols,
  s_th2_thm = mycols
)




newnames <- lapply(
  rownames(df),
  function(x) bquote(italic(.(x))))

# rename cells
#colnames(df_annot) <- c("Th1vsTh0", "Th1vsTh2", "Th1vsTh1/2", "Th2vsTh0", "Th2vsTh1/2")
p <- pheatmap(df,
              breaks = breaks,
              color = heatmap_colors,
              cluster_cols = F,
              cluster_rows = T,
              show_rownames = F,
              show_colnames = F,
              gaps_col = c(10,20,30),
              treeheight_row = 0,
              #labels_row = as.expression(newnames),
              fontsize = 8,
              fontsize_row = 8,
              scale = "row",
              annotation_row = df_switch,
              annotation_colors = ann_colors,
              cellwidth = cellw,
              cellheight = cellh,
              labels_col = labels_col,
              annotation_legend = F,
              width = width,
              height = height,
              legend = T,
              filename =  paste0(sname, ".png")
)

save_pheatmap(p, paste0(sname, ".svg"), width = width, height = height)