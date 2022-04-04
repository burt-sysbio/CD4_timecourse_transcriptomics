require(readr)
require(dplyr)
require(ggplot2)
require(tidyr)
require(cowplot)
# load data from correlation clustering and make barplot for assigned categories

df_th0 <- read.csv("data/correlation/hclust_output/hclust_corr_th0.csv")
df_thm <- read.csv("data/correlation/hclust_output/hclust_corr_thm.csv")

df_th0$contrast <- "Th0"
df_thm$contrast <- "Th1/2"

df_th0 <- df_th0[,c("contrast", "category", "gene")]
df_thm <- df_thm[,c("contrast", "category", "gene")]

ann_df <- bind_rows(df_th0, df_thm)

# store ann_df for further workings
write.csv(ann_df, "data/correlation/hclust_output/hclust_correlation_summary.csv", row.names = F)

# make output gene lists for sebastian
clusters <- c("Th1like", "Th2like", "Superpos", "Ind")
for(cl in clusters){
  df1 <- ann_df %>% filter(category == cl, contrast == "Th0") %>% 
    separate(gene, into = c("gene", "isoform"), sep = "\\.")  %>% select(gene)
  
  df2 <- ann_df %>% filter(category == cl, contrast == "Th1/2") %>% 
    separate(gene, into = c("gene", "isoform"), sep = "\\.")  %>% select(gene)
  
  df1 <- unique(df1)
  df2 <- unique(df2)

  mydir <- paste0("genesets_output/hclust_correlation/")
  file1 <- paste0(mydir, "hclust_corr_th0_", cl, ".txt")
  file2 <- paste0(mydir, "hclust_corr_thm_", cl, ".txt")
  write.table(df1, file1, row.names = F, col.names = F, quote = F)
  write.table(df2, file2, row.names = F, col.names = F, quote = F)
}
