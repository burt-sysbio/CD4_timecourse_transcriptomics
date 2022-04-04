require(readr)
require(tidyr)
require(dplyr)

# take the genes that were assigned to different clusters and remove the isoform information
# needed to do pathway analysis to see pathways for individual clusters
df <- read.csv("data/clustered_genes/kinetic_genes_20210714_3clust.csv")

df_wide <- df %>% pivot_wider(names_from = celltype, values_from = annot)

df_clust1 <- df_wide[rowSums(df_wide[,2:5]==1)==4,]
df_clust2 <- df_wide[rowSums(df_wide[,2:5]==2)==4,]
df_clust3 <- df_wide[rowSums(df_wide[,2:5]==3)==4,]
gene_vec <- list(df_clust1, df_clust2, df_clust3)

sdir <- paste0("genesets_output/pathway_analysis/")

# df for individual clusters
for(i in seq(1:3)){
  df2 <- df %>% filter(annot == i)
  df2 <- df2 %>% separate(gene, c("gene", "isoform"), sep ="\\.") %>% select(gene)
  df2 <- unique(df2)
  sname <- paste0(sdir, "kinetic_genes_clust", i, ".txt")
  
  write.table(df2, sname, row.names = F, col.names = F, quote = F)
}

for(i in seq(1:3)){
  df3 <- gene_vec[[i]]
  df3 <- df3 %>% separate(gene, c("gene", "isoform"), sep ="\\.") %>% select(gene)
  df3 <- unique(df3)
  sname <- paste0(sdir, "kinetic_genes_clust_new", i, ".txt")
  
  write.table(df3, sname, row.names = F, col.names = F, quote = F)
}