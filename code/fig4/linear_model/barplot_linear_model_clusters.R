require(dplyr)
require(ggplot2)
require(tidyr)

df_clust = read.csv("data/clustered_genes/kinetic_genes_clusters_renamed_20211214.csv")
mygenes <- read.csv("data/linear_model/linear_model_processed.csv")


df <- left_join(mygenes, df_clust)
df <- df %>% filter(model == "Th1/2") %>% filter(celltype == "Th1/2")
g <- ggplot(data = df, aes(x= annot, fill = category)) +
  geom_bar() + theme_bw(base_size = 12) + theme(panel.grid.major = element_blank(),
                                                                             panel.grid.minor = element_blank())

ggsave("figures/fig4/linear_model/barplot_kinetic_clusters.png", width = 4, height = 2)