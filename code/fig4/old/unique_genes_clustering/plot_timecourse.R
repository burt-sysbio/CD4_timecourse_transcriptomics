
# plot gene timecourse for different pathways and categories
require(readr)
require(dplyr)
require(tidyr)
require(ggplot2)
require(cowplot)

df <- readRDS("data/peine_tidy.RDS")

# get genes for unique/hybrid ThM (this is an updated object might need to check downstream code)
genes_clust <- readRDS("r_objects/fig2/unique_hybrid.RDS")
genes_clust$gene <- genes_clust$Symbol


gene_dir <- "gene_sets/references/"
genes_cytos <- read.csv(paste0(gene_dir, "cytokines.csv"), header = F, col.names = "gene")
genes_th1 <- read.csv(paste0(gene_dir, "stubbington_", "th1.csv"), header = F, col.names = "gene")
genes_th2 <- read.csv(paste0(gene_dir, "stubbington_", "th2.csv"), header = F, col.names = "gene")
genes_tf <- read.csv(paste0(gene_dir, "stubbington_", "tfs.csv"), header = F, col.names = "gene")
genes_receptors <- read.csv(paste0(gene_dir, "stubbington_", "receptors.csv"), header = F, col.names = "gene")

genes_cytos$category <- "Cytokine"
genes_tf$category <- "TF"
genes_receptors$category <- "Cyto_Receptor"

# can only do analysis with TH1/TH2 or cyto recep. TFs because they overlap
mygenes <- rbind(genes_cytos, genes_tf, genes_receptors)


#df <- df %>% filter(gene %in% mygenes$gene)


# first focus on hybrid
df_h <- genes_clust %>% filter((clust == "Hybrid") & (contrast == "Thm"))
df_u <- genes_clust %>% filter((clust == "Unique") & (contrast == "Thm"))

df_h <- df %>% filter(gene %in% df_h$Symbol)
df_u <- df %>% filter(gene %in% df_u$Symbol)


p1 <- ggplot(df_h, aes(time, value, color = celltype))
p1 + geom_line() +
  facet_wrap(~gene, scales = "free_y")
ggsave2("figures/fig2/timecourse_hybrid_all.png")

p2 <- ggplot(df_u, aes(time, value, color = celltype))
p2 + geom_line() +
  facet_wrap(~gene, scales = "free_y")
ggsave2("figures/fig2/timecourse_unique_all.png")

# choose some genes based on previous plot
genes_h <- c("Tbx21", "Stat1", "Ccr8")
genes_u <- c("Il21r", "Hsf2", "Twist1")

df_h2 <- df %>% filter(gene %in% genes_h)
df_u2 <- df %>% filter(gene %in% genes_u)


p3 <- ggplot(df_h2, aes(time, value, color = celltype))
p3 <- p3 + geom_line(size = 1.2) + theme_bw()+
  facet_wrap(~gene, scales = "free_y", ncol = 3)+
  ylab("expr.") +
  theme(text = element_text(size = 15), panel.grid.minor = element_blank(), legend.position = "none")


p4 <- ggplot(df_u2, aes(time, value, color = celltype))
p4 <- p4 + geom_line(size = 1.2) + theme_bw()+
  facet_wrap(~gene, scales = "free_y", ncol = 3)+
  ylab("expr.") +
  theme(text = element_text(size = 15), panel.grid.minor = element_blank(), legend.position = "none")

plot_grid(p3, p4, nrow = 2,align = "hv")
ggsave2("figures/fig2/timecourse_cluster_subset.png", width = 6, height = 4)

