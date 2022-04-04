require(dplyr)
require(tidyr)
require(readr)
require(pheatmap)
source("code/paper_theme.R")

# load data and add minimal of both pvalues for sorting
x1 <- read.csv("genesets_output/quantDEGS/degs_Th1_ThMix_cor0.3.csv")
x2 <- read.csv("genesets_output/quantDEGS/degs_Th2_ThMix_cor0.3.csv")

intersect_genes <- x1$gene[x1$gene %in% x2$gene]

# load timecourse data and sort by pnew
df <- read_csv("data/peine_logfc_tidy.csv")

df <- df[df$gene %in% intersect_genes,]


df <- df %>% filter(gene != "Nrp1.1")
#df <- df %>% separate(gene, into = c("gene", "iso"), sep = "\\.")

g <- ggplot(data = df, aes(x=time, y=value, colour = celltype))
g2 <- g + geom_point() + geom_line() + theme_bw(base_size = 10) + facet_wrap(~gene, scales = "free_y", ncol = 6) +
  ylab("logFC(day0") + xlab("time (h)")

sname <-"figures/fig3/timecourse_venn_deg_candidates.png"
ggsave(sname, bg = "white", width = 10, height = 5)