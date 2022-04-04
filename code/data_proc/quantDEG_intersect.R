# take DEGs from Th1 Thm and Th2 Thm and intersect
require(readr)
require(tidyr)

contrast <- "ThMix"

cor_thres <- "0.3"
loaddir <- paste0("genesets_output/quantDEGS/", "degs_cleaned_")
df_th1 <- read.csv(paste0(loaddir,"Th1_", contrast, "_cor", cor_thres, ".txt"), header = F, col.names = "gene")
df_th2 <- read.csv(paste0(loaddir,"Th2_", contrast, "_cor", cor_thres, ".txt"), header = F, col.names = "gene")

genes_th1 <- df_th1$gene
genes_th2 <- df_th2$gene

genes_superpos <- genes_th1[genes_th1 %in% genes_th2]

th2_like <- genes_th1[!(genes_th1 %in% genes_th2)]
th1_like <- genes_th2[!(genes_th2 %in% genes_th1)]


f1 <- paste0(loaddir, "degs_th1_like_", contrast, ".txt")
f2 <- paste0(loaddir, "degs_th2_like_", contrast, ".txt")
f3 <- paste0(loaddir, "degs_superpos_", contrast, ".txt")

write.table(th1_like, f1, row.names = F, col.names = F, quote = F)
write.table(th2_like, f2, row.names = F, col.names = F, quote = F)
write.table(genes_superpos, f3, row.names = F, col.names = F, quote = F)

