require(readr)
require(dplyr)
require(stringr)
require(tidyr)

# load the chipseq genes
loaddir <- "genesets_literature/references/"

stat1 <- read.table(paste0(loaddir, "stat1_chipseq_mouse_tsskb.tsv"), header = T)

genes_stat1 <- stat1$Target_genes

# filter colnames
idx1 <- grepl(colnames(stat1), pattern = "CD4")
idx2 <- grepl(colnames(stat1), pattern = "Th1")

df1 = stat1[,idx1]
df2 = stat1[,idx2]

genes1 <- genes_stat1[df2$SRX183954.Th1_Cells > 150]
stat1_th1_targets <- as.data.frame(genes1) 
colnames(stat1_th1_targets) <- c("gene")

write.csv(stat1_th1_targets, "genesets_literature/references/stat1_th1_targets.csv", row.names = F)

# genes that are expressed above 200 in at least 5 studies
#genes2 <- genes_stat1[rowSums(df1>300)>4]
#stat1_CD4_targets <- as.data.frame(genes2)
#colnames(stat1_CD4_targets) <- c("gene")
#write.csv(stat1_CD4_targets, "genesets_literature/references/stat1_CD4_targets.csv", row.names = F)
