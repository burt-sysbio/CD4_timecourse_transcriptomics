require(readr)
require(dplyr)
require(stringr)
require(tidyr)

# load the chipseq genes
loaddir <- "genesets_literature/references/"

stat4 <- read.table(paste0(loaddir, "stat4_chipseq_mouse_tsskb.tsv"), header = T)
genes_stat4 <- stat4$Target_genes

# I looked it up. this should be
#SRX021624
#GSM550303: Stat4WTTh1 from Wei et al paper 
df <- stat4$SRX021624.Th1_Cells
# genes that are expressed above x in at least n studies

genes <- genes_stat4[df>10]

stat4_th1_targets <- as.data.frame(genes) 
colnames(stat4_th1_targets) <- c("gene")

write.csv(stat4_th1_targets, "genesets_literature/references/stat4_th1_targets.csv", row.names = F)
