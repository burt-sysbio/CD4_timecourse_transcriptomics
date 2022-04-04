require(readr)
require(dplyr)
require(stringr)
require(tidyr)

# load the chipseq genes
loaddir <- "genesets_literature/references/"

stat6 <- read.table(paste0(loaddir, "stat6_chipseq_mouse_tsskb.tsv"), header = T)
genes_stat6 <- stat6$Target_genes

# I looked it up. this should be
#SRX021632
#GSM550311: Stat6WTTh2 
df = stat6$SRX021632.Th1_Cells

# genes that are expressed above x in at least n studies
genes <- genes_stat6[df>10]

stat6_th1_targets <- as.data.frame(genes) 
colnames(stat6_th1_targets) <- c("gene")
write.csv(stat6_th1_targets, "genesets_literature/references/stat6_th2_targets.csv", row.names = F)
