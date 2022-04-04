require(readr)
require(dplyr)
require(stringr)
require(tidyr)

# load the chipseq genes
loaddir <- "genesets_literature/references/"

tbet_pos <- read.csv(paste0(loaddir, "zhu_tbet_targets_pos.csv"), col.names = "Tbet_pos", header = FALSE)
tbet_neg <- read.csv(paste0(loaddir, "zhu_tbet_targets_neg.csv"), col.names = "Tbet_neg", header = FALSE)


stat4_genes <- read.csv(paste0(loaddir, "stat4_targets_wei_2010.csv"), col.names = c("X", "Stat4"), header = FALSE, sep = "") %>% select(Stat4)
stat6_genes <- read.csv(paste0(loaddir, "stat6_targets_wei_2010.csv"), col.names = c("X", "Stat6"), header = FALSE, sep = "") %>% select(Stat6)

gata_th1 <- read.csv(paste0(loaddir, "wei_gata_targets_th1.csv"), col.names = "Gata_Th1", header = FALSE)
gata_th2 <- read.csv(paste0(loaddir, "wei_gata_targets_th2.csv"), col.names = "Gata_Th2", header = FALSE)

# from chipseqdb
stat1_genes <- read.csv("genesets_literature/references/stat1_th1_targets.csv", header = T, col.names = "Stat1")
#stat4_genes <- read.csv("genesets_literature/references/stat4_th1_targets.csv", header = T, col.names = "Stat4")
#stat6_genes <- read.csv("genesets_literature/references/stat6_th2_targets.csv", header = T, col.names = "Stat6")

#stat1_genes2 <- read.csv("genesets_literature/references/stat1_CD4_targets.csv", header = T, col.names = "Stat1_CD4")
tbet_all <- read.csv(paste0(loaddir, "zhu_tbet_all_wduplicates.csv"), col.names = "Tbet_large", header = FALSE) %>% unique()


chipseq_genes <- list(tbet_pos, tbet_neg, gata_th1, gata_th2, stat1_genes, stat4_genes, stat6_genes, tbet_all)

myfun <- function(df){
  df$gs_name <- colnames(df)[1]
  colnames(df) <- c("gene_symbol", "gs_name")
  df <- df[c("gs_name", "gene_symbol")]
  return(df)
}

chipseq_genes <- lapply(chipseq_genes, myfun)
chipseq_genes <- bind_rows(chipseq_genes)
chipseq_genes <- na.omit(chipseq_genes)

write.csv(chipseq_genes, "genesets_literature/references/chipseq_genes_all.csv", row.names = F)
