# plot kinetic genes including class switches and generate gene sets for sebastian
require(venn)
require(VennDiagram)
require(dplyr)
require(tidyr)
source("code/fig3/utils.R")

# suppress log files from venn diag
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")


df2 <- readRDS("r_objects/fig3/cluster_assignment.RDS")
# save gene sets for sebastian
clusters <- seq(1:4)
cells <- c("Th1", "Th2","ThMix", "Th0")

for(cell in cells){
  for(c in clusters){
    df <- df2 %>% filter(celltype == cell, annot == c) %>% select(Symbol)
    filename <- paste0("output_gene_lists/clustering/kinetic_genes_", cell, "_cluster", c, ".txt")
    write.table(df, file = filename, col.names = F, row.names = F, quote = F)
  }
}

# make venn only for th1 and th2 + get class switches
th1 <- df2 %>% filter(celltype == "Th1")
th2 <- df2 %>% filter(celltype == "Th2")
df <- rbind(th1,th2)
df <- df %>% pivot_wider(values_from = annot, names_from = celltype)
df$switch <- df$Th1 != df$Th2
n_switch <- sum(df$switch, na.rm = T)

switch_genes <- df$Symbol[which(df$switch == T)]
shared_genes <- df$Symbol[which(!is.na(df$switch))]

# extract symbols
th1 <- th1$Symbol
th2 <- th2$Symbol
th1_unique <- th1[!(th1 %in% shared_genes)]
th2_unique <- th2[!(th2 %in% shared_genes)]

write.table(switch_genes, "output_gene_lists/clustering/kinetic_genes_th1_th2_switch.txt", row.names = F, col.names = F, quote = F)
write.table(shared_genes, "output_gene_lists/clustering/kinetic_genes_th1_th2_shared.txt", row.names = F, col.names = F, quote = F)
write.table(th1_unique, "output_gene_lists/clustering/kinetic_genes_th1_unique.txt", row.names = F, col.names = F, quote = F)
write.table(th2_unique, "output_gene_lists/clustering/kinetic_genes_th2_unique.txt", row.names = F, col.names = F, quote = F)

png(filename = "figures/fig3/venn_th1_th2_kinetic.png", width = 3, height = 3, units = "in", res = 600)
venn(x= list(th1,th2), snames = c("Th1", "Th2"), box = F, ilcs = 1.0)
dev.off()


# plot venn diag and get gene lists that switch class between celltypes
cell1 <- "Th1"
cell2 <- "Th2"
cell3 <- "ThMix"
l_thm <- pipeline(df2, cell1, cell2, cell3)


cell1 <- "Th1"
cell2 <- "Th2"
cell3 <- "Th0"
l_th0 <- pipeline(df2, cell1, cell2, cell3)

