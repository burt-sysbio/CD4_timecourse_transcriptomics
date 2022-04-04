require(pheatmap)
require(dplyr)
require(RColorBrewer)
require(stringr)
require(stringi)
require(cowplot)
require(VennDiagram)
require(venn)
source("code/paper_theme.R")
# define unqiue and hybrid genes for thm and th0
# approach: take th1 th2 genes as background
# approach2: filter genes with high cor across all replicates
# then do hclust and assign cluster names manually

# get data and degs th1 vs th2
df_cor <- readRDS("data/correlation_df.RDS")

mincor_df <- read.csv("data/minimal_correlation.csv")
mincor_df <- mincor_df %>% filter(sig == "sig")

# get pvals from masigpro run, use both thm and th0 as contrast
res <- readRDS("code/fig2/mspro_tstep_contrast_Th1.RDS")
res2 <- readRDS("code/fig2/mspro_tstep_contrast_Th2.RDS")

# get degs vs thm (or Th0)
sigs <- get.siggenes(res, rsq = 0.7, vars = "groups")
sigs2 <- get.siggenes(res2, rsq = 0.7, vars = "groups")

th1_thm <- sigs$summary$ThMixvsTh1
th2_thm <- sigs2$summary$ThMixvsTh2

th1_th0 <- sigs$summary$Th0vsTh1
th2_th0 <- sigs2$summary$Th0vsTh2

# reduce to keep only degs that are also different in Th1/Th2 + correlation filter
th1_thm <- th1_thm[th1_thm %in% mincor_df$gene]
th2_thm <- th2_thm[th2_thm %in% mincor_df$gene]
th1_th0 <- th1_th0[th1_th0 %in% mincor_df$gene]
th2_th0 <- th2_th0[th2_th0 %in% mincor_df$gene]

# suppress log files from venn diag
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

png(filename = "figures/fig3/venn_thm_mincor.png", width = 3, height = 3, units = "in", res = 600)
venn(x= list(th1_thm,th2_thm), snames = c("vsTh1", "vsTh2"), box = F, ilcs = 1.0)
dev.off()

png(filename = "figures/fig3/venn_th0_mincor.png", width = 3, height = 3, units = "in", res = 600)
venn(x= list(th1_th0,th2_th0), snames = c("vsTh1", "vsTh2"), box = F, ilcs = 1.0)
dev.off()

# note: not intuitive: genes that are specific to th1_thm are th2 like because its DEGS
th2_like <- th1_thm[!(th1_thm %in% th2_thm)]
th1_like <- th2_thm[!(th2_thm %in% th1_thm)]
thm_intersect <- th1_thm[th1_thm %in% th2_thm]

th2_like2 <- th1_th0[!(th1_th0 %in% th2_th0)]
th1_like2 <- th2_th0[!(th2_th0 %in% th1_th0)]
th0_intersect <- th1_th0[th1_th0 %in% th2_th0]


sdir <- "pathway_analysis/my_genelists/venn_mincor/"
write.table(th2_like, paste0(sdir, "th2_like", ".txt"), row.names = F, col.names = F, quote = F)
write.table(th1_like, paste0(sdir, "th1_like", ".txt"), row.names = F, col.names = F, quote = F)
write.table(th1_like2, paste0(sdir, "th1_like2", ".txt"), row.names = F, col.names = F, quote = F)
write.table(th2_like2, paste0(sdir, "th2_like2", ".txt"), row.names = F, col.names = F, quote = F)
write.table(th2_like, paste0(sdir, "th2_like", ".txt"), row.names = F, col.names = F, quote = F)
write.table(thm_intersect, paste0(sdir, "thm_intersect", ".txt"), row.names = F, col.names = F, quote = F)
write.table(th0_intersect, paste0(sdir, "th0_intersect", ".txt"), row.names = F, col.names = F, quote = F)

