# get unique and shared thm and th0 genes by running masigpro
# plot venn diagram th1vsthm th2vsthm
require(venn)

df <- readRDS("r_objects/fig2/degs_cor_masigpro.RDS")

# split into th0 and thm contrast

th1_thm <- df$Symbol[df$cell == "Th1_ThMix"]
th2_thm <- df$Symbol[df$cell == "Th2_ThMix"]

th1_th0 <- df$Symbol[df$cell == "Th1_Th0"]
th2_th0 <- df$Symbol[df$cell == "Th2_Th0"]

# th1 th2 vs thm
png(file="figures/fig2/venn_masigpro_cor_thm.png")
venn(list(th1_thm, th2_thm),
     snames = c("Th1vsThM", "Th2vsThM"),
     ilcs = 2.0,
     sncs = 3.0,
     box = F)
dev.off()


# th1 th2 vs thm
png(file="figures/fig2/venn_masigpro_cor_th0.png")
venn(list(th1_th0, th2_th0),
     snames = c("Th1vsTh0", "Th2vsTh0"),
     ilcs = 2.0,
     sncs = 3.0,
     box = F)
dev.off()

# save output to use as gene lists
thm_unique <- th1_thm[th1_thm %in% th2_thm]
th1_like <- th1_thm[!(th1_thm %in% th2_thm)]
th2_like <- th2_thm[!(th2_thm %in% th1_thm)]

#
th0_unique <- th1_th0[th1_th0 %in% th2_th0]
th1_like2 <- th1_th0[!(th1_th0 %in% th2_th0)]
th2_like2 <- th2_th0[!(th2_th0 %in% th1_th0)]

out_list <- list(thm_unique, th1_like, th2_like, th0_unique, th1_like2, th2_like2)
names <- c("thm_unique", "th1_thm", "th2_thm", "th0_unique", "th1_th0", "th2_th0")

for(i in seq_along(names)){
  out <- out_list[[i]]
  name <- names[i]
  file <- paste0("output_gene_lists/unique_genes/", name, ".txt")
  write.table(out, file = file, col.names = F, row.names = F, quote = F)
}