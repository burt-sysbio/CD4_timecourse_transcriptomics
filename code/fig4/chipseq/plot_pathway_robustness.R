require(dplyr)
require(tidyr)
require(ggplot2)

threshold_array <- seq(0,0.5,0.01)

fnames <- paste0("code/fig4/chipseq/chipseq_pathway_results_degs_Th1_Th2_cor", threshold_array, ".csv")

filelist <- lapply(fnames, read.csv)

for(i in seq_along(filelist)){
  filelist[[i]]$cor <- threshold_array[i] 
}

approach <- "linear model"
df <- bind_rows(filelist)
df <- df[df$approach == approach,]
df <- df %>% mutate(pval2 = -log10(pvalue))

df$ID[df$ID== "Tbet_large"] <- "Tbet"

df <- df[!(df$ID %in% c("Tbet_neg", "Tbet_pos")),]

g <- ggplot(data = df, aes(x = cor, y = pval2))

g2 <-  g + geom_line(aes(colour=ID), size = 1.0) +
  xlab("Correlation filter") +
  ylab("-log10(pvalue)") +
  facet_wrap(~model, ncol = 4) +
  theme_bw(base_size = 12)

print(g2)


sname <- paste0("figures/fig4/chipseq/chipseq_pathways_robustness")
ggsave(paste0(sname, ".png"), width = 8, height = 2, bg = "white")
ggsave(paste0(sname, ".svg"), width = 8, height = 2, bg = "white")

