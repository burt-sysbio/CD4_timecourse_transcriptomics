require(dplyr)
require(tidyr)
require(pheatmap)
source("code/paper_theme.R")

# read timecourse df, chipseq results and linear model results
df_long <- read.csv("data/peine_logfc_tidy.csv")

df_chipseq <- read.csv("code/fig4/chipseq/chipseq_pathway_results.csv")

# only focus on stat6 and independent genes
model <- "Ind"
ID <- "Stat4"
mygenes <- df_chipseq %>% filter(ID == ID) %>% filter(model == model) %>% select(geneID)
mygenes <- strsplit(mygenes$geneID, "/")[[1]]

mydf <- df_long[df_long$gene %in% mygenes,]

sname <-paste0("figures/fig4/linear_model/timecourse_linear_model_", model, "_", ID, ".png")

# reorder df based on car gene sort


g <- ggplot(data = mydf, aes(x=time, y=value, colour = celltype))
g2 <- g + geom_point() + geom_line() + theme_bw(base_size = 10) + facet_wrap(~gene, scales = "free_y", ncol = 7)
#sname <-paste0("figures/fig4/linear_model/timecourse_linear_model_", modeltype, ".png")
ggsave(sname, bg = "white", width =12, height = 20)