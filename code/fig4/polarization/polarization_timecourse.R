
require(ggplot2)
require(readr)
require(dplyr)

source("code/paper_theme.R")


df <- readRDS("data/peine_logfc_tidy.RDS")  


df$celltype[df$celltype == "ThMix"] = "Th1/2"

genes1 <- c("Cxcr3", "Ly6c1", "Il12rb2")
genes2 <- c("Eomes", "Ifng", "Il18r1")
genes3 <- c("Ccr4", "Ifng", "Il18rap")



df <- df %>% filter(gene %in% genes1)


p <- ggplot(data = df, aes(time, value, color = celltype))
p + geom_line(size = 1) + mytheme() +
  xlab("time (h)") + ylab("expr. norm.") +
  facet_wrap(~gene, scales = "free_y", ncol = 3)

sdir <- "figures/fig3/timecourse_thm_th1like.png"
ggsave(sdir, width = 7.5, height = 2.2)
