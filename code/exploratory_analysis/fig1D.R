
require(ggplot2)
require(readr)
require(dplyr)
source("code/paper_theme.R")


df <- read.csv("data/peine_logfc_tidy.csv")  
df$celltype[df$celltype == "ThMix"] = "Th1/2"


genes <- c("Ifng", "Gata3", "Tbx21", "Il4")

df <- df %>% filter(gene %in% genes)

# reorder 
df$gene <- factor(df$gene, levels = c("Tbx21", "Gata3", "Ifng", "Il4"))

p <- ggplot(data = df, aes(time, value, color = celltype))
p + geom_line(size = 0.5, show.legend = FALSE) + 
  geom_point(size = 2, show.legend = FALSE) +
  mytheme(base_size = 12) +
  xlab("time (h)") + 
  ylab("log2FC(day0)") +
  facet_wrap(~gene, ncol = 4) +
  theme(plot.margin=unit(c(0.2,0.5,0.1,0.1),"cm"),strip.text = element_text(face = "italic"))


ggsave("figures/Fig1D.png", width = 6, height = 1.8)

