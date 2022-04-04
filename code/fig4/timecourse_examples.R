
# timecourse for class switches
require(ggplot2)
require(dplyr)
require(tidyr)
require(readr)
source("code/paper_theme.R")
  
# change this between df trans and df down
# load tidy data
df_tidy <- read.csv("data/peine_logfc_tidy.csv")

# load info
df_categories = read.csv("data/linear_model/linear_model_processed_degs_Th1_Th2_cor0.3.csv")


# if you wan tto check candidates use this, otherwise df_down/trans$Symbol for filtering df_tidy to get overview
mygenes <- c("Ly6a", "Eomes", "Ifng", "Ube2l6")
categories <- c("Th1-like", "Th2-like", "Superposition", "Independent")

annot <- data.frame(mygenes, categories)
colnames(annot) <- c("gene", "category")

df_tidy <- df_tidy %>% filter(gene %in% mygenes)
df_tidy$time <- as.numeric(df_tidy$time)


# add category
df_tidy <- left_join(df_tidy, annot, on = "gene")

p1 <- ggplot(df_tidy, aes(time, value, color = celltype))
p1 + 
  geom_line(size = 0.5, show.legend = F) + 
  geom_point(size = 1.5, show.legend = F) +
  facet_wrap(~category, ncol = 4, scales = "free_y") +
  xlab("time (h)") +
  ylab("log2FC(d0)") +
  scale_x_continuous(breaks = c(0,60,120)) +
  mytheme() +
  theme(plot.margin=unit(c(0.2,0.5,0.1,0.1),"cm"))

ggsave("figures/fig4/timecourse_examples.png", width = 8, height = 2.1, bg = "white")
#
