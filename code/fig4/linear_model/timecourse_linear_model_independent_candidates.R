require(dplyr)
require(tidyr)
require(readr)
library(Cairo)

require(pheatmap)
source("code/paper_theme.R")


# load data and add minimal of both pvalues for sorting
mygenes <- c("Nrp1", "Serpina3f", "Ube2l6", "Oasl1.2", "Usp18", "Bst2")


df <- read_csv("data/peine_logfc_tidy.csv")
df <- df[df$gene %in% mygenes,]

df <- df %>% separate(gene, into = c("gene", "ISO"), sep ="\\.")

g <- ggplot(data = df, aes(x=time, y=value, colour = celltype))
g2 <- g + geom_point(size = 1) + geom_line(size=0.5) + theme_bw(base_size = 10) + facet_wrap(~gene, scales = "free_y", ncol = 6) +
  ylab("log2FC(day0)") + xlab("time (h)") +
  theme(strip.text = element_text(face = "italic"))

sname <-paste0("figures/fig4/linear_model/timecourses/timecourse_linmo_ind_candidates.png")
ggsave(sname, bg = "white", width =9, height = 1.5, type = "cairo-png")


# load data and add minimal of both pvalues for sorting
mygenes <- c("Nrp1", "Serpina3f", "Ube2l6", "Bst2")


df <- df[df$gene %in% mygenes,]

df <- df %>% separate(gene, into = c("gene", "ISO"), sep ="\\.")

g <- ggplot(data = df, aes(x=time, y=value, colour = celltype))
g2 <- g + geom_point(size = 1) + geom_line(size=0.5) + theme_bw(base_size = 10) + facet_wrap(~gene, scales = "free_y", ncol = 4) +
  ylab("log2FC(day0)") + xlab("time (h)") +
  theme(strip.text = element_text(face = "italic"))

sname <-paste0("figures/fig4/linear_model/timecourses/timecourse_linmo_ind_candidates2.png")
ggsave(sname, bg = "white", width =7, height = 1.7, type = "cairo-png")
