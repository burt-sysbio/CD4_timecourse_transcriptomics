
require(readr)
require(dplyr)
require(tidyr)
require(ggplot2)
source("code/paper_theme.R")
require(svglite)
require(cowplot)
# load tidy data
df_tidy <- read.csv("data/peine_logfc_tidy.csv")

# if you wan tto check candidates use this, otherwise df_down/trans$Symbol for filtering df_tidy to get overview
mygenes <- c("Eomes")


sname1 <- paste0("figures/fig3/class_switch_cartoon.png")
sname2 <- paste0("figures/fig3/class_switch_cartoon.svg")

  
df_tidy <- df_tidy %>% filter((gene %in% mygenes))
df_tidy$time <- as.numeric(df_tidy$time)

p1 <- ggplot(df_tidy, aes(time, value, color = celltype))
p1 + 
  geom_point(size = 0.8) +
  geom_line(size = 0.2) + 
  theme_bw(base_size = 10) +
  xlab("time (h)") + ylab("expr. norm.") + 
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_blank(),
        axis.title = element_text(size = 10),
        plot.margin=unit(c(0.2,0.5,0.1,0.1),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
        coord_cartesian(xlim = c(0,120), expand = T, ylim = c(-0.5,3)) +
  scale_x_continuous(breaks = c(0,60,120))


ggsave(sname1, width = 1.2, height = 1.0, bg = "white")
ggsave(sname2, width = 1.2, height = 1.0, bg = "white")

#
