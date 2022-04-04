
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
mygenes <- c("Eomes", "Ifng", "Batf", "Il12rb2", "Ly6c1", "Il4", "Il18rap",
             "S100a1", "Gata3", "Hopx", "Stat1", "Il18r1", "Cxcr3", "Ccr6")


sname1 <- paste0("figures/supp/fig_S4/class_switch_candidates.png")
sname2 <- paste0("figures/supp/fig_S4/class_switch_candidates.svg")

  
df_tidy <- df_tidy %>% filter((gene %in% mygenes))
df_tidy$time <- as.numeric(df_tidy$time)

p1 <- ggplot(df_tidy, aes(time, value, color = celltype))
p1 + 
  geom_line(size = 0.2) + 
  geom_point(size = 1.0) +
  theme_bw(base_size = 8) +
  xlab("time (h)") + ylab("expr. norm.") + 
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        panel.spacing.x = unit(1, "mm")) +
  facet_wrap(~gene, ncol = 5, scales = "free_y")


ggsave(sname1, width = 14, height = 7, bg = "white", dpi = 600, units = "cm")
ggsave(sname2, width = 14, height = 7, bg = "white", units = "cm")

#
