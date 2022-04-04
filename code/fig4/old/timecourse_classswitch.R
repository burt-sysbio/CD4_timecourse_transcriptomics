
# timecourse for class switches
require(ggplot2)
require(dplyr)
require(tidyr)
df2 <- readRDS("r_objects/fig3/cluster_assignment.RDS")
df_cor <- readRDS("r_objects/correlation_df.RDS")

# check for class switches
df2 <- df2 %>% pivot_wider(names_from = celltype, values_from = annot)

#
df2$s_th1_th2 <- df2$Th1 != df2$Th2
df2$s_th1_thm <- df2$Th1 != df2$ThMix
df2$s_th2_thm <- df2$Th2 != df2$ThMix

df2$s_th1_th0 <- df2$Th1 != df2$Th0
df2$s_th2_th0 <- df2$Th2 != df2$Th0
df2$s_th0_thm <- df2$Th0 != df2$ThMix

# kick out genes for which there are no class switches
df3 <- df2[,6:ncol(df2)]
df4 <- df2[rowSums(df3,na.rm = T) > 0,]

# assign candidates
df_cor <- inner_join(df4, df_cor)

mydf <- df_cor[,2:5]

# for each row, check if any cluster is 4(down reg)
down_switch <-  apply(mydf, 1, function(x) (any(x == 4)))
# for each row check if any cluster is 2 while another cluster is 3
trans_switch <- apply(mydf, 1, function(x) (any(x == 3) & any(x== 2)))


df_trans <- df_cor[which(trans_switch == T),]
df_down <-  df_cor[which(down_switch == T),]


# change this between df trans and df down
# load tidy data
df_tidy <- readRDS("data/peine_logfc_tidy.RDS")

# if you wan tto check candidates use this, otherwise df_down/trans$Symbol for filtering df_tidy to get overview
mygenes <- c("Il8rb", "Slfn1")
df_tidy <- df_tidy %>% filter((gene %in% mygenes))
df_tidy$time <- as.numeric(df_tidy$time)

p1 <- ggplot(df_tidy, aes(time, value, color = celltype))
p1 + geom_line(size = 1) + 
  theme_bw() + 
  facet_wrap(~gene, scales = "free") +
  xlab("time (h)") +
  ylab("logFC(expr.)") +
  theme(legend.position = "none")
#ggsave("figures/fig3/class_switch_candidates.png", width = 3.5, height = 2)
#
