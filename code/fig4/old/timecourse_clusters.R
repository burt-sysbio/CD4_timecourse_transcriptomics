require(ggplot2)
source("code/paper_theme.R")
# use color scheme wo th0
set_colors(F)

df1 <- readRDS("r_objects/fig3/timecourse_clusterTh1.RDS")
df2 <- readRDS("r_objects/fig3/timecourse_clusterTh2.RDS")
df3 <- readRDS("r_objects/fig3/timecourse_clusterThMix.RDS")
df4 <- readRDS("r_objects/fig3/timecourse_clusterTh0.RDS")

n_th1 <- nrow(df1)
n_th2 <- nrow(df2)
n_thm <- nrow(df3)
n_th0 <- nrow(df4)

df_cor <- readRDS("r_objects/correlation_df.RDS")

# reordering of clusters
ann_th2 <- data.frame(annot = seq(1:4), ann_new = c(1,4,2,3)) 
ann_thm <- data.frame(annot = seq(1:4), ann_new = c(1,4,3,2)) 
ann_th0 <- data.frame(annot = seq(1:4), ann_new = c(1,4,2,3)) 


ann_th2$annot <- as.factor(ann_th2$annot)
ann_th2$ann_new <- as.factor(ann_th2$ann_new)

ann_thm$annot <- as.factor(ann_thm$annot)
ann_thm$ann_new <- as.factor(ann_thm$ann_new)


ann_th0$annot <- as.factor(ann_th0$annot)
ann_th0$ann_new <- as.factor(ann_th0$ann_new)


df2 <- left_join(df2, ann_th2) %>% select(-annot) 
df2 <- rename(df2, annot = ann_new)

df3 <- left_join(df3, ann_thm) %>% select(-annot) 
df3 <- rename(df3, annot = ann_new)


df4 <- left_join(df4, ann_th0) %>% select(-annot) 
df4 <- rename(df4, annot = ann_new)



df1$celltype <- "Th1"
df2$celltype <- "Th2"
df3$celltype <- "ThMix"
df4$celltype <- "Th0"

df <- rbind(df1,df2,df3, df4)
annot <- df_cor[c("Symbol", "id")]
df <- left_join(df, annot)
df <- df %>% select(-id)

# plot time course of averages per cluster
p1 <- ggplot(df, aes(time, avg, color = celltype))
p1 + geom_line(size = 1) + 
  facet_wrap(~annot) +
  theme_bw(base_size = 15) +
  ylab("expr fc") + xlab("time (h)")
#ggsave("figures/fig3/timecourse.png", width = 5, height = 3)
#annot$id <- rownames(annot)

# plot gebes per cluster
df2 <- df %>% select(Symbol, annot, celltype)
df2 <- unique(df2)

# summarize number of genes in barplot
df_bar <- df2 %>% group_by(celltype, annot) %>% count() %>% ungroup()
df_bar <- df_bar %>% group_by(celltype) %>% mutate(ngenes = sum(n)) %>% ungroup()
df_bar <- df_bar %>% mutate(n_perc = 100*(n/ngenes))


g <- ggplot(df_bar, aes(annot, n_perc, fill = celltype))
g + 
  geom_bar(position = position_dodge(), width = 0.8, stat = "identity") + 
  theme_bw(base_size = 15) +
  xlab("cluster") + ylab("genes in cluster (%)")

# make th1 and th2 barplot but also add information about shared genes
df_t12 <- df2 %>% filter(celltype %in% c("Th1", "Th2")) %>% pivot_wider(names_from = celltype, values_from = annot)
df_t12$new <- "switch"
df_t12$new[df_t12$Th1 == df_t12$Th2] <- "shared"

df_t12$new[(is.na(df_t12$Th1)) & (!(is.na(df_t12$Th2)))] <- "th2_u"
df_t12$new[(is.na(df_t12$Th2)) & (!(is.na(df_t12$Th1)))] <- "th1_u"

df_t12 <- df_t12 %>% pivot_longer(-c(Symbol, new)) %>% na.omit()
df_shared <- df_t12 %>% filter(new == "shared")

g <- ggplot(df_t12, aes(value, fill = name))
g+ geom_bar(position = position_dodge()) + mytheme() +
  geom_bar(data = df_shared, color = "black", aes(value, fill = name), position = position_dodge())

#ggsave("figures/fig3/barplot.png", width = 4, height = 3)

saveRDS(df2, "r_objects/fig3/cluster_assignment.RDS")