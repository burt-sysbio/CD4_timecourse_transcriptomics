require(readr)
require(ggplot2)
require(tidyr)
require(dplyr)
require(maSigPro)
require(cowplot)

options(ggplot2.discrete.color = c("grey", "red"))
options(ggplot2.discrete.fill = c("grey", "red"))

# plot volcano plots individual cell substs
# plot volcano plots with miniaml correaltion
# plot barplot for candidates polarization

savedir <- "figures/fig3/"

# get correlation and tidy up
df_tidy <- read.csv("data/peine_logfc_tidy.csv")
df_cor <- read.csv("data/correlation_df.csv")


df_cor <- df_cor %>% select(gene, cor_th1_th2)
df_tidy <- left_join(df_tidy, df_cor)


# get pvals from masigpro run, use both thm and th0 as contrast
res <- readRDS("code/fig2/mspro_tstep_contrast_Th1.RDS")
res2 <- readRDS("code/fig2/mspro_tstep_contrast_Th2.RDS")

sigs2 <- get.siggenes(res2, rsq = 0, vars = "groups")

rsq_thres <- 0.7

out <- sigs$sig.genes$Th2vsTh1$sig.pvalues
out$gene <- rownames(out)

out2 <- sigs2$sig.genes$Th1vsTh2$sig.pvalues
out2$gene <- rownames(out2)


# only use th1 th2 degs as background!
rsq_vals <- out %>% select("gene", "p-value", "R-squared")
colnames(rsq_vals) <- c("gene", "p_adj", "rsq")
df_join <- dplyr::inner_join(df_tidy, rsq_vals)
df_join <- df_join %>% mutate(cor_idx = 1 - cor)

# assign signif genes based on rsq and correlation threshold
cor_thres <- 0.3
df_join <- df_join %>% 
  mutate(sig = ifelse((rsq >= rsq_thres)&(cor_idx >= cor_thres), "sig", "ns"))

# assign candidates

candidates <- read_csv("genes_caro.csv", col_names = "gene")

df_annot <- df_join %>% filter(gene %in% candidates$gene)
df_annot <- df_annot %>% filter(rsq >= rsq_thres, cor_idx >= cor_thres)

# plot individual volcano plots rsq vs corr for each celltype (for Th1Th2 Degs)
g <- ggplot(df_join, aes(rsq, cor_idx))
g + geom_point(alpha = 0.8, size = 0.5, aes(color = sig)) +
  geom_hline(yintercept=cor_thres, linetype="dashed") +
  geom_vline(xintercept=rsq_thres, linetype="dashed") +
  geom_text(data = df_annot, label = df_annot$gene, check_overlap = T) +
  facet_wrap(~cell, nrow = 1, ncol = 5) +
  coord_cartesian(xlim = c(0.2, NA)) +
  ylab("correlation index") +
  xlab(expression(R^{2})) +
  theme_bw(base_size = 15) +
  theme(legend.position = "none") 

ggsave("figures/fig4/volcano_plots.png", width =9, height = 2.5)


# get minimal correlation across cell comparisons and plot one volcano plot for all
# this subset is used for downstream clustering and identifying unique genes etc
mincor <- df_join %>% group_by(gene) %>% summarize(mincor = min(cor))
mincor <- mincor %>% mutate(mincor_idx = 1- mincor)
mincor <- dplyr::left_join(df_join, mincor)

mincor <- mincor %>% 
  mutate(sig = ifelse((rsq >=rsq_thres)&(mincor_idx >= cor_thres), "sig", "ns"))


# 
sname <- paste0("data/minimal_correlation_rsq", rsq_thres, ".csv")
write.csv(mincor, sname, row.names = F)

df_annot2 <- mincor %>% filter(gene %in% candidates$gene)
df_annot2 <- df_annot2 %>% filter(rsq >=rsq_thres, mincor_idx >= cor_thres)


g <- ggplot(mincor, aes(rsq, mincor_idx, )) +
  geom_point(alpha = 0.5, size = 0.1, aes(color = sig)) +
  geom_hline(yintercept=cor_thres, linetype="dashed") +
  geom_vline(xintercept=rsq_thres, linetype="dashed") +
  geom_text(data = df_annot2, label = df_annot2$gene, check_overlap = T) +
  coord_cartesian(xlim = c(0.2, NA)) +
  ylab("cor. idx min.") +
  xlab(expression(R^{2})) +
  theme_bw(base_size = 15) +
  theme(legend.position = "none") 

ggsave("figures/fig4/volcano_plot_minimal_correlation.png", width =5, height = 4)

