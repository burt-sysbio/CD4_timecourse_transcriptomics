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

savedir <- "figures/fig4/"

compare <- "Th2vsTh1"
# get correlation add to tidy data
df_cor <- read.csv("data/correlation_df.csv")
df_cor <- df_cor %>% pivot_longer(-gene, values_to = "value", names_to = "correlation")


# get pvals from masigpro run, use both thm and th0 as contrast
res <- readRDS("data/masigpro_output/tfit_fullrun_kinbackground_a0.01")
pvals <- readRDS("data/masigpro_output/pvec_fullrun_kinbackground_a0.01")
p_adj <- pvals$p.adjusted
gene_names <- rownames(res$sol)


sigs <- get.siggenes(res, rsq = 0.6, vars = "groups")
sigs <- sigs$summary



rsq_thres <- 0.7

if (compare == "Th2vsTh1"){
  out <- sigs$sig.genes$Th2vsTh1$sig.pvalues
  df_cor <- df_cor %>% filter(correlation == "cor_th1_th2")
  rsq_vals <- data.frame(out$p.valor_timePointxTh2, out$p.valor_timePoint2xTh2,
                         out$p.valor_timePoint3xTh2, out$p.valor_timePoint4xTh2)
  
}

rsq_vals$pmax <- apply(rsq_vals, 1, max, na.rm = T)
rsq_vals$gene <- rownames(out)
rsq_vals$padj_anova <- out$`p-value`
rsq_vals$rsq <- out$`R-squared`

df_join <- dplyr::inner_join(df_cor, rsq_vals)
df_join <- df_join %>% mutate(cor_idx = 1 - value)

# assign signif genes based on rsq and correlation threshold
cor_thres <- 0.3
df_join <- df_join %>% 
  mutate(sig = ifelse((cor_idx >= cor_thres), "sig", "ns"))

# convert p values to -log10
df_join <- df_join %>% mutate(padj_log10 = -log10(pmax))
# assign candidates

candidates <- read_csv("genes_caro.csv", col_names = "gene")

df_annot <- df_join %>% filter(gene %in% candidates$gene)
df_annot <- df_annot %>% filter(rsq >= rsq_thres, cor_idx >= cor_thres)

# plot individual volcano plots rsq vs corr for each celltype (for Th1Th2 Degs)
g <- ggplot(df_join, aes(padj_log10, cor_idx))
g + geom_point(alpha = 0.8, size = 0.5, aes(color = sig)) +
  geom_hline(yintercept=cor_thres, linetype="dashed") +
  geom_text(data = df_annot, label = df_annot$gene, check_overlap = T) +
  coord_cartesian(xlim = c(1.2, 10)) +
  ylab("correlation index") +
  xlab("-log10(p-adj)") +
  theme_bw(base_size = 15) +
  theme(legend.position = "none") 

fname <- paste0("figures/fig4/volcano_plot", compare, ".png")
ggsave(fname, width =3, height = 2.5)

