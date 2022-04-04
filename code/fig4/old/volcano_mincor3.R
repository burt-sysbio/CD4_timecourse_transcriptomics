require(readr)
require(ggplot2)
require(tidyr)
require(dplyr)
require(maSigPro)
require(cowplot)

options(ggplot2.discrete.color = c("grey", "darkgreen", "red"))
options(ggplot2.discrete.fill = c("grey", "darkgreen", "red"))


# plot volcano plots individual cell substs
# plot volcano plots with miniaml correaltion
# plot barplot for candidates polarization

savedir <- "figures/fig4/"

# get correlation add to tidy data
df_cor <- read.csv("data/correlation_df.csv")
df_cor <- df_cor %>% pivot_longer(-gene, values_to = "value", names_to = "correlation")



# get DEGS from masig run
res <- readRDS("data/masigpro_output/tfit_fullrun_kinbackground_a0.01")
sigs <- get.siggenes(res, rsq = 0.6, vars = "groups")
sigs <- sigs$summary


################ info here which celltypes to compare #############################
comparison <- "Th1_Th0"

if(comparison == "Th1_Th2"){
  sigs <- sigs$Th2vsTh1
  myfilter <- "cor_th1_th2"
} else if(comparison == "Th1_ThMix") {
  sigs <- sigs$ThMixvsTh1
  myfilter <- "cor_th1_thmix"
} else if(comparison == "Th1_Th0"){
  sigs <- sigs$Th0vsTh1
  myfilter <- "cor_th1_th0"
}

############ end info celltypes to compare ########################################

df_cor <- df_cor %>% filter(correlation == myfilter)

# get pvals from masigpro run
pvals <- readRDS("data/masigpro_output/pvec_fullrun_kinbackground_a0.01")
p_adj <- pvals$p.vector
gene_names <- rownames(pvals$p.vector)
df_pvals <- data.frame(gene_names, p_adj)
colnames(df_pvals) <- c("gene", "padj")

# combine correlation and pvals
df_join <- dplyr::inner_join(df_cor, df_pvals)
df_join <- df_join %>% mutate(cor_idx = 1 - value)

# assign signif genes based on pval and correlation threshold
cor_thres <- 0.15
df_join <- df_join %>% 
  mutate(sig = ifelse((cor_idx >= cor_thres) & (padj < 0.05), "sig", "ns"))


# add color based on DEG information from masigpro full run
df_join$deg <- "ns"
df_join$deg[df_join$gene %in% sigs] <- "qual. DEG"
df_join$deg[(df_join$gene %in% sigs) & (df_join$sig == "sig")] <- "quant. DEG"


# convert p values to -log10
df_join <- df_join %>% mutate(padj_log10 = -log10(padj))

# assign candidates
candidates <- read_csv("genes_caro.csv", col_names = "gene")
df_annot <- df_join %>% filter(gene %in% candidates$gene)
df_annot <- df_annot %>% filter(padj < 0.05, cor_idx >= cor_thres)

# plot individual volcano plots pval vs corr for each celltype
g <- ggplot(df_join, aes(cor_idx, padj_log10))
g + geom_point(alpha = 0.5, size = 1.0, aes(color = deg)) +
  geom_vline(xintercept=cor_thres, linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_text(data = df_annot, label = df_annot$gene, check_overlap = T) +
  #coord_cartesian(xlim = c(1.2, 10)) +
  xlab("correlation index") +
  ylab("-log10(p-adj)") +
  theme_bw(base_size = 12) +
  theme(legend.title=element_blank()) 

fname <- paste0("figures/fig4/volcano_plot_", comparison, ".png")
ggsave(fname, width =4, height = 2.5)

# get DEGS
degs <- df_join %>% filter(deg == "quant. DEG") %>% select(gene)

