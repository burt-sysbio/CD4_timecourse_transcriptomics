require(readr)
require(ggplot2)
require(tidyr)
require(dplyr)
require(maSigPro)
require(ggthemes)
require(ggpubr)
options(ggplot2.discrete.color = c("grey", "red"))
options(ggplot2.discrete.fill = c("grey", "red"))


savedir <- "figures/fig2/"

# get pvals from masigpro run, use both thm and th0 as contrast
tstep_th0 <- readRDS("r_objects/fig2/msigpro_Th1_Th2_Th0_Q0.05a0.05.RDS")
tstep_thm <- readRDS("r_objects/fig2/msigpro_Th1_Th2_ThMix_Q0.05a0.05.RDS")

sigs_thm <- get.siggenes(tstep_thm, vars = "groups")
sigs_th0 <- get.siggenes(tstep_th0, vars = "groups")

p1 <- sigs_th0$sig.genes$Th1vsTh0$sig.pvalues
p2 <- sigs_th0$sig.genes$Th2vsTh0$sig.pvalues
p3 <- sigs_thm$sig.genes$Th1vsThMix$sig.pvalues
p4 <- sigs_thm$sig.genes$Th2vsThMix$sig.pvalues


p1$id <- rownames(p1)
p1$cell <- "Th1_Th0"

p2$id <- rownames(p2)
p2$cell <- "Th2_Th0"

p3$id <- rownames(p3)
p3$cell <- "Th1_ThMix"

p4$id <- rownames(p4)
p4$cell <- "Th2_ThMix"

p <- list(p1,p2,p3,p4)

p <- lapply(p, function(x) x[c("p-value", "id", "cell")])
p <- do.call("rbind", p)
p <- p %>% rename(p_val = 'p-value')


gene_names <- readRDS("data/peine_featuredata.RDS")

# get correlation and tidy up
df_cor <- readRDS("r_objects/fig2/correlation_df.RDS")

tidy_cor <- function(df_cor){
  df_cor <- df_cor[, c("Symbol", "id", "cor_th1_thm", "cor_th2_thm", "cor_th1_th0", "cor_th2_th0")]
  colnames(df_cor) <- c("Symbol", "id", "Th1_ThMix", "Th2_ThMix", "Th1_Th0", "Th2_Th0")
  
  df_cor <- df_cor %>% pivot_longer(-c(Symbol,id), names_to = "cell", values_to = "cor")
  return(df_cor)
}


# combine p values with correlation
df_tidy <- tidy_cor(df_cor)
df_join <- dplyr::inner_join(df_tidy, p)
df_join <- df_join %>% mutate(log10padj = -log10(p_val), cor_idx = 1-cor)



# assign candidates
my_tf <- c("Id3", "Lef1", "Tcf7", "Klf4", "Gata3", "Tbx21", "Stat1", "Stat2", "Stat4", "Irf8", "Id2",
           "Stat3", "Klf6", "Maf", "Prdm1", "Foxp3", "Klf2")
my_cyto <- c("Il4", "IFNg", "Il2", "Il6st", "Il2ra", "Il10", "Il12rb2", "Icos", "Tnfsf8")
my_rec <- c("Cxcr3, Cxcr4, Cxcr5", "Cxcr6", "Ccr7", "Ccr2", "Ccl3")
candidates <- c(my_tf, my_cyto, my_rec)

# add candidates
cor_thres <- 0.3
df_annot <- df_join %>% filter(Symbol %in% candidates)
df_annot <- df_annot %>% filter(p_val <= 0.05, cor_idx >= cor_thres)

df_join <- df_join %>% 
  mutate(sig = ifelse((p_val <=0.05)&(cor_idx >= cor_thres), "sig", "ns"))


g <- ggplot(df_join, aes(log10padj, cor_idx))
g + geom_point(alpha = 0.8, size = 0.5, aes(color = sig)) +
  geom_hline(yintercept=cor_thres, linetype="dashed") +
  geom_text(data = df_annot, label = df_annot$Symbol, check_overlap = T) +
  facet_wrap(~cell, nrow = 2, ncol = 2) +
  coord_cartesian(xlim = c(0, NA)) +
  ylab("correlation index") +
  xlab("-log10(padj)") +
  theme_bw() +
  theme(legend.position = "none") 

ggsave("figures/fig2/volcano2.png", width = 5, height = 4)


# save output of df_cor for hits with adj p and cor thres
genes_out <- df_join %>% filter(sig == "sig")
genes_out <- genes_out %>% select(Symbol, cell, id)

saveRDS(genes_out, "r_objects/fig2/degs_cor_masigpro.RDS")