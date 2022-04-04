
# spread cor1 and cor2, need to kick out cols that are different between th1 and th2
# except for corr col to pivot

# only keep genes that are significant (rsq and not all correlated)
# note that sig col is different from sig col in df join

require(ggplot2)
require(readr)
require(dplyr)
require(tidyr)
require(cowplot)

rsq_thres <- 0.6
loadname <- paste0("data/minimal_correlation_rsq", rsq_thres, ".csv")
mincor <- read.csv(loadname)
mincor <- mincor %>% select(-X) %>% filter(sig == "sig")

# kick out Th1_Th2 corr

mincor <- mincor %>% filter(cell != "Th1_Th2")

df2 <- mincor %>% select(gene, cell, cor_idx) %>% pivot_wider(names_from = cell, values_from = cor_idx)
df2$polarization_thm <- log(df2$Th1_ThMix/df2$Th2_ThMix, base = 2)
df2$polarization_th0 <- log(df2$Th1_Th0/df2$Th2_Th0, base = 2)

myfun <- function(df, contrast){
  if(contrast == "Th0"){
    g <- ggplot(data = df2, aes(x=reorder(gene, polarization_th0), y= polarization_th0))
  } else if(contrast == "Thm"){
    g <- ggplot(data = df2, aes(x=reorder(gene, polarization_thm), y= polarization_thm))
  }
  g1 <- g +
    geom_bar(stat = "identity", fill = "grey") +
    theme_half_open(font_size = 20)+
    xlab("genes ranked") +
    ylab("log2(cor.idx Th1 / cor.idx Th2)") +
    geom_hline(yintercept = 1) + 
    geom_hline(yintercept = -1) +
    scale_x_discrete(labels = NULL, breaks = NULL) +
    coord_flip(ylim = c(-5, 5))
  sname <- paste0("figures/fig3/polarization/polarization_nofilter_rsq", rsq_thres, contrast,".png")
  ggsave(sname, width =4, height = 5)
}

myfun(df2, "Th0")
myfun(df2, "Thm")

# summarize which genes are above or below threshold
df3 <- df2 %>% pivot_longer(cols = starts_with("polarization_"))
df3$category <- "superpos."
df3$category[df3$value < -1] <- "Th1-like"
df3$category[df3$value > 1] <- "Th2-like"

df4 <- df3 %>% group_by(name, category) %>% summarize(n = n())

g <- ggplot(data = df4, aes(category, n, fill = name))
g + geom_bar(stat = "identity", position = position_dodge2()) +
  theme_half_open()+
  xlab("") +
  ylab("n genes") +
  scale_y_continuous( expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c("grey", "purple")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
sname <- paste0("figures/fig3/polarization/polarization_barplot_rsq", rsq_thres,".png")
ggsave(sname, width =4, height = 3)

# polarization for signature genes
mygenes <- read_csv("genes_caro.csv", col_names = "gene")
df_mygenes <- df3 %>% filter(gene %in% mygenes$gene)

g <- ggplot(data = df_mygenes, aes(x= value, y = gene, fill = name))
g1 <- g+
  geom_bar(stat= "identity", alpha = 0.9)+
  theme_half_open(font_size = 10)+
  xlab("log2(cor.idx Th1 / cor.idx Th2)") +
  ylab("") +
  coord_cartesian(xlim= c(-12,12)) +
  scale_fill_manual(values = c("grey", "purple")) +
  scale_x_continuous(breaks = c(-10,-5,0,5,10))

sname <- paste0("figures/fig3/polarization/polarization_caro_rsq", rsq_thres,".png")
ggsave(sname, width = 4, height = 3)

# 
# plot_polarization <- function(df, contrast, candidates, sname){
#   if(contrast == "ThMix"){
#     myvec <- c("Th1_ThMix", "Th2_ThMix")
#   }else{
#     myvec <- c("Th1_Th0", "Th2_Th0")
#   } 
#   df_barplot <- df %>% filter(cell %in% myvec)
#   
#   df_barplot2 <- df_barplot %>% select(gene, cell, cor_idx) %>%
#     pivot_wider(names_from = cell, values_from = cor_idx)
#   
#   df_barplot2$ratio <- df_barplot2[[2]] / df_barplot2[[3]]
#   df_barplot2$pol1 <- log2(df_barplot2$ratio)
#   
#   df_barplot2 <- na.omit(df_barplot2)
#   
#   
#   if(!is.null(candidates)){
#     df_barplot3 <- df_barplot2 %>% filter(gene %in% candidates$gene)
#     df_barplot3 <- na.omit(df_barplot3)
#     g <- ggplot(data = df_barplot3, aes(x = reorder(gene, pol1), y= pol1))
#     g2 <- g + geom_bar(stat= "identity") + 
#       coord_flip() +
#       theme_half_open(font_size = 25)+
#       xlab("") +
#       ylab("Polarization") 
#     
#     ggsave(paste0("figures/fig3/polarization/barplot_polarization_", sname, "_", contrast, ".png"), 
#            width = 5, height = 8)
#     
#     print(g2)
#     return(list(df_barplot3, g2))
#     
#   } else{
#     return(list(df_barplot2, "none"))
#   }
# }
# 
# plot_histograms <- function(df, sname, candidates = NULL){
#   
#   # get polarization data for thm and th0 and bind together
#   out1 <- plot_polarization(df, "ThMix", candidates, sname)
#   out2 <- plot_polarization(df, "Th0", candidates, sname)
#   df1 <- out1[[1]] %>% select(gene, pol1)
#   df2 <- out2[[1]] %>% select(gene, pol1)
#   
#   df1$cell <- "ThMix"
#   df2$cell <- "Th0"
#   
#   
#   df <- bind_rows(df1, df2)
#   df$cell <- as.factor(df$cell)
#   g <- ggplot(data = df, aes(pol1, fill = cell))
#   g2 <- g + geom_histogram(binwidth = 0.5, alpha = 0.5, position = "identity") + theme_minimal()
#   print(g2)
# }
# 
# c2 <- read_csv("pathway_analysis/gene_sets/references/stubbington_th1.csv")
# c3 <- read_csv("pathway_analysis/gene_sets/references/stubbington_th2.csv", col_names =  "gene")
# c4 <- read_csv("pathway_analysis/gene_sets/references/stubbington_tfs.csv")
# c5 <- read_csv("pathway_analysis/gene_sets/references/stubbington_cytokines.csv")
# c6 <- read_csv("pathway_analysis/gene_sets/references/stubbington_receptors.csv")
# c7 <- read_csv("pathway_analysis/gene_sets/references/wei_gata_targets_th1.csv", col_names =  "gene")
# c8 <- read_csv("pathway_analysis/gene_sets/references/wei_gata_targets_th2.csv", col_names =  "gene")
# c9 <- read_csv("pathway_analysis/gene_sets/references/zhu_tbet_targets_neg.csv", col_names =  "gene")
# c10 <- read_csv("pathway_analysis/gene_sets/references/zhu_tbet_targets_pos.csv", col_names =  "gene")
# colnames(c2) <- "gene"
# colnames(c3) <- "gene"
# colnames(c4) <- "gene"
# colnames(c5) <- "gene"
# colnames(c6) <- "gene"
# 
# sname <- "cytokines"
# candidates <- c5
# 
# candidate_list <- c(c2,c3,c4,c5,c6,c7,c8,c9,c10)
# names(candidate_list)<- c("th1", "th2", "tf", "cytos", "receptors", "gata_target_th1", "gata_target_th2",
#                           "tbet_target_neg", "tbet_target_pos")

#df <- plot_histograms(mincor, "test")
#g <- plot_polarization(mincor, contrast = "ThMix", candidates = c10, sname = "tbet2")

#candidates <- c5
#name <- "cytos"
#plot_polarization(mincor, contrast = "Th0", candidates = candidates, sname = name)
#plot_polarization(mincor, contrast = "ThMix", candidates = candidates, sname = name)



#candidates <- bind_rows(c7, c8,c9,c10)