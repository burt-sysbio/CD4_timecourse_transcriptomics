require(readr)
require(cowplot)
require(ggplot2)
require(tidyr)
require(dplyr)
require(maSigPro)
require(ggrepel)
require(svglite)
require(stringr)

options(ggplot2.discrete.colour = c("grey", "black", "red"))

# plot volcano plots individual cell substs
# plot volcano plots with miniaml correaltion
# plot barplot for candidates polarization
################ info here which celltypes to compare #############################

get_sigs <- function(comparison, sigs, sigs2){
  if(comparison == "Th1_Th2"){
    sigs <- sigs$Th2vsTh1
    myfilter <- "cor_th1_th2"
  } else if (comparison == "Th1_ThMix"){
    sigs <- sigs$ThMixvsTh1
    myfilter <- "cor_th1_thmix"
  } else if (comparison == "Th1_Th0"){
    sigs <- sigs$Th0vsTh1
    myfilter <- "cor_th1_th0"
  } else if (comparison == "Th2_Th0"){
    sigs <- sigs2$Th0vsTh2
    myfilter <- "cor_th2_th0"
  } else if (comparison == "Th2_ThMix"){
    sigs <- sigs2$ThMixvsTh2
    myfilter <- "cor_th2_thmix"
  }
  return(list(sigs, myfilter))
}

############ end info celltypes to compare ########################################
prep_data <- function(comparison, pvals, df_cor, sigs, sigs2, candidates, cor_thres){
  
  # helper function to get comparisons right
  comp_info <- get_sigs(comparison, sigs, sigs2)
  myfilter <- comp_info[[2]]
  sigs <- comp_info[[1]]
  
  # filter data
  df_cor <- df_cor %>% filter(correlation == myfilter)
  df_pvals <- prep_pvals(pvals)
  
  # combine with pvals
  df_join <- dplyr::inner_join(df_cor, df_pvals)
  
  # add some output annotations
  out <- proc_data(df_join, sigs, candidates, cor_thres)
  return(out)
}

# process pvalues into data frame
prep_pvals <- function(pvals){
  p_adj <- pvals$p.vector
  gene_names <- rownames(pvals$p.vector)
  df_pvals <- data.frame(gene_names, p_adj)
  colnames(df_pvals) <- c("gene", "padj")
  return(df_pvals)
}

# add annotation dfs and additional readouts for colouring volcano plot
proc_data <- function(df_join, sigs, candidates, cor_thres){
  
  df_join <- df_join %>% mutate(cor_idx = 1 - value)
  
  # assign signif genes based on pval and correlation threshold
  df_join <- df_join %>% mutate(sig = ifelse((cor_idx >= cor_thres) & (padj < 0.05), "sig", "ns"))
  
  # add color based on DEG information from masigpro full run
  df_join$deg <- "ns"
  # if found in sig array, annotate
  df_join$deg[df_join$gene %in% sigs] <- "qual. DEG"
  df_join$deg[(df_join$gene %in% sigs) & (df_join$sig == "sig")] <- "quant. DEG"
  print(df_join$correlation[1])
  # convert p values to -log10
  df_join <- df_join %>% mutate(padj_log10 = -log10(padj))

  # add vector length criterion
  df_join <- df_join %>% 
    mutate(veclen1 = padj_log10 / max(padj_log10), 
           veclen2 = cor_idx / max(cor_idx),
           veclen3 = sqrt(veclen1**2 + veclen2**2))
  
  # find matches for candidates (including probe annotations with .1, .2 etc, needs str_detect or similar)
  pattern <- paste(candidates$gene , collapse = "|")
  df_annot <- df_join[str_detect(df_join$gene, pattern),]
  #df_annot <- df_join %>% filter(gene %in% candidates$gene)
  df_annot <- df_annot %>% filter(deg == "quant. DEG")

  # use the top right corner for annotation by vector length
  n_top = 10
  df_annot2 <- df_join %>% filter(deg == "quant. DEG") %>% slice_max(veclen3, n = n_top) 
  
  # just try to get the top right corner by manual thresholds
  comp <- df_join$correlation[1]

  if(comp == "th1_th0"){
    thres_cor <- 1.2
    thres_padj <- 12
  } else if(comp == "cor_th1_th0"){
    thres_cor <- 0.5
    thres_padj <- 4
  } else if(comp == "cor_th1_th2"){
    thres_cor <- 1.2
    thres_padj <- 12
  } else if(comp == "cor_th1_thmix"){
    thres_cor <- 1.2
    thres_padj <- 12
  } else if(comp == "cor_th2_th0"){
    thres_cor <- 1.2
    thres_padj <- 12
  } else if(comp == "cor_th2_thmix"){
    thres_cor <- 1.2
    thres_padj <- 12
  }
  
  df_annot3 <- df_join %>% filter((cor_idx > thres_cor) & (padj_log10 > thres_padj) & (deg == "quant. DEG"))
  df_annot4 <- bind_rows(df_annot3, df_annot) %>% unique()
  df_annot4 <- slice_sample(df_annot4, n = 8)
  
  return(list(df_join, df_annot, df_annot2, df_annot3, df_annot4))
}

save_output <- function(df_join, data_label, comparison, cor_thres){
  
  # plot individual volcano plots pval vs corr for each celltype
  # reorder data frame based on degs
  ptsize = 0.6
  alpha = 0.7
  g <- ggplot(df_join, aes(cor_idx, padj_log10))
  g1 <- g + 
    geom_point(data = subset(df_join, deg == "ns"), alpha = alpha, size = ptsize, color = "grey") +
    geom_point(data = subset(df_join, deg == "qual. DEG"), alpha = alpha, size = ptsize, color = "black") +
    geom_point(data = subset(df_join, deg == "quant. DEG"), alpha = alpha, size = ptsize, color = "red2") +
    geom_vline(xintercept=cor_thres, linetype="dashed") +
    geom_hline(yintercept=-log10(0.05), linetype="dashed") +
    geom_text_repel(data = data_label, 
                    label = data_label$gene,
                    size = 3, fontface = "italic", 
                    min.segment.length = 0,
                    max.overlaps = Inf
                    ) +
    labs(x= "correlation index", y = "-log10(FDR)") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none") +
    coord_cartesian(xlim = c(NA, 1.8))
  
  fname1 <- paste0("figures/fig3/volcano_plots/volcano_", comparison, "_cor", cor_thres, ".png")
  fname2 <- paste0("figures/fig3/volcano_plots/volcano_", comparison, "_cor", cor_thres, ".svg")
  ggsave(fname1, width =3, height = 2.5)
  ggsave(fname2, width =3, height = 2.5)
  
  # save DEGS
  degs <- df_join %>% filter(deg == "quant. DEG") %>% select(gene)
  
  # for pathway analysis get degs without isoforms
  degs2 <- df_join %>% filter(deg == "quant. DEG") %>% 
    separate(gene, c("gene", "isoform"), sep = "\\.") %>% select(gene) %>% unique()
  
  sdir <- paste0("genesets_output/quantDEGS/")
  sname <- paste0(sdir, "degs_", comparison, "_cor", cor_thres, ".csv")
  
  # additional textfile specific without isoforms
  sname2 <- paste0(sdir, "degs_cleaned_", comparison, "_cor", cor_thres, ".txt")
  
  write.csv(degs, sname, row.names = F)
  write.table(degs2, sname2, row.names = F, col.names = F, quote = F)
  
  # get genes with high correlation but not significant, high pval
  ns_genes <- df_join %>% filter((cor_idx > 0.5)&(padj > 0.05))
  sname3 <- paste0("data/nonsignif_genes/ns_genes_", comparison, "_cor", cor_thres, ".csv")
  write.csv(ns_genes, sname3, row.names = F)
  
  return(g1)
}

volcano_pipeline <- function(comparison, pvals, df_cor, sigs, sigs2, candidates, cor_thres = 0.5,
                             labeltype = "combined"){
  out <- prep_data(comparison, pvals, df_cor, sigs, sigs2, candidates, cor_thres)
  df_join <- out[[1]]
  df_annot <- out[[2]]
  df_annot2 <- out[[3]]
  df_annot3 <- out[[4]]
  df_annot4 <- out[[5]]
  
  if(labeltype == "candidates"){
    data_label <- df_annot
  } else if(labeltype == "vector"){
    data_label <- df_annot2
  } else if(labeltype == "upper_right"){
    data_label <- df_annot3
  } else if(labeltype == "combined"){
    data_label <- df_annot4
  }
  
  save_output(df_join, data_label, comparison, cor_thres)
  return(out)
}

# assign candidates
candidates <- read.csv("genesets_literature/candidate_genes_sorted_short.csv", header = TRUE) %>% select(gene)

pvals <- readRDS("data/masigpro_output/pvec_fullrun_kinbackground_a0.01")
res <- readRDS("data/masigpro_output/tfit_fullrun_kinbackground_a0.01")
res2 <- readRDS("data/masigpro_output/tfit_fullrun_kinbackgroundTH2_a0.01")
df_cor <- read.csv("data/correlation/correlation_df.csv")

# get correlation add to tidy data
df_cor <- df_cor %>% pivot_longer(-gene, values_to = "value", names_to = "correlation")

# get DEGS from masig run
rsq <- 0.6
sigs <- get.siggenes(res, rsq = rsq, vars = "groups")
sigs <- sigs$summary
sigs2 <- get.siggenes(res2, rsq = rsq, vars = "groups")
sigs2 <- sigs2$summary


comparisons <- list("Th1_Th2", "Th1_ThMix", "Th1_Th0", "Th2_Th0", "Th2_ThMix")


 
# run volcano threshold analysis for multiple thresholds
#threshold_array <- seq(0.24,0.5,0.01)

#for(cthres in threshold_array){
#  plots <- lapply(comparisons, volcano_pipeline, pvals, df_cor, sigs, sigs2, candidates, cthres)
#}

cor_thres <- 0.3
plots <- lapply(comparisons, volcano_pipeline, pvals, df_cor, sigs, sigs2, candidates, cor_thres)

#volcanos <- plot_grid(plotlist = plots, align = c("hv"), ncol = 5)

data <- bind_rows(lapply(plots, function(x) x[[1]]))
data_label <- bind_rows(lapply(plots, function(x) x[[2]]))

# remove rows with duplicate isoform entries in each comparison
data_label <- data_label %>% separate(gene, into = c("gene", "ISO"), sep = "\\.")
data_label <- data_label %>% group_by(correlation) %>% distinct(gene, .keep_all = TRUE)

label_names <- c(
  `cor_th2_th0` = "Th2vTh0",
  `cor_th2_thmix` = "Th2vTh1/2",
  `cor_th1_th2` = "Th1vTh2",
  `cor_th1_thmix` = "Th1vTh1/2",
  `cor_th1_th0` = "Th1vTh0"
)

ptsize = 0.6
alpha = 0.7
g <- ggplot(data, aes(cor_idx, padj_log10))
g1 <- g + 
  geom_point(data = subset(data, deg == "ns"), alpha = alpha, size = ptsize, color = "grey") +
  geom_point(data = subset(data, deg == "qual. DEG"), alpha = alpha, size = ptsize, color = "black") +
  geom_point(data = subset(data, deg == "quant. DEG"), alpha = alpha, size = ptsize, color = "red2") +
  geom_vline(xintercept=cor_thres, linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_text_repel(data = data_label, 
                  label = data_label$gene,
                  size = 3, fontface = "italic", 
                  min.segment.length = 0,
                  max.overlaps = Inf
  ) +
  labs(x= "correlation index", y = "-log10(FDR)") +
  theme_bw() +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(NA, 1.8)) +
  facet_wrap(~correlation, ncol = 5, labeller = as_labeller(label_names))

ggsave(paste0("figures/fig3/volcano_plots/volcanos_cor", cor_thres, ".png"), width = 9.0, height = 1.9, bg = "white")
