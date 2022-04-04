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
  df_join$deg[df_join$gene %in% sigs] <- "quant. DEG"
  df_join$deg[(df_join$gene %in% sigs) & (df_join$sig == "sig")] <- "qual. DEG"

  # convert p values to -log10
  df_join <- df_join %>% mutate(padj_log10 = -log10(padj))
  
  return(df_join)
}

# assign candidates
candidates <- read.csv("genesets_literature/volcano_candidates.csv", header = TRUE) 

candidates$correlation <- paste0("cor_", candidates$correlation)


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


cor_thres <- 0.3
plots <- lapply(comparisons, prep_data, pvals, df_cor, sigs, sigs2, candidates, cor_thres)

data <- bind_rows(plots)

# split gene for left join with candidates
data <- data %>% separate(gene, into = c("candidate", "iso"), sep = "\\.", remove = F)
# remove isoforms from core volcano frame
data <- data %>% group_by(correlation, deg) %>% distinct(candidate, .keep_all = TRUE)

data_label <- data %>% filter(deg == "qual. DEG")
data_label <- left_join(candidates, data_label, by = c("candidate", "correlation"))

# remove rows with duplicate isoform entries in each comparison
data_label <- data_label %>% group_by(correlation) %>% distinct(candidate, .keep_all = TRUE)

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
  geom_point(data = subset(data, deg == "quant. DEG"), alpha = alpha, size = ptsize, color = "black") +
  geom_point(data = subset(data, deg == "qual. DEG"), alpha = alpha, size = ptsize, color = "red2") +
  geom_vline(xintercept=cor_thres, linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_text_repel(data = data_label, 
                  label = data_label$candidate,
                  size = 3, fontface = "italic", 
                  min.segment.length = 0,
                  max.overlaps = Inf) +
  labs(x= expression(paste("correlation index ", (1-R^{2}))), y = "-log10(FDR)") +
  theme_bw() +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(NA, 1.8)) +
  facet_wrap(~correlation, ncol = 5, labeller = as_labeller(label_names))

ggsave(paste0("figures/fig3/volcano_plots/volcanos_new_annot.png"), width = 9.0, height = 1.9, bg = "white")
