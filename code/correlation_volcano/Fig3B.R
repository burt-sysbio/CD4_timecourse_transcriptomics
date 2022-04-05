require(readr)
require(ggplot2)
require(tidyr)
require(dplyr)
require(ggrepel)
require(stringr)
source("code/utils.R")

options(ggplot2.discrete.colour = c("grey", "black", "red"))


# plot correlation volcano plots 

# assign candidates
candidates <- read.csv("genesets_input/volcano_candidates.csv", header = TRUE) 

candidates$correlation <- paste0("cor_", candidates$correlation)

# load pvalues and correlation values from masigpro run
pvals <- readRDS("data/masigpro_output/pvec_fullrun_kinbackground_a0.01")
res <- readRDS("data/masigpro_output/tfit_fullrun_kinbackground_a0.01")
res2 <- readRDS("data/masigpro_output/tfit_fullrun_kinbackgroundTH2_a0.01")
df_cor <- read.csv("data/correlation/correlation_df.csv")

# get correlation add to tidy data
df_cor <- df_cor %>% pivot_longer(-gene, values_to = "value", names_to = "correlation")

# get DEGS from masigro run
rsq <- 0.6
sigs <- get.siggenes(res, rsq = rsq, vars = "groups")
sigs <- sigs$summary
sigs2 <- get.siggenes(res2, rsq = rsq, vars = "groups")
sigs2 <- sigs2$summary


comparisons <- list("Th1_Th2", "Th1_ThMix", "Th1_Th0", "Th2_Th0", "Th2_ThMix")

# add correlation and pvalue thresholds and annotate DEGs accordingly
cor_thres <- 0.3
plots <- lapply(comparisons, prep_volcano_data, pvals, df_cor, sigs, sigs2, candidates, cor_thres)

data <- bind_rows(plots)

# split gene for left join with candidates
data <- data %>% separate(gene, into = c("candidate", "iso"), sep = "\\.", remove = F)
# remove isoforms from core volcano frame
data <- data %>% group_by(correlation, deg) %>% distinct(candidate, .keep_all = TRUE)

# add labels
data_label <- data %>% filter(deg == "qual. DEG")
data_label <- left_join(candidates, data_label, by = c("candidate", "correlation"))

# remove rows with duplicate isoform entries in each comparison
data_label <- data_label %>% group_by(correlation) %>% distinct(candidate, .keep_all = TRUE)


# plotting params
ptsize = 0.6
alpha = 0.7
labelsize = 3
yintercept = -log10(0.05)
xlabel = expression(paste("correlation index ", (1-R^{2})))
ylabel = "-log10(FDR)"

label_names <- c(
  `cor_th2_th0` = "Th2vTh0",
  `cor_th2_thmix` = "Th2vTh1/2",
  `cor_th1_th2` = "Th1vTh2",
  `cor_th1_thmix` = "Th1vTh1/2",
  `cor_th1_th0` = "Th1vTh0"
)

# plot
g <- ggplot(data, aes(cor_idx, padj_log10))
g1 <- g + 
  geom_point(data = subset(data, deg == "ns"), alpha = alpha, size = ptsize, color = "grey") +
  geom_point(data = subset(data, deg == "quant. DEG"), alpha = alpha, size = ptsize, color = "black") +
  geom_point(data = subset(data, deg == "qual. DEG"), alpha = alpha, size = ptsize, color = "red2") +
  geom_vline(xintercept=cor_thres, linetype="dashed") +
  geom_hline(yintercept=-log10(yintercept), linetype="dashed") +
  geom_text_repel(data = data_label, 
                  label = data_label$candidate,
                  size = labelsize, fontface = "italic", 
                  min.segment.length = 0,
                  max.overlaps = Inf) +
  labs(x= xlabel, y = ylabel) +
  theme_bw() +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(NA, 1.8)) +
  facet_wrap(~correlation, ncol = 5, labeller = as_labeller(label_names))

ggsave(paste0("figures/fig3B.png"), width = 9.0, height = 1.9, bg = "white")
