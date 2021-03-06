require(dplyr)
require(tidyr)
require(ggplot2)
source("code/paper_theme.R")
require(RColorBrewer)



plot_pathways <- function(baseline, approach = "linear model", remove_category = c("Tbet_pos", "Tbet_neg")){
  df <- read.csv(paste0("code/linear_model/chipseq_pathway_results/chipseq_pathway_results_", baseline, ".csv"))
  
  # rename labels
  df$ID[df$ID=="Tbet_large"] <- "Tbet"
  df <- df[df$approach == approach,]
  
  # kick out unwanted categories
  df <- df[!(df$ID %in% remove_category),]
  
  # transform pvalues
  df <- df %>% mutate(pval2 = -log10(pvalue))
  
  # plotting params
  pal <- brewer.pal(name ="Dark2", n = 8)
  mycols <- c(pal[2], pal[1], pal[6], pal[8])
  
  # plot
  g <- ggplot(data = df, aes(x = ID, y = pval2))
  
  g2 <-  g + geom_bar(aes(fill=model), stat = "identity", position = position_dodge()) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_y_continuous(limits = c(0,NA), expand = expansion(mult = c(0, .1))) +
    geom_hline(yintercept = -log10(0.05))+
    labs(x = "", y = "-log10(pvalue)") +
    scale_fill_manual(values=mycols)
  
  ggsave("figures/Fig4C_2.png", width = 3, height = 2, bg = "white")
  
}

baseline <- "degs_Th1_Th2_cor0.3"

plot_pathways(baseline)