require(dplyr)
require(tidyr)
require(readr)
require(pheatmap)
source("code/paper_theme.R")


plot_timecourse <- function(baseline){
  modeltype <- "Ind"
  
  # load data and add minimal of both pvalues for sorting
  mygenes <- read.csv(paste0("data/linear_model/linear_model_processed_", baseline, ".csv"))
  pvals <- mygenes[c("p_th1", "p_th2")]
  pvals$pnew <- apply(pvals, 1, min) 
  mygenes <- left_join(mygenes, pvals)
  mygenes <- mygenes %>% filter(category == modeltype, model == "Th1/2") %>%
    select(gene, category, model, pnew)
  
  
  # only focus on first 30 genes or so
  n <- 40
  mygenes <- mygenes %>% filter(gene != "Nrp1.1") %>% filter(gene != "H2-D1.1") %>% 
    filter(gene != "Oas1g.3") %>% filter(gene != "Oasl2.1")
  mygenes <- mygenes %>% slice_max(pnew, n = n)
  
  
  # load timecourse data and sort by pnew
  df <- read_csv("data/peine_logfc_tidy.csv")
  
  df <- inner_join(df, mygenes)
  
  # only keep genes that have fold-change > n
  #df2 <- df %>% filter(celltype=="Th1/2") %>% filter(abs(value) > 0.5)
  
  #df <- df[df$gene %in% df2$gene,]
  
  
  df$gene <- reorder(df$gene, -df$pnew)
  
  
  
  
  df <- df %>% separate(gene, into = c("gene", "iso"), sep = "\\.")
  
  g <- ggplot(data = df, aes(x=time, y=value, colour = celltype))
  g2 <- g + geom_point() + geom_line() + theme_bw(base_size = 10) + facet_wrap(~gene, scales = "free_y", ncol = 6) +
    ylab("logFC(day0") + xlab("time (h)") 
  
  sname <-paste0("figures/fig4/linear_model/timecourses/timecourse_linear_model_", baseline, "_", modeltype, ".png")
  ggsave(sname, bg = "white", width = 12, height = 8)
  
}
