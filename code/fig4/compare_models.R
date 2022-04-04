require(dplyr)
require(tidyr)

dir1 = "data/correlation/hclust_output/hclust_corr_thm.csv"
dir2 = "data/correlation/hclust_output/hclust_corr_th0.csv"

baseline = "degs_Th1_Th2_cor0.3"
dir3 = paste0("data/linear_model/linear_model_processed_" , baseline , ".csv")


hclust_thm <- read.csv(dir1)
hclust_th0 <- read.csv(dir2)

df3 <- read.csv(dir3)

linmo_th0 <- df3 %>% filter(model == "Th0")
linmo_thm <- df3 %>% filter(model == "Th1/2")


myfun <- function(category, hclust_df, linmo_df){
  hclust_genes <- hclust_df$gene[hclust_df$category == category]
  linmo_genes <- linmo_df$gene[linmo_df$category == category]
  
  overlap <- intersect(hclust_genes, linmo_genes)
  combined <- union(hclust_genes, linmo_genes)
  out <- length(overlap) / length(combined)
  return(out)
}

test <- myfun("Superpos", hclust_thm, linmo_thm)
print(test)