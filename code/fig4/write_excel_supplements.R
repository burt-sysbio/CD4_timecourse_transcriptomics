require(xlsx)
require(dplyr)
require(tidyr)

baseline = "degs_Th1_Th2_cor0.3"
dir3 = paste0("data/linear_model/linear_model_processed_" , baseline , ".csv")

df <- read.csv(dir3)


df_thm <- df %>% filter(model == "Th1/2")


myfun <- function(df){
  # clean gene names
  df <- df %>% select(gene, category) %>% 
    separate(gene, into = c("gene", "Iso"), sep = "\\.") %>% 
  select(-Iso) %>%
  distinct()
  
  nrow <- df %>% group_by(category) %>% count() 
  nrow <- max(nrow$n)
  
  
  col1 <- df$gene[df$category == "Th1like"]
  col2 <- df$gene[df$category == "Th2like"]
  col3 <- df$gene[df$category == "Ind"]
  col4 <- df$gene[df$category == "Superpos"]
  
  length(col1) <- nrow
  length(col2) <- nrow
  length(col3) <- nrow
  length(col4) <- nrow
  
  out <- data.frame("Th1like" = col1, "Th2like" = col2, "Independent" = col3, "Superposition" = col4)
  return(out)
}

out <- myfun(df_thm)

write.xlsx(out, "figures/supp/table_S6_superposition.xlsx", sheetName = "Gene categories Th1_2 hybrid cells", row.names = FALSE, showNA = FALSE)