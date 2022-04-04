require(readr)
require(dplyr)

make_correlation_df <- function(df){
  df_th0 <- df[,1:10]
  df_th1 <- df[,11:20]
  df_th2 <- df[,21:30]
  df_thmix <- df[,31:40]
  
  cor_th1_th2 <- diag(cor(t(df_th1), t(df_th2)))
  cor_th1_th0 <- diag(cor(t(df_th1), t(df_th0)))
  cor_th1_thmix <- diag(cor(t(df_th1), t(df_thmix)))
  cor_th2_th0 <- diag(cor(t(df_th2), t(df_th0)))
  cor_th2_thmix <- diag(cor(t(df_th2), t(df_thmix)))
  cor_th0_thmix <- diag(cor(t(df_th0), t(df_thmix)))
  
  df_cor <- data.frame(cor_th1_th2,cor_th1_th0,cor_th1_thmix,cor_th2_th0,cor_th2_thmix,cor_th0_thmix)
  df_cor$gene <- df$gene
  
  write.csv(df_cor, "data/correlation_df.csv", row.names = F)
  
  
}