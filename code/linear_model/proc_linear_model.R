require(dplyr)
require(tidyr)
require(readr)
require(stringr)

# get the p values for Th1 and Th2 regression coefficients 
# then assign categories (superpos --> th1 and th2 coeff are both significant)
proc_output <- function(df, use_fdr = T){
  
  df <- df %>% filter(covariate != "(Intercept)")
  # first get value for the Th1 and Th2 estimates
  df_est <- df %>% select(estimate, gene, covariate) %>% pivot_wider(names_from = covariate, values_from = estimate)
  colnames(df_est) <- c("gene", "beta_th1", "beta_th2")
  
  df <- df %>% select(pvalue, gene, covariate)
  df <- df %>% pivot_wider(names_from = covariate, values_from = pvalue)
  colnames(df) <- c("gene", "p_th1", "p_th2")
  
  if(use_fdr){
    df$padj_th1 <- p.adjust(df$p_th1, method = "fdr")
    df$padj_th2 <- p.adjust(df$p_th2, method = "fdr")   
  } else {
    df$padj_th1 <- df$p_th1 
    df$padj_th2 <- df$p_th2 
  }

  df$th1_sig <- "ns"
  df$th1_sig[df$padj_th1<=0.05] <- "sig"
  df$th2_sig <- "ns"
  df$th2_sig[df$padj_th2<=0.05] <- "sig"
  
  df$category <- "Ind"
  df$category[(df$th1_sig) == "ns" & (df$th2_sig) == "sig"] <- "Th2like"
  df$category[(df$th2_sig) == "ns" & (df$th1_sig) == "sig"] <- "Th1like"
  df$category[(df$th2_sig) == "sig" & (df$th1_sig) == "sig"] <- "Superpos"
  
  df <- inner_join(df_est, df)
  
  return(df)
}


proc_model <- function(baseline){
  #"kinetic_th12"
  sdir <-"data/linear_model/"
  fit_res <- read.csv(paste0(sdir, "linear_model_fitresults_", baseline, ".csv"))
  aov_res <- read.csv(paste0(sdir, "linear_model_ftest_", baseline, ".csv"))
  
  fit_res1 <- fit_res %>% filter(fit_name == "mixedmodel", model == "Th1/2")
  fit_res2 <- fit_res %>% filter(fit_name == "mixedmodel", model == "Th0")
  
  fit_res1 <- proc_output(fit_res1)
  fit_res2 <- proc_output(fit_res2)
  
  fit_res1$model <- "Th1/2"
  fit_res2$model <- "Th0"
  out <- bind_rows(fit_res1, fit_res2)
  
  
  write.csv(out, paste0(sdir, "linear_model_processed_", baseline, ".csv"))
  return(out)
}
