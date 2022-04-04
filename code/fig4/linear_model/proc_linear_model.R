require(tidyr)
require(dplyr)
require(ggplot2)

# get the p values for Th1 and Th2 regression coefficients 
# then assign categories (super pos is th1 and th2 coeff are both signif)
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
    df$padj_th1 <- df$p_th1 #p.adjust(df$p_th1, method = "fdr")
    df$padj_th2 <- df$p_th2 #p.adjust(df$p_th2, method = "fdr")
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
  
  #test <- testfun(fit_res[fit_res$model == "Th1/2",], aov_res[aov_res$model== "Th1/2",])
  
  fit_res1 <- proc_output(fit_res1)
  fit_res2 <- proc_output(fit_res2)
  
  fit_res1$model <- "Th1/2"
  fit_res2$model <- "Th0"
  out <- bind_rows(fit_res1, fit_res2)
  
  
  write.csv(out, paste0(sdir, "linear_model_processed_", baseline, ".csv"))
  return(out)
}

testfun <- function(df, aov){
  # use ftest from mixed model to identify superposition genes
  df1 <- df %>% filter(fit_name == "Th1model", covariate == "Th1")
  df2 <- df %>% filter(fit_name == "Th2model", covariate == "Th2")
  
  genes_ind1 <- df1$gene[df1$pvalue>0.05]
  genes_ind2 <- df2$gene[df2$pvalue>0.05]
  genes_ind <- genes_ind1[genes_ind1 %in% genes_ind2]

  # these are the genes that are at least significant for one fit
  genes_th1 <- df1$gene[df1$pvalue<0.05]
  genes_th2 <- df2$gene[df2$pvalue<0.05]
  
  # these genes are significant in both fits
  genes_superpos <- genes_th1[genes_th1 %in% genes_th2]
  
  # these genes are exclusively significant in only one fit
  genes_th1 <- genes_th1[!(genes_th1 %in% genes_superpos)]
  genes_th2 <- genes_th2[!(genes_th2 %in% genes_superpos)]
  
  # get the independent genes for which a mixed model did not improve the fit substantially
  superpos_aov <- aov %>% filter(pvalue > 0.05) %>% select(gene) %>% unique()
  
  # these are the genes that are not signif for either fit and also not for mixedmodel fit
  genes_ind <- genes_ind[!(genes_ind %in% superpos_aov$gene)]
  
  out = list(ind = genes_ind, th1like = genes_th1, th2like = genes_th2, superpos = union(superpos_aov$gene, genes_superpos))
  return(out)
}
