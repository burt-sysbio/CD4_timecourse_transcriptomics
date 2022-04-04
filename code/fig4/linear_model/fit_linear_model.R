
# process df, then apply linear model and kmeans clustering
pipeline <- function(df, use_thm){
  if (use_thm){
    remove <- "Th0"
    use <- "Th1/2"
  } else {
    remove <- "Th1/2"
    use <- "Th0"
  }
  df <- df %>% filter(celltype != remove)
  
  # kick out th0 or thm and then make wider
  df2 <- df %>% pivot_wider(names_from = celltype, values_from = value)
  
  # new names for better indexing
  colnames(df2)[colnames(df2) == use] <- "Y"
  
  # fit linear models to each gene
  df_list <- df2 %>% group_by(gene) %>% group_split()
  model_output <- lapply(df_list, myfun)
  
  # collapse output from multiple genes into one dataframe
  model_fits <- bind_rows(lapply(model_output, `[[`, 1))
  model_fits$model <- use
  
  anova_results <- bind_rows(lapply(model_output, `[[`, 2))
  anova_results$model <- use
  # return ftest results for model comparisons and individual model fits
  return(list("anova_results" = anova_results, "fit_results" = model_fits))
}

# take data frame and apply linear model on it
myfun <- function(df){
  fit1 <- lm(Y ~ Th1, data = df)
  fit2 <- lm(Y ~ Th2, data = df)
  fit3 <- lm(Y ~ Th1+Th2, data = df)
  
  mygene <- unique(df$gene)
  sfit1 <- proc_fit(fit1, "Th1model", mygene)
  sfit2 <- proc_fit(fit2, "Th2model", mygene)
  sfit3 <- proc_fit(fit3, "mixedmodel", mygene)
  
  # run f test against mixed model
  aov13 <- anova(fit1, fit3, test = "F")
  aov23 <- anova(fit2, fit3, test = "F")
  aov_results <- data.frame(ftest = c("Th1vsMixed", "Th2vsMixed"),
                            pvalue = c(aov13$`Pr(>F)`[2], aov23$`Pr(>F)`[2]),
                            gene = c(mygene, mygene))
  fit_summary <- bind_rows(sfit1, sfit2, sfit3)
  out <- list("lm_fits" = fit_summary, "Ftest" = aov_results)
  return(out)
}

# take linear model output and grab coefficients and gene 
proc_fit <- function(fit_result, fit_name, gene){
  out <- summary(fit_result)
  coeffs <- as.data.frame(out$coefficients)
  colnames(coeffs) <- c("estimate", "std_error", "t_value", "pvalue")
  coeffs$gene <- gene
  coeffs$covariate <- rownames(coeffs)
  coeffs$fit_name <- fit_name
  return(coeffs)
}


fit_model <- function(baseline){

  if(baseline == "kinetic_th12"){
    df_kin <- read.csv("data/kinetic_genes/kinetic_genes_combined.csv")
    df_kin <- df_kin[df_kin$kinetic_thmix,]
    baseline_df <- df_kin

  } else {
    baseline_df <- read.csv(paste0("genesets_output/quantDEGS/", baseline, ".csv"))
  }
  
  # take dataframe and filter by baseline genes
  df <- read.csv("data/peine_logfc_tidy.csv")
  df <- df %>% filter(gene %in% baseline_df$gene)
  
  # run linear models and ftests
  df1 <- pipeline(df, use_thm = T)
  df2 <- pipeline(df, use_thm = F)
  
  fitres_1 <- df1$fit_results
  fitres_2 <- df2$fit_results
  
  aov1 <- df1$anova_results
  aov2 <- df2$anova_results
  
  # store output
  fit_res <- bind_rows(fitres_1, fitres_2)
  aov <- bind_rows(aov1, aov2)
  
  sdir <- "data/linear_model/"
  write.csv(fit_res, paste0(sdir, "linear_model_fitresults_", baseline, ".csv"), row.names = F)
  write.csv(aov, paste0(sdir, "linear_model_ftest_", baseline, ".csv"), row.names = F)
  print(aov)
  
}
