require(dplyr)
require(tidyr)


# all genes that switch
df_switch <- read.csv("data/classswitch_quantDEGS/class_switch_summary_3clust.csv")

# all kinetic
df_kin <- read.csv("data/clustered_genes/kinetic_genes_20210714_3clust.csv")

df_kin_wide <- df_kin %>%  pivot_wider(names_from = celltype, values_from = annot)


df_cor <- read.csv("data/correlation/correlation_df.csv")

# p values anova
pvals <- readRDS("data/masigpro_output/pvec_fullrun_kinbackground_a0.01")

df_qual <- read.csv("data/classswitch_quantDEGS/quantDEGS_summary_cor0.3.csv")
df_quant <- read.csv("data/classswitch_quantDEGS/quantDEGS_summary_cor0.csv")

pval_df <- as.data.frame(pvals$p.vector)
pval_df$gene <- rownames(pval_df)


df <- df_cor[df_cor$gene %in% df_kin$gene, ]
df <- left_join(df, pval_df)

myfun <- function(df, contrast){
  # for a given contrast (DEG comparison) combine all metric into large dataframe

  # select relevant columns
  df1 <- df[c(paste0("cor_", contrast), "gene", "p.value")]
  colnames(df1) <- c("correlation", "gene", "pvalue_anova")
  
  # add cluster information
  if(contrast == "th1_th2"){
    mycols <- c("gene", "Th1", "Th2")
    
  } else if (contrast == "th1_th0"){
    mycols <- c("gene", "Th1", "Th0")
    
  } else if (contrast == "th2_th0"){
    mycols <- c("gene", "Th2", "Th0")
    
  } else if (contrast == "th1_thmix"){
    mycols <- c("gene", "Th1", "Th1/2")
    
  } else if (contrast == "th2_thmix"){
    mycols <- c("gene", "Th2", "Th1/2")
  }
  
  df_kin_wide <- df_kin_wide[mycols]
  colnames(df_kin_wide) <- c("gene", "clust_A", "clust_B")
  df1 <- left_join(df1, df_kin_wide)
  
  # check if class switch occurs based on cluster annotation
  df1$switch <- df1$clust_A != df1$clust_B

  # get qualitattive and quantitative genes for given comparison
  myquant <- df_quant$gene[df_quant$deg == contrast]
  myqual <- df_qual$gene[df_qual$deg == contrast]
  
  # if gene is qualitative (or quantitative) assign False (needed because of grouping, will be later reversed)
  df1$qual_deg <- T
  df1$quant_deg <- T
  df1$quant_deg[df1$gene %in% myquant] <- F
  df1$qual_deg[df1$gene %in% myqual] <- F
  
  # order columns
  df1 <- df1 %>% select(gene, correlation, pvalue_anova, quant_deg, qual_deg, everything())
  
  # clean isoforms and sort by correlation
  df1 <- df1 %>% separate(gene, into = c("gene", "iso"), sep = "\\.") %>% select(-iso) %>% distinct(gene, .keep_all = T)
  df1 <- df1 %>% group_by(qual_deg) %>% arrange(correlation, .by_group = TRUE)
  
  # reverse boolean annotation 
  df1$qual_deg<- !df1$qual_deg
  df1$quant_deg<- !df1$quant_deg
  df1 <- as.data.frame(df1)
  
  sname <- paste0("figures/supp/table_S1_", contrast, ".xlsx")
  write.xlsx(df1, sname, row.names = F)
  
  return(df1)
}

arr <- c("th1_th2", "th1_th0", "th2_th0", "th1_thmix", "th2_thmix")
for(contrast in arr){
  myfun(df, contrast)
}
