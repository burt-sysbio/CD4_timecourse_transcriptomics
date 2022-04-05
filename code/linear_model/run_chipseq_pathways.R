require(clusterProfiler)
require(msigdbr)
require(dplyr)
require(tidyr)
require(readr)
require(stringr)

run_pathways <- function(baseline){
  sdir <- "data/linear_model/"
  
  mygenes <- read.csv(paste0(sdir, "linear_model_processed_", baseline, ".csv"))
  mygenes <- mygenes %>% filter(model == "Th1/2")
  
 
  background <- mygenes %>% separate(gene, sep = "\\.", into = c("gene", "ISO", remove = F))
  universe <- as.character(background$gene)
  universe <- unique(universe)
  
  # clean isoforms for each model
  myfun <- function(df, model = "Ind"){
    df <- df[df$category %in% model,]
    df <- df %>% separate(gene, sep = "\\.", into = c("gene", "ISO"), remove = FALSE)
    out <- unique(df$gene)
    return(out)
  }
  
 
  chipseq_genes <- read.csv("genesets_input/chipseq_genes_all.csv")
  chipseq_genes <- as.data.frame(chipseq_genes)
  
  mynames <- list("Ind", "Superpos", "Th1like", "Th2like")
  
  # run pathway analysis for each model
  myfun2 <- function(name){
    s1 <- myfun(mygenes, model = name)

    res_s1 <- enricher(s1, 
                       TERM2GENE = chipseq_genes, 
                       universe = universe, 
                       pAdjustMethod = "fdr", 
                       maxGSSize = 5000, 
                       pvalueCutoff = 1, 
                       qvalueCutoff = 1)
    res1 <- res_s1@result
    res1$model <- name
    res1$approach <- "linear model"
    
    return(res1)
  }
  
  df <- lapply(mynames,myfun2)
  df_out <- bind_rows(df)
  write.csv(df_out, paste0("code/linear_model/chipseq_pathway_results/chipseq_pathway_results_", baseline, ".csv"), row.names = F)
  return(df_out)
}
