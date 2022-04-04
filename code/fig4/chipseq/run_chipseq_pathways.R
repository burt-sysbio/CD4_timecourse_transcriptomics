
run_pathways <- function(baseline){
  sdir <- "data/linear_model/"
  
  mygenes <- read.csv(paste0(sdir, "linear_model_processed_", baseline, ".csv"))
  mygenes <- mygenes %>% filter(model == "Th1/2")
  
 
  background <- mygenes %>% separate(gene, sep = "\\.", into = c("gene", "ISO", remove = F))
  universe <- as.character(background$gene)
  universe <- unique(universe)
  
  
  mygenes_hclust <- read.csv("data/correlation/hclust_output/hclust_corr_thm_no_superpos.csv")
  mygenes_hclust2 <- mygenes_hclust %>% filter(category == "Superpos") %>% mutate(category = "Superpos_Ind")
  mygenes_hclust3 <- mygenes_hclust %>% filter(category == "Ind") %>% mutate(category = "Superpos_Ind")
  
  mygenes_hclust <- bind_rows(mygenes_hclust, mygenes_hclust2, mygenes_hclust3)
  
  
  background_hclust <- mygenes_hclust %>% separate(gene, sep = "\\.", into = c("gene", "ISO", remove = F)) 
  background_hclust <- unique(background_hclust$gene)
  
  
  # clean isoforms for each model
  myfun <- function(df, model = "Ind"){
    df <- df[df$category %in% model,]
    df <- df %>% separate(gene, sep = "\\.", into = c("gene", "ISO"), remove = FALSE)
    out <- unique(df$gene)
    return(out)
  }
  
  #TFT_genes <- msigdbr(species = "Mus musculus", category = "C3", subcategory = "TFT:TFT_Legacy")
  #TFT_genes <- as.data.frame(TFT_genes) %>% select(gs_name, gene_symbol)
  
  chipseq_genes <- read.csv("genesets_literature/references/chipseq_genes_all.csv")
  chipseq_genes <- as.data.frame(chipseq_genes)
  
  mynames <- list("Ind", "Superpos", "Th1like", "Th2like")
  
  # run pathway analysis for each model
  myfun2 <- function(name){
    s1 <- myfun(mygenes, model = name)
    s2 <- myfun(mygenes_hclust, model = name)
    
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
    
    res_s2 <- enricher(s2, 
                       TERM2GENE = chipseq_genes, 
                       universe = background_hclust, 
                       pAdjustMethod = "fdr", 
                       maxGSSize = 5000, 
                       pvalueCutoff = 1, 
                       qvalueCutoff = 1)
    res2 <- res_s2@result
    res2$model <- name
    res2$approach <- "hclust"
    
    res <- bind_rows(res1,res2)
    return(res)
  }
  
  df <- lapply(mynames,myfun2)
  df_out <- bind_rows(df)
  write.csv(df_out, paste0("code/fig4/chipseq/chipseq_pathway_results_", baseline, ".csv"), row.names = F)
  return(df_out)
}
