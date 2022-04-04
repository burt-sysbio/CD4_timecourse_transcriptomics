
# run pathways but based on superposition by angle cluster criterion

sdir <- "data/linear_model/"
baseline <- "degs_Th1_Th2_cor0"
mygenes <- read.csv(paste0(sdir, "linear_model_processed_", baseline, ".csv"))
mygenes <- mygenes %>% filter(model == "Th1/2")


background <- mygenes %>% separate(gene, sep = "\\.", into = c("gene", "ISO", remove = F))
universe <- as.character(background$gene)
universe <- unique(universe)


# clean isoforms for each model
myfun <- function(df, model = "Ind"){
  df <- df[df$clust_new %in% model,]
  df <- df %>% separate(gene, sep = "\\.", into = c("gene", "ISO"), remove = FALSE)
  out <- unique(df$gene)
  return(out)
}


chipseq_genes <- read.csv("genesets_literature/references/chipseq_genes_all.csv")
chipseq_genes <- as.data.frame(chipseq_genes)

mynames <- list("Superpos", "Th1like", "Th2like")

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
