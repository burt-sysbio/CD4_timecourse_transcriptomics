require(readxl)
require(dplyr)
require(stringr)
require(ggplot2)

# function that checks if pathway info is available and only then returns read in file, otherwise dummy df
myfun <- function(name){
  
  dirs <- str_split(name, pattern = "/")
  category <- dirs[[1]][5]

  out <- tibble()
  if(file.exists(name)){
    df <- read_xlsx(name)
    if("pvalue" %in% colnames(df)){
      out <- df
      out$DB <- category
    }
  }
  return(out)
}

# function to take comparison, read file and combine data bases
myfun2 <- function(comparison, mydir, categories){
  mytable <- "table.xlsx"
  fname <- paste0("degs_cleaned_20210825_", comparison, ".txt")
  name <- paste(mydir, fname, categories, mytable, sep = "/")
  
  # read all files that are found in databases (databases specifief in category)
  files <- lapply(name, myfun)
  out <- bind_rows(files)
  out$name <- comparison
  return(out)
}


proc_dotplot <- function(category, df, filter, alpha){
  
  if(category== "noGO"){
    df <- df[df$DB != "GO_BP", ]
  } else if(category != "all"){
    df <- df[df$DB == category, ]
  } 

  # alpha: keep only FDR below 0.05; filter: keep only IDs who turn up more than n times
  df_signif <- df[df$qvalue < alpha,]
  # count how often pathway is significant
  counts <- df_signif %>% group_by(ID) %>% summarize(n = n())
  myselec <- counts$ID[counts$n>=filter]
  
  # reduce df
  df_selec <- df_signif %>% filter(ID %in% myselec)
  
  # order df selec
  df_selec <- df_selec %>% mutate(new_ID = paste(DB, Description, sep = "_"))
  df_selec$new_ID <- as.factor(df_selec$new_ID)

  
  sname <- paste0("dotplot_", category, "_alpha", alpha, "_filter", filter, ".csv")
  write.csv(df_selec, paste0("data/pathway_analysis/sebastian_processed/", sname),
            row.names = F)
  return(df_selec)
}


mydir <- "data/pathway_analysis/pathways_sebastian"
categories <- c("GO_BP", "KEGG", "Msigdb_C2_Wikipathways", "Msigdb_H", "Msigdb_TFT", "Reactome")
comparisons <- c("Th1_Th0", "Th1_Th2", "Th1_ThMix", "Th2_Th0", "Th2_ThMix")

# for all databases and comparisons, read data
df_list <- lapply(comparisons, myfun2, mydir, categories)
df <- bind_rows(df_list)
df$fdrlog10 <- -log10(df$qvalue)


categories_dotplot <- c("KEGG", "Msigdb_C2_Wikipathways", "Msigdb_H", "Msigdb_TFT", "Reactome", "all", "noGO")
out <- lapply(categories_dotplot, proc_dotplot, df, 2, 0.01)
out <- lapply(categories_dotplot, proc_dotplot, df, 2, 0.05)
out <- lapply(categories_dotplot, proc_dotplot, df, 1, 0.01)
out <- lapply(categories_dotplot, proc_dotplot, df, 1, 0.05)


