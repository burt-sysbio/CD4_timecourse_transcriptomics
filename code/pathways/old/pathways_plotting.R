require(readr)
require(dplyr)
require(ggplot2)
require(cowplot)
require(stringr)

# read all files from pathway analysis
mypath <- "data/Ora_output/ORA_burt/"
filenames <- list.files(path = mypath, pattern = ".txt", full.names = T)
filenames2 <- list.files(path = mypath, pattern = ".txt", full.names = F)

filenames_save <- str_sub(filenames2, 1, -5)
files <- lapply(filenames, read.table, header = T)


# for each file make barplot of top fdr log10 transformed pathways
for(i in seq_along(files)){
  
  df <- files[[i]]
  
  # from string extract the database
  sname <- filenames_save[[i]]
  category <- str_split(sname, "_")
  category <- category[[1]][2]
  sname <- paste0("figures/pathway_results/", category, "/", sname, ".png")
  print(sname)

  df$fdrlog10 <- -log10(df$qvalue)
  
  df_sort <- df %>% arrange(desc(fdrlog10)) %>% slice_max(fdrlog10, n = 30)
  
  # process category to remove hallmark etc from names
  df_sort <- df_sort %>% mutate(Description = str_sub(Description, (nchar(category)+2), -1))
  
  
  g <- ggplot(data = df_sort, aes(fdrlog10, reorder(Description, fdrlog10)))
  g + geom_col() +
    xlab("-log10(padj)") +
    ylab("") +
    theme_cowplot()
  
  ggsave(sname, width = 15, height = 10, bg = "white")
}

