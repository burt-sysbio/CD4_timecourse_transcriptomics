require(dplyr)
require(tidyr)
require(ggplot2)
require(stringr)

mydir <- "data/pathway_analysis/sebastian_processed/"

files <- list.files(mydir, ".csv", full.names = T)
files2 <- list.files(mydir, ".csv", full.names = F)
files2 <- str_sub(files2, 1, -5)
snames <- paste0("figures/pathway_results/sebastian/", files2, ".png")


files <- lapply(files, read.csv)


plot_dotplot <- function(df, fname){
  #df$DB <- as.factor(df$DB)
  S1 <- ggplot(df, aes(x= name, y=new_ID, size=Count, color=fdrlog10)) + geom_point(alpha = 0.8) + 
    theme_bw(base_size = 15)
  
  S1 = S1+scale_color_gradient(low = "lightgrey",  high = "black", space = "Lab")
  S2 <- S1+scale_size(range = c(2, 8)) + xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #+
  #scale_y_discrete(labels = reorder(df$Description,df$new_ID))#

  if(grepl("dotplot_all", fname)){
    width <- 18
    height <- 25
    #S3 <- S2 + facet_wrap(~DB, nrow = 5, scales = "free_y")
  } else {
    width <- 10
    height <- 7
  }
  ggsave(fname, width = width, height = height, bg = "white")
  
}

replace_database <- function(df, db_name, new_name){
  df$new_ID <- str_replace(df$new_ID, db_name, new_name)
  df$DB <- str_replace(df$DB, db_name, new_name)
  return(df)
}

proc_df <- function(df){
  df <- replace_database(df, "Msigdb_TFT", "TFT")
  df <- replace_database(df, "Msigdb_C2_Wikipathways", "WIKIPATH")
  df <- replace_database(df, "Reactome", "REACTOME")
  
  # manually change some names
  df$new_ID[grepl("TUBERCULOSIS", df$new_ID)] <- "WIKIPATH_Immune Response to Tuberculosis"
  # manually change some names
  df$new_ID[grepl("SELECTIVE", df$new_ID)] <- "WIKIPATH_Expression of Chemokine Receptor"
  # manually change some names
  df$new_ID[grepl("HOSTPATHOGEN", df$new_ID)] <- "WIKIPATH_Hostpathogen Interaction Coronavirus"
  # manually change some names
  df$new_ID[grepl("Antiviral mechanism", df$new_ID)] <- "REACTOME_Antiviral IFN-stimulated genes"
  
  # remove reduntant category
  df <- df[!grepl("STTT", df$new_ID),]  
  return(df)
}


out <- lapply(files, proc_df)
mapply(plot_dotplot, out, snames)
