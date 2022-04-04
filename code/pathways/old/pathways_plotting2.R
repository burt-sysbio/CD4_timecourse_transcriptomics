# take multiple pathways and plot them together in dotplot

require(readr)
require(tidyr)
require(ggplot2)
require(stringr)
require(dplyr)

mypath <- "data/ORA_output/ORA_burt/"
filenames <- list.files(mypath, pattern = ".txt", full.names = T)
filenames2 <- list.files(mypath, pattern = ".txt", full.names = F)


# get only the deg files (but not for Th0)
filenames <- filenames[grepl("cor0.5", filenames)]
filenames2 <- filenames2[grepl("cor0.5", filenames2)]

# read files
files <- lapply(filenames, read.table, header = T)

myfun1 <- function(df, fname){
  # add the filename as a column, transform filename first
  name <- str_split(fname, "_")
  name <- paste0(name[[1]][5], "vs", name[[1]][6])
  df$name <- name
  return(df)
}
out <- mapply(myfun1, files, filenames2, SIMPLIFY = F)
out <- bind_rows(out)

# add database as column
out <- out %>% separate(ID, into = c("DB"), remove = F)

# kick out GO?
out <- out %>% filter(DB != "GOBP")

# remove nonsignif categories
myfun2 <- function(df, alpha = 0.05, filter = 2){
  df <- df %>% filter(qvalue <= alpha)
  counts <- df %>% group_by(ID) %>% count()
  myselec <- counts$ID[counts$n>=filter]
  df <- df %>% filter(ID %in% myselec)  
  return(df)
}

out <- myfun2(out)

# remove the category stuff from plot
#out <- out %>% mutate(ID = str_sub(ID, (nchar(category)+2),-1))
#S1<- ggplot(out, aes(x= name, y=ID, size=Count, color=fdrlog10, group=name)) + geom_point(alpha = 0.8) + 
#  theme_bw(base_size = 15)

#S1 = S1+scale_color_gradient(low = "lightgrey",  high = "black", space = "Lab")
#S1+scale_size(range = c(2, 8)) + xlab("") + ylab("") +
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#sname <- paste0("figures/pathway_results/pathway_dotplot_", category, ".png")
#ggsave(sname, bg = "white", width = 12, height = 8)