require(readr)
require(dplyr)

# function that separates isoforms denoted by .1, .2 etc
clear_isoforms <- function(x){
  x <- data.frame(gene = x)
  x <- x %>% separate(gene, into = c("gene", "isoform"), sep = "\\.")
  x <- unique(x$gene)
  return(x)
}

# prepare class switches and class switches for DEGs for pathway analysis
mydate <- "20210825"
mydir <- "data/classswitch_quantDEGS/"
df_switch <- read.csv(paste0(mydir, "class_switch_summary_", mydate, "_3clust.csv"))


df_deg <- read.csv(paste0(mydir, "quantDEGS_summary_", mydate, "_3clust.csv"))

# get switches for DEGs, not that DEGs can also be in the same cluster!

use_deg <- F
if(use_deg == T){
  df_switch <- df_switch[df_switch$gene %in% df_deg$gene,]
}


# switches 1 -> 3
df_switch13 <- df_switch$gene[df_switch$switch1_3]

# switches 1 -> 2
df_switch12 <- df_switch$gene[df_switch$switch1_2]

# switches 2 -> 3
df_switch23 <- df_switch$gene[df_switch$switch2_3]

# switches th1 th2
df_switch_th1_th2 <- df_switch$gene[df_switch$s_th1_th2]

# switches th1 th0
df_switch_th1_th0 <- df_switch$gene[df_switch$s_th1_th0]

# switches th1 thmix
df_switch_th1_thmix <- df_switch$gene[df_switch$s_th1_thm]

# switches th2 thmix
df_switch_th2_thmix <- df_switch$gene[df_switch$s_th2_thm]

# switches th2 th0
df_switch_th2_th0 <- df_switch$gene[df_switch$s_th2_th0]



mylist <- list(df_switch13, df_switch12,
            df_switch23, df_switch_th1_th2,
            df_switch_th1_th0, df_switch_th1_thmix,
            df_switch_th2_thmix, df_switch_th2_th0)

fnames <- c("classswitch_clust13", "classswitch_clust12",
            "classswitch_clust23", "classswitch_th1_th2",
            "classswitch_th1_th0", "classswitch_th1_thmix",
            "classswitch_th2_thmix", "classswitch_th2_th0")

if(use_deg == T){
  fnames <- paste0(fnames, "deg.txt")
}else{
  fnames <- paste0(fnames, ".txt")
}

fnames <- paste0("genesets/", mydate, "/", mydate, "_classswitches/", fnames)

# clear isoform names and get unique genes
mylist2 <- lapply(mylist, clear_isoforms)

# store output
mapply(write.table, mylist2, fnames, MoreArgs = list(row.names = F, col.names = F, quote = F))
