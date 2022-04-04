
#prepare data frame for class switches
require(dplyr)
require(tidyr)


mydir <- paste0("genesets_output/quantDEGS/")

threshold_array <- as.character(seq(0,0.5,0.01))

for(cor_thres in threshold_array){
  th1_th0 <- read.csv(paste0(mydir, "degs_Th1_Th0_cor", cor_thres, ".csv"))
  th1_th2 <- read.csv(paste0(mydir, "degs_Th1_Th2_cor", cor_thres, ".csv"))
  th1_thmix <- read.csv(paste0(mydir, "degs_Th1_ThMix_cor", cor_thres, ".csv"))
  th2_th0 <- read.csv(paste0(mydir, "degs_Th2_Th0_cor", cor_thres, ".csv"))
  th2_thmix <- read.csv(paste0(mydir, "degs_Th2_ThMix_cor", cor_thres, ".csv"))
  
  th1_th0$deg <- "th1_th0"
  th1_th2$deg <- "th1_th2"
  th1_thmix$deg <- "th1_thmix"
  th2_th0$deg <- "th2_th0"
  th2_thmix$deg <- "th2_thmix"
  
  df <- bind_rows(th1_th0, th1_th2, th1_thmix, th2_th0, th2_thmix)
  
  writedir <- paste0("data/classswitch_quantDEGS/quantDEGS_summary_cor", cor_thres, ".csv")
  write.csv(df, writedir, row.names = F)
}


# load information about degs