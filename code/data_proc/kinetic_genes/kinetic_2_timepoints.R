# run masigpro, create design matrix on the way
require(dplyr)
require(pheatmap)
require(readr)
source("code/data_proc/utils_data_proc.R")

df <- read_csv("data/peine_data_log2fc_annot.csv")

df_th0 <- df[,1:10]
df_th1 <- df[,11:20]
df_th2 <- df[,21:30]
df_thmix <- df[,31:40]

thres <- 0.5

kin_th0 <- get_kin_genes(df_th0, thres = thres)
kin_th1 <- get_kin_genes(df_th1, thres = thres)
kin_th2 <- get_kin_genes(df_th2, thres = thres)
kin_thmix <- get_kin_genes(df_thmix, thres = thres)

df_kinetic <- tibble(gene = df$gene, 
                     kinetic_th0 = kin_th0$kinetic,
                     kinetic_th1 = kin_th1$kinetic, 
                     kinetic_th2 = kin_th2$kinetic,
                     kinetic_thmix = kin_thmix$kinetic)

sdir <- "data/kinetic_genes/"
sname <- paste0(sdir, "kinetic_genes_2_timepoints_logfc", thres, ".csv")

write.csv(df_kinetic, sname, row.names = F)
