
# combine output from masigpro into one data set indicating whether gene is kinetic or not kinetic
require(readr)
require(maSigPro)
source("code/data_proc/utils_data_proc.R")
require(tidyr)


df <- read_csv("data/peine_data_log2fc_annot.csv")

# load Th1 masigpro fit
Q = "0.05"
a = "0.05"


kin_th1 <- read.csv("data/masigpro_output/mspro_genes_contrast_Th1.csv")
kin_th2 <- read.csv("data/masigpro_output/mspro_genes_contrast_Th2.csv")
kin_th0 <- read.csv("data/masigpro_output/mspro_genes_contrast_Th0.csv")
kin_thmix <- read.csv("data/masigpro_output/mspro_genes_contrast_ThMix.csv")


df_kin <- tibble(gene = df$gene, 
                     kinetic_th1 = F, 
                     kinetic_th2 = F,
                     kinetic_th0 = F,
                     kinetic_thmix = F)


df_kin$kinetic_th1[df$gene %in% kin_th1$x] <- T
df_kin$kinetic_th2[df$gene %in% kin_th2$x] <- T
df_kin$kinetic_th0[df$gene %in% kin_th0$x] <- T
df_kin$kinetic_thmix[df$gene %in% kin_thmix$x] <- T

sdir <- "data/kinetic_genes/"
sname <- paste0(sdir, "kinetic_genes_masigpro_summary.csv")

write.csv(df_kin, sname, row.names = F)
