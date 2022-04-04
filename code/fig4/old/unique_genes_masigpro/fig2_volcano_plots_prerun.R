require(readr)
require(ggplot2)
require(tidyr)
require(dplyr)
require(maSigPro)
source("code/fig2/utils.R")

fdata <- readRDS("data/peine_featuredata.RDS")
df <- read_csv("data/peine_data_median_threshold_annot.csv")

gene <- df$gene
df <- df %>% select(-gene)

n_degree = 4
Q = 0.05
alpha = 0.05

cell1 <- "Th1_Th2_Th0"
cell2 <- "Th1_Th2_ThMix"

fit1 <- get_degs(df, fdata, cell1)
fit2 <- get_degs(df, fdata, cell2)

saveRDS(fit1, "r_objects/fig2/mspro_pvec_th0.RDS")
saveRDS(fit2, "r_objects/fig2/mspro_pvec_thm.RDS")
