require(readr)
require(maSigPro)

source("code/fig2/unique_genes_masigpro/make_design_matrix_new.R")

# load data
df <- read_csv("data/peine_data_median_threshold_annot.csv")
df <- as.data.frame(df)
rownames(df) <- df$gene
df$gene <- NULL

# set params
n_degree = 4
Q = 0.05
alpha = 0.05

# make design matrix
fdata <- readRDS("data/peine_featuredata.RDS")
designInfo <- make_des_mat(fdata)
designMatrix <- make.design.matrix(designInfo, degree=n_degree)

# run masigpro
myFit <- p.vector(df, designMatrix, Q = Q, MT.adjust = "BH")
myTstep <- T.fit(myFit, step.method = "backward", alfa = alpha)

saveRDS(myFit, "r_objects/fig2/mspro_pvec.RDS")
saveRDS(myTstep, "r_objects/fig2/mspro_tstep.RDS")
