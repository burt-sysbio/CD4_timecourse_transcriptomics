require(readr)
require(maSigPro)

# get pvalues for each comparison with modified pvcetor script to access all pvals


# set wd explicitly in case script is used in server call
setwd("~/Documents/projects/2020/kinetic_paper/")
source("code/data_proc/utils_data_proc.R")

# load data
df <- read_csv("data/peine_data_log2fc_annot.csv")
df <- as.data.frame(df)
rownames(df) <- df$gene
df$gene <- NULL

# set params
n_degree = 4
Q = 0.05
alpha = 0.05

# make design matrix
fdata <- readRDS("data/peine_featuredata.RDS")


cell <- c("Th1", "Th2")

# set cell to "none" as default
designInfo <- make_des_mat(fdata, cell = cell)
print(designInfo)
designMatrix <- make.design.matrix(designInfo, degree=n_degree)

# run masigpro
#myFit1 <- p.vector_modified(df, designMatrix, Q = Q, MT.adjust = "BH")

# check if output from Tstep is the same as the pvals < 0.05 (for rsq=0)
myFit2 <- p.vector(df, designMatrix, Q = Q, MT.adjust = "BH")
myTfit = T.fit(myFit2, step.method = "backward", alfa = alpha)

# run tstep to get rsq value

