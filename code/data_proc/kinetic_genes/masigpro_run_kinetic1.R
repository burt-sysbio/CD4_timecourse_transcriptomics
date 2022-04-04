require(readr)
require(maSigPro)

# get and save all the kinetic genes as a table

source("code/data_proc/utils_data_proc.R")

# load data
df <- read_csv("data/data_log2norm/data_norm_d0.csv")
df <- as.data.frame(df)
rownames(df) <- df$gene
df$gene <- NULL

# set params
n_degree = 4
Q = 0.05
alpha = 0.05

# make design matrix for single cell analysis (no comparisons, just check for each 
# cell individually
fdata <- readRDS("data/misc/peine_featuredata.RDS")

cells <- c("Th1", "Th2", "ThMix", "Th0")
for (cell in cells){

  # set cell to "none" as default
  designInfo <- make_des_mat(fdata, cell = cell)

  designMatrix <- make.design.matrix(designInfo, degree=n_degree)
  
  # run masigpro
  myFit <- p.vector(df, designMatrix, Q = Q, MT.adjust = "BH")
  myTstep <- T.fit(myFit, step.method = "backward", alfa = alpha)
  
  # get degs
  degs <- get.siggenes(myTstep, significant.intercept = "all", vars = "all", rsq = 0)
  degs_summary <- degs$summary
  # get kinetic genes
  sname <- paste0("data/masigpro_output/mspro_genes_contrast_", cell, ".csv")
  write.csv(degs_summary, sname, row.names = F) 
}