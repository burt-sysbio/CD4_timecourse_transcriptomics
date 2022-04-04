require(readr)
require(maSigPro)

# take established kinetic genes and run full masigpro but only on these genes


# get and save all the kinetic genes as a table
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

# make design matrix for single cell analysis (no comparisons, just check for each 
# cell individually
fdata <- readRDS("data/peine_featuredata.RDS")


df_kin <- read.csv("data/kinetic_genes/kinetic_genes_combined.csv")
# if a gene is kinetic in any cell type, keep it
kinetic_genes <- df_kin$gene[apply(df_kin[,2:5], 1, any)]


df <- df[rownames(df) %in% kinetic_genes,]
  
# set cell to "none" as default
cells <- c("Th2", "Th1", "Th0", "ThMix")
designInfo <- make_des_mat(fdata, cell = cells)

designMatrix <- make.design.matrix(designInfo, degree=n_degree)

# run masigpro
myFit <- p.vector(df, designMatrix, Q = Q, MT.adjust = "BH")
myTstep <- T.fit(myFit, step.method = "backward", alfa = alpha)

sname1 <- paste0("tfit_fullrun_kinbackgroundTH2_a", alpha)
sname2 <- paste0("pvec_fullrun_kinbackgroundTH2_a", alpha)

sdir <- "data/masigpro_output/"

saveRDS(myTstep, paste0(sdir, sname1))
saveRDS(myFit, paste0(sdir, sname2))

