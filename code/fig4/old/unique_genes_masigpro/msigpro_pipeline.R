# run masigpro, create design matrix on the way
library(maSigPro)
source("code/fig2/make_design_matrix_new.R")
# load replicate 1 data
load("data/4_ready_to_use/allTh_Rep1.Rdata")
load("data/4_ready_to_use/allTh_names_Rep1.Rdata")

# remove lowly/ no-express genes
threshold <- median(allTh_1)
rowmax <- apply(allTh_1, 1, max)
allTh_1 <- allTh_1[which(rowmax >= threshold),]
dim(allTh_1) #18284 

# geneNames for probes IDs
geneIds <- readRDS("annotation/burt_annot.rds")

### MaSigPro - per cell subset with reduced FDR
cell = "Th1"
stopifnot(cell %in% c("Th1", "Th2", "ThMix", "Th0", "Th1_Th2_ThMix", "Th1_Th2_Th0", "Th1_Th2"))
# masigpro params
n_degree = 4
Q = 0.05
alpha = 0.05

# time points to analyze
keep_3 <- c(0, 48, 120)
keep_4 <- c(0, 24, 73, 120)
keep_5 <- c(0, 12, 35, 73, 120)
keep_6 <- c(0, 6, 24, 48, 96, 120)
keep_7 <- c(0, 3, 12, 35, 73, 96, 120)
keep_8 <- c(0, 3, 6, 12, 35, 48, 96, 120)
keep_9 <- c(0, 3, 6, 12, 24, 35, 73, 96, 120)
keep_all <- NULL
timepoints <- c(0,3,6,12,24,35,48,73,96,120)
tp <- rep(timepoints, 4)

# assign which timepoints to keep here
runs <- list(keep_6)

for(keep_timepoints in runs){
  
  df_data <- allTh_1

  # kick out not used columns if keep timepoints is not Null
  if(length(keep_timepoints)!=0){
    stopifnot(n_degree <= length(keep_timepoints)-1)
    keep_cols <- tp %in% keep_timepoints
    df_data <- df_data[, keep_cols]
  }

  # only keep columns in data that match time points to be kept
  # get design info for cell type
  # check if Th0 are needed and get appropriate design
  if (cell %in% c("Th1", "Th2", "ThMix", "Th1_Th2_ThMix", "Th1_Th2")){
    use_th0 = F
    toRemove <- grep("Th0",colnames(df_data))
    
    # remove th0 columns
    df_data <- df_data[,-toRemove]
    designInfo <- make_des_mat(allTh_names_1, cell, use_th0 = use_th0, keep_timepoints = keep_timepoints)  
  } else {
    use_th0 = T
    designInfo <- make_des_mat(allTh_names_1, cell, use_th0 = use_th0, keep_timepoints = keep_timepoints)  
  }
  
  # reduce df and design to appropriate cell type
  if ((cell!="Th1_Th2_ThMix") & (cell!="Th1_Th2_Th0") & (cell!= "Th1_Th2")){
    # subset expr data
    keep <- grep(cell, colnames(df_data))
    df_data <- df_data[,keep]
  }
  
  # run masigpro
  designMatrix <- make.design.matrix(designInfo, degree=n_degree)

  myFit <- p.vector(df_data, designMatrix, Q = Q, MT.adjust = "BH")
  
  # only continue analysis if DEGs are found at all
  print("n signif genes found")
  print(myFit$i)
  if(myFit$i>0){
    myTstep <- T.fit(myFit, step.method = "backward", alfa = alpha)
    
    if (length(keep_timepoints)>0){
      saveRDS(myTstep, paste0("r_objects/", "msgpro_",cell,"_keep_", length(keep_timepoints), "_Q",Q,"a",alpha,".RDS"))
    } else {
      saveRDS(myTstep, paste0("r_objects/", "msgpro_",cell,"_Q",Q,"a",alpha,".RDS"))
    }    
  } else {
    print("timepoints")
    print(keep_timepoints)
  }
}
