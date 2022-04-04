# create design info matrix for masigpro from zuzanna data, see masigpro user guide how design matrix should look like

make_des_mat <- function(df){

  designInfo <- df[,1:3]
  colnames <- c("ThMix","Th1","Th2", "Th0")

  designInfo <- as.data.frame(designInfo)
  # kick out timepoints that sould be ommitted
  if(length(keep_timepoints) != 0){
    designInfo <- designInfo[designInfo[,"timePoint"] %in% keep_timepoints, ]
  }  

  # add new column per cell type filled with zeros, then add 1 if celltype row matches column name
  for(i in colnames){
    # add new column with title and fill with zeros
    designInfo[i] <- 0
    # add 1 if coltitle pops up in first column
    designInfo[which(designInfo[,1] == i), i] <- 1
  }

  # kick out celltype column now
  designInfo <- designInfo[,-1]
  # reorder if th0 is comparison
  

  colorder <- c("timePoint", "replicate", "Th0", "Th1", "Th2", "ThMix")
  designInfo <- designInfo[,colorder]

  # fill replicate column with sequence 1,2,.. but omit places with time point zero
  # this need to be done at last step, otherwise numbering will come out wrong
  sequence <- seq(1:length(designInfo[,1]))
  designInfo[,"replicate"] = 1
  # indexing vector for sequence (2 is correct)
  j <- 2
  for (i in seq_along(sequence)){
    if (designInfo[i,"timePoint"] != 0){
      designInfo[i,"replicate"] <- j
      j <- j+1
    } else {
      # for a commong starting point at t=0 set all replicate information to 1
      designInfo[i, 2:ncol(designInfo)] <- 1
    }
  }  
 
  designInfo <- as.matrix(designInfo)
  designInfo <- apply(designInfo,1:2,as.numeric)
  
  return(designInfo)
}