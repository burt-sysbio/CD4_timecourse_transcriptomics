load("data/replicate_2/allTh_names_Rep2.Rdata")
load("data/replicate_2/allTh_Rep2.Rdata")
require(dplyr)
require(tidyr)

# geneNames for probes IDs
geneIds <- readRDS("r_objects/annotation/burt_annot.rds")

# remove lowly/ no-express genes
threshold <- median(allTh_2)
rowmax <- apply(allTh_2, 1, max)
allTh_2 <- allTh_2[which(rowmax >= threshold),]
dim(allTh_2) #18284 

# annotate df with genes names
allTh_2 <- as.data.frame(allTh_2)
allTh_2$id <- rownames(allTh_2)
geneIds$id <- rownames(geneIds)
df <- inner_join(allTh_2, geneIds, by = "id")
rownames(df) <- df$Symbol

df <- df[,1:40]
# compute log2fc
# this works because for day 0 all columns should be identic
df_log2fc <- df / df[,1]
df_log2fc <- log2(df_log2fc)
df_log2fc$gene <- rownames(df)
df$gene <- rownames(df)

write.csv(df_log2fc, "data/replicate_2/peine_data_log2fc_annot_rep2.csv",row.names = F)
write.csv(df, "data/replicate_2/peine_data_median_threshold_annot_rep2.csv",row.names = F)


# make tidy
tidyfun <- function(df, sname){
  df <- pivot_longer(df, -gene)
  df <- separate(df, col = name, into = c("celltype", "time", NA), sep = "_")
  df$time <- as.numeric(df$time)
  write.csv(df, paste0("data/replicate_2/peine_tidy_rep2_", sname, ".csv"))
}

#tidyfun(df, "median_")
tidyfun(df_log2fc, "log2fc")