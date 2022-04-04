
require(dplyr)

# load replicate 1 data
load("data/raw_data/4_ready_to_use/allTh_Rep1.Rdata")
load("data/raw_data/4_ready_to_use/allTh_names_Rep1.Rdata")

# remove lowly/ no-express genes
threshold <- median(allTh_1)
rowmax <- apply(allTh_1, 1, max)
allTh_1 <- allTh_1[which(rowmax >= threshold),]
dim(allTh_1) #18284 

# geneNames for probes IDs
geneIds <- readRDS("data/annotation/burt_annot.rds")


# annotate df with genes names
allTh_1 <- as.data.frame(allTh_1)
allTh_1$id <- rownames(allTh_1)
geneIds$id <- rownames(geneIds)
df <- inner_join(allTh_1, geneIds, by = "id")
rownames(df) <- df$Symbol

df <- df[,1:40]

#write.csv(df, "data/peine_data_median_threshold_annot.csv")

df_cleaned <- df
df_cleaned$mymax <- rowSums(df_cleaned)

df_cleaned$gene <- rownames(df)
df_cleaned <- df_cleaned %>% separate(gene, into = c("gene", "iso"), sep ="\\.")
df_cleaned <- df_cleaned %>% group_by(gene) %>% filter(mymax == max(mymax))
df_cleaned <- as.data.frame(df_cleaned)
rownames(df_cleaned) <- df_cleaned$gene


df_log <- as.matrix(df_cleaned[,1:40])
df_log <- log2(df_log / df_log[,1])
df_log <- as.data.frame(df_log)
df_log$gene <- rownames(df_log)


write.csv(df_log, "data/cleaned_isoforms/data_log2FC.csv", row.names = F)


df_log_tidy <- df_log %>% pivot_longer(-gene)
df_log_tidy <- df_log_tidy %>% separate(name, into = c("celltype", "time", NA), sep = "_")
df_log_tidy$time <- as.numeric(df_log_tidy$time)

write.csv(df_log_tidy, "data/cleaned_isoforms/data_log2FC_tidy.csv", row.names = F)
