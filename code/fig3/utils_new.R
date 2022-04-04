# write function for manual clustering as it is done in masigpro
# make clustered heatmaps using masigpro for each cells kinetic genes
require(readr)
require(maSigPro)
require(tidyr)

clusterfun <- function(clusterdata, k){
  dcorrel <- matrix(rep(1, nrow(clusterdata)^2), 
                    nrow(clusterdata), nrow(clusterdata)) - 
    cor(t(clusterdata), use = "pairwise.complete.obs")
  clust <- hclust(as.dist(dcorrel), method = "ward.D")
  
  cut <- cutree(clust, k = k)
  
  OUTPUT <- list(cut, clust$order)
  names(OUTPUT) <- c("cut", "order")
  OUTPUT
}

tidy_df <- read_csv("data/peine_data_log2fc_annot.csv")

# get correlation and tidy up
df_cor <- readRDS("r_objects/fig2/correlation_df.RDS")


df <- df_cor[,1:40]
rownames(df) <- df_cor$Symbol

# get pvals from masigpro run, use both thm and th0 as contrast
res <- readRDS("code/fig2/mspro_tstep_contrast_Th1.RDS")
sigs <- get.siggenes(res, rsq = 0.6, vars = "all")
kinetic_genes <- sigs$summary

df <- df[rownames(df) %in% kinetic_genes,]

# check if caro genes are in here

caro_genes <- read_csv("genes_caro.csv", col_names = "gene")
missing <- caro_genes$gene[!(caro_genes$gene %in% df_cor$Symbol)]
print(sum(rownames(df) %in% caro_genes$gene))

load("data/raw_data/4_ready_to_use/allTh_names_Rep1.Rdata")