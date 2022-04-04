# write function for manual clustering as it is done in masigpro
# make clustered heatmaps using masigpro for each cells kinetic genes
require(readr)
require(maSigPro)
require(tidyr)
require(dplyr)
require(pheatmap)
source("code/paper_theme.R")
source("code/utils.R")



# get correlation and tidy up
df <- read.csv("data/peine_data_log2fc_annot.csv")
df <- drop_genecol(df)

# get kinetic genes (combined criterion)
#res <- readRDS("code/fig2/mspro_tstep_contrast_Th1.RDS")
#sigs <- get.siggenes(res, rsq = 0.0, vars = "all")
#kinetic_genes <- sigs$summary

df_kin <- read.csv("data/kinetic_genes/kinetic_genes_combined.csv")
# if a gene is kinetic in any cell type, keep it
kinetic_genes <- df_kin$gene[apply(df_kin[,2:5], 1, any)]


df <- df[rownames(df) %in% kinetic_genes,]

df1 <- df[,1:10]
df2 <- df[,11:20]
df3 <- df[,21:30]
df4 <- df[31:40]

nclust <- 3

make_heatmap <- function(df, out, nclust, sname, cbar = F){

  # reorder df based on correlation clustering
  df <- df[out$order, ]

  annot <- as.data.frame(out$cut)
  colnames(annot) <- "c"
  annot$c <- as.character(annot$c)
  timepoints <- c(0,3,6,12,24,35,48,73,96,120)
  timepoints <- as.character(timepoints)
  
  # induce some gaps where clusters are
  # get the indices, convert to df for this, this counts number of observations
  test <- out$cut[out$order]
  idx <- which(diff(test)!=0)

  print("adding gaps")
  print(idx)
  pheatmap(df,
           scale = "none",
           breaks = seq(-2,2, length.out = 101),
           color = heatmap_colors,
           cluster_rows = F,
           cluster_cols = F,
           show_rownames = F,
           cellheight = 0.03,
           cellwidth = 12,
           #annotation_row = annot,
           labels_col = timepoints,
           gaps_row = idx,
           annotation_legend = F,
           legend = cbar,
           filename = paste0("figures/fig2/heatmaps_kinetic_clusters/", sname, "_", nclust,"clust.png"))
}

out1 <- clusterfun(df1, nclust)
out2 <- clusterfun(df2, nclust)
out3 <- clusterfun(df3, nclust)
out4 <- clusterfun(df4, nclust)

# plot the heatmaps
make_heatmap(df1,out1,nclust, "hm_Th0")
make_heatmap(df2,out2,nclust, "hm_Th1")
make_heatmap(df3,out3,nclust, "hm_Th2")
make_heatmap(df4,out4,nclust, "hm_Thmix")

annot1 <- as.data.frame(out1$cut)
annot2 <- as.data.frame(out2$cut)
annot3 <- as.data.frame(out3$cut)
annot4 <- as.data.frame(out4$cut)

annot1$gene <- rownames(annot1)
annot2$gene <- rownames(annot2)
annot3$gene <- rownames(annot3)
annot4$gene <- rownames(annot4)


# rename so stupid
colnames(annot1)[1] <- "annot"
colnames(annot2)[1] <- "annot"
colnames(annot3)[1] <- "annot"
colnames(annot4)[1] <- "annot"

annot1$celltype <- "Th0"
annot2$celltype <- "Th1"
annot3$celltype <- "Th2"
annot4$celltype <- "Th1/2"

# the cluster of Th1 is different, need to manually swap annotations
#clust2 <- annot2$annot == 2
#clust3 <- annot2$annot == 3
#annot2$annot[clust2] <- 3
#annot2$annot[clust3] <- 2
#
#clust2 <- annot4$annot == 2
#clust3 <- annot4$annot == 3
#annot4$annot[clust2] <- 3
#annot4$annot[clust3] <- 2


df <- bind_rows(annot1,annot2,annot3, annot4)
sname <- paste0("data/clustered_genes/kinetic_genes_20210714_",nclust, "clust.csv")
write.csv(df,sname, row.names = FALSE)