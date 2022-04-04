require(maSigPro)
require(ggplot2)
require(dplyr)
require(readr)

# get DEGs for Th1 vs Th2 then store
res <- readRDS("code/fig2/mspro_tstep_contrast_Th1.RDS")
# choose all to get all significant genes
# choose "groups" for th1 th2 degs

sigs <- get.siggenes(res, rsq = 0.7, vars = "all")
siggenes <- sigs$sig.genes
#siggenes <- sigs$sig.genes$Th2vsTh1

k <- 6
# use different clustering metrics (raw, coeffs, tscores for cluster.data) pass to hclust
out1 <- see.genes(siggenes, k=k, cluster.data = 1)
out2 <- see.genes(siggenes, k=k, cluster.data = 2)
out3 <- see.genes(siggenes, k=k, cluster.data = 3)

# these are cluster assignments for all genes with rsq < 0.6
df_clust <- data.frame(c1 = out1$cut, c2 = out2$cut, c3 = out3$cut)
df_clust$gene <- rownames(df_clust)

sname <- paste0("code/fig2/peine_kinetic_deg_", k, "_clusters")
#saveRDS(df_clust, paste0(sname, ".RDS"))
write.csv(df_clust, paste0(sname, ".csv"))
