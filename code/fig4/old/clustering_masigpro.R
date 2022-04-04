
# make clustered heatmaps using masigpro for each cells kinetic genes
require(readr)
require(pheatmap)
require(maSigPro)
source("code/fig3/utils.R")
require(tidyr)

# load Th1 masigpro fit
cell = "Th1"
Q = "0.05"
a = "0.05"

stopifnot(cell %in% c("Th1", "Th2", "ThMix", "Th0"))

tfit <- readRDS(paste0("r_objects/msgpro_", cell,"_Q", Q, "a", a, ".RDS"))

# extracting significat profiles differences
allSigs <- get.siggenes(tfit, rsq = 0.6, vars = "all")

kin_th1 <- allSigs$summary
# clustering
n_clust <- 4
myclust <- msclust(allSigs$sig.genes, 
                        show.fit = T, 
                        dis =tfit$dis,
                        cluster.method="hclust",
                        cluster.data = 1, 
                        k = n_clust, 
                        newX11=FALSE)


df <- tfit$dat

df <- as.data.frame(df)
df <- df[rownames(df) %in% kin_th1,]
annot <- myclust$cut
annot <- as.data.frame(annot)


#df <- df %>% arrange(annot) %>% select(-annot)

timepoints <- c(0,3,6,12,24,35,48,73,96,120)
timepoints <- as.character(timepoints)
order <- myclust$order
df <- df[order,]
pheatmap(df, 
         scale = "row",
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = F,
         labels_col = timepoints,
         annotation_row = annot,
         width = 3,
         height = 4,
         filename = paste0("figures/fig3/heatmap_",cell,".png"))

#
df <- log2(df/df[,1])
df$id <- rownames(df)

df <- df %>% pivot_longer(-id, names_to = "time", values_to = "expr")
df <- df %>% separate(time, into = c(NA, "time", NA), sep ="_")
df$time <- as.numeric(df$time)

annot$id <- rownames(annot)
df <- left_join(df, annot)

df <- df %>% group_by(time, annot) %>% mutate(avg = mean(expr), sd = sd(expr)) %>% ungroup()
df <- df %>% mutate(upper = avg+sd, lower = avg-sd)

df$annot <- as.factor(df$annot)

saveRDS(df, paste0("r_objects/fig3/timecourse_cluster", cell, ".RDS"))