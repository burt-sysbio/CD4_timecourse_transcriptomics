
#prepare data frame for class switches
require(dplyr)
require(tidyr)

nclust <- 3

df <- read.csv("data/clustered_genes/kinetic_genes_20210714_3clust.csv")

mydate <- "20210825"
# check for class switches
df2 <- df %>% pivot_wider(names_from = celltype, values_from = annot)

#
df2$s_th1_th2 <- df2$Th1 != df2$Th2
df2$s_th1_thm <- df2$Th1 != df2$`Th1/2`
df2$s_th2_thm <- df2$Th2 != df2$`Th1/2`

df2$s_th1_th0 <- df2$Th1 != df2$Th0
df2$s_th2_th0 <- df2$Th2 != df2$Th0
df2$s_th0_thm <- df2$Th0 != df2$`Th1/2`

# kick out genes for which there are no class switches
df3 <- df2[rowSums(df2[,6:ncol(df2)],na.rm = T) > 0,]


# get switches from clust 1-->3
df3$switch1_3 <- apply(df3[,2:5] == 1, 1, any) & apply(df3[,2:5] == 3, 1, any)

# switches from clust 2 --> 3
df3$switch2_3 <- apply(df3[,2:5] == 2, 1, any) & apply(df3[,2:5] == 3, 1, any)

# switches from clust 1-- > 2
df3$switch1_2 <- apply(df3[,2:5] == 1, 1, any) & apply(df3[,2:5] == 2, 1, any)

sname <- paste0("data/classswitch_quantDEGS/class_switch_summary_", mydate, "_", nclust, "clust.csv")
write.csv(df3, sname, row.names = F)



