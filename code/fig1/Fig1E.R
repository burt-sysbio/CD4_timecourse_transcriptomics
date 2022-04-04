
# plot  PCA
require(readr)
require(dplyr)
require(ggplot2)
require(stringr)


df <- read_csv("data/peine_data_log2fc_annot.csv")


df <- df %>% select(-gene)

mat <- t(df)
pca <- prcomp(mat)

out <- pca$x
out <- as.data.frame(out)

cells <- rep(c("Th0","Th1", "Th2", "Th1/2"), each = 10)
time <- c(0,3,6,12,24,35,48,73,96,120)
#time <- seq(1:10)
time <- rep(time, 4)

out$time <- time
out$cell <- cells

# get PC var percentage
var <- summary(pca)$importance
var <- as.data.frame(var)
var1 <- round(var$PC1[2]*100,1)
var2 <- round(var$PC2[2]*100,1)

#
g <- ggplot(data = out, aes(PC1, PC2, color = time, shape = cell))
g + geom_point(size = 2.5) + 
  theme_bw(base_size = 18) +
  xlab(paste0("PC1 (", var1, "%)"))+
  ylab(paste0("PC2 (", var2, "%)"))

ggsave("figures/fig1E.png", width = 5, height = 3.3)