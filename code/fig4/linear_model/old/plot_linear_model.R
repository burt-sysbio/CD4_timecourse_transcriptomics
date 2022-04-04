# plot results from linear model
require(tidyr)
require(ggplot2)
require(dplyr)
source("code/paper_theme.R")
df <- read.csv("data/peine_logfc_tidy.csv")
df_out <- read.csv("data/linear_model.csv")


limits <- c(0,1)
df_out$clust <- as.character(df_out$clust)
g <- ggplot(data = df_out, aes(x = Th1, y = Th2, colour = clust))
g1 <- g + geom_point() + theme_bw() + facet_wrap(~model) +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_hline(yintercept=0, linetype="dashed")
print(g1)


# get top th2 genes
df_out2 <- df_out %>% filter(model == "Th0") # %>%
  #filter((Th2 > 0.5) & (Th1 > 0.5))



df_red <- inner_join(df, df_out2[c("gene", "Th1", "Th2", "clust")], by = c("gene"))



annot <- df_red[c("gene", "Th1", "Th2")]
annot <- unique(annot)
annot <- annot %>% mutate(title = paste("a", round(Th1,2), "b", round(Th2,2)))
annot$x <- 50
annot$y <- 5

g <- ggplot(data = df_red, aes(x = time, y = value))
g1 <- g + geom_line(aes(colour = celltype, linetype = gene)) + facet_wrap(~clust)# + 
  #geom_label(data = annot, aes(x=x, y =y, label = title))
print(g1)

