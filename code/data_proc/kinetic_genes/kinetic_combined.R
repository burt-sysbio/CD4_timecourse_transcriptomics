# load both approaches how to define kinetic genes and combine
require(readr)
require(dplyr)
require(tidyr)

df1 <- read_csv("data/kinetic_genes/kinetic_genes_2_timepoints_logfc1.csv")
df2 <- read_csv("data/kinetic_genes/kinetic_genes_masigpro_summary.csv")

df1 <- pivot_longer(df1, -gene, values_to = "val1")
df2 <- pivot_longer(df2, -gene, values_to = "val2")

df <- inner_join(df1, df2)
# check for each row if there is at least one True value
df$val3 <- apply(df[,2:3]== TRUE, 1, any)

# only keep last entry this is the kinetic definition: either masig or 2 timepoitns

df <- df %>% select(gene, name, val3)
df <- df %>% pivot_wider(names_from = name, values_from = val3)

write.csv(df, "data/kinetic_genes/kinetic_genes_combined.csv", row.names = F)