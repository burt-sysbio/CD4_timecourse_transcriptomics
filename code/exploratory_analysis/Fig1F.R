require(ggplot2)
require(tidyr)
require(dplyr)
require(tibble)
require(cowplot)
require(readr)
## make data

df1 <- read_csv("data/kinetic_genes/kinetic_genes_combined.csv")

df1 <- df1 %>% separate(gene, into = c("gene", "iso"), sep = "\\.")

kinetic_th0 <- unique(df1$gene[df1$kinetic_th0])
kinetic_th1 <- unique(df1$gene[df1$kinetic_th1])
kinetic_th2 <- unique(df1$gene[df1$kinetic_th2])
kinetic_thmix <- unique(df1$gene[df1$kinetic_thmix])



dat <- matrix(seq(1,16,1), ncol=4)
dat[upper.tri(dat)] <- NA

cellnames <- c("Th0", "Th1", "Th1/2", "Th2")
rownames(dat) <- cellnames
## reshape data (tidy/tall form)
dat2 <- dat %>%
  as_tibble() %>%
  rownames_to_column('Var1') 

dat2$Var1 <- cellnames
colnames(dat2) <- c("Var1", cellnames)

dat2 <- dat2 %>%
  gather(Var2, value, -Var1)

dat2 <- na.omit(dat2)

# this is a manual addition based on the combinations of var1 and var2
value2 <- c(sum(kinetic_th0 %in% kinetic_th0),
            sum(kinetic_th1 %in% kinetic_th0),
            sum(kinetic_thmix %in% kinetic_th0),
            sum(kinetic_th2 %in% kinetic_th0),
            sum(kinetic_th1 %in% kinetic_th1),
            sum(kinetic_thmix %in% kinetic_th1),
            sum(kinetic_th2 %in% kinetic_th1),
            sum(kinetic_thmix %in% kinetic_thmix),
            sum(kinetic_th2 %in% kinetic_thmix),
            sum(kinetic_th2 %in% kinetic_th2))


dat2$value2 <- value2
## plot data
ggplot(dat2, aes(Var1, Var2, fill = value2)) +
  geom_tile(color = "black") + 
  geom_text(aes(label = round(value2, 1)), size = 3.5) +
  scale_fill_gradient(low = "white", high = "red", name="shared genes") +
  theme_minimal(base_size = 14)+ 
  labs(x = "", y ="")+
  coord_fixed() 


ggsave("figures/fig2/tilemap_kinetic_genes.png", bg = "white", width = 4, height = 4)