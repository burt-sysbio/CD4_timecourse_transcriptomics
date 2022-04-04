# plot venn diagram for two given lists
require(venn)
require(dplyr)
require(tidyr)
require(svglite)
require(readxl)
require(ggpubr)

source("code/paper_theme.R")

plot_venn <- function(x1, x2, contrast, cor_thres){
  
  x1 <- na.omit(x1)
  x2 <- na.omit(x2)
  #mycolors <- c("blue", "red", "green")
  sname <- paste0("figures/fig3/venn_quantDEGS/venn_quantDEGs_", contrast, "_cor", cor_thres, ".svg")
  
  if(contrast == "ThMix"){
    names = c("Th1vsTh1/2", "Th2vsTh1/2")
  } else {
    names = c("Th1vsTh0", "Th2vsTh0")
  }
  
  svg(filename = sname, width = 3, height = 3)
  venn(x=  list(x1, x2), borders = FALSE, box = F, ilcs = 1.0, snames = names, lty = 1)
  
  # add some annotations?
  # x <- c(50, 300)
  # y <- c(1,1)
  # mylabels <- c("a","b")
  # mytext <- paste(mylabels, c(1,2))
  # text(x, y, mytext, cex = 0.5)
  dev.off()
}

col_types <- c("text", "numeric", "numeric", "logical", "logical", "numeric", "numeric", "logical")
df1 <- read_xlsx("figures/supp/table_S1_th1_thmix.xlsx", col_types = col_types)
df2 <- read_xlsx("figures/supp/table_S1_th2_thmix.xlsx", col_types = col_types)
df3 <- read_xlsx("figures/supp/table_S1_th1_th2.xlsx", col_types = col_types)
df4 <- read_xlsx("figures/supp/table_S1_th1_th0.xlsx", col_types = col_types)
df5 <- read_xlsx("figures/supp/table_S1_th2_th0.xlsx", col_types = col_types)


myfun <- function(df, comparison){
  qual_deg <- length(df$gene[df$qual_deg])
  quant_deg <- length(df$gene[df$quant_deg])
  
  qual_switch <- length(df$gene[df$qual_deg & df$switch])
  quant_switch <- length(df$gene[df$quant_deg & df$switch])
  
  out1 <- paste0(quant_deg, " (", quant_switch, ")")
  out2 <- paste0(qual_deg, " (", qual_switch, ")")
  out <- data.frame(comparison = comparison, quant_deg = out1, qual_deg = out2)
  return(out)
}

qual_th1_th12 <- df1$gene[df1$qual_deg]
qual_th2_th12 <- df2$gene[df2$qual_deg]

quant_th1_th12 <- df1$gene[df1$quant_deg]
quant_th2_th12 <- df2$gene[df2$quant_deg]

plot_venn(qual_th1_th12,qual_th2_th12, contrast = "ThMix", "0.3")
plot_venn(qual_th1_th12,qual_th2_th12, contrast = "Th0", "0.3")

plot_venn(quant_th1_th12,quant_th2_th12, contrast = "ThMix", "0")
plot_venn(quant_th1_th12,quant_th2_th12, contrast = "Th0", "0")


out1 <- myfun(df1, "Th1vsTh1/2")
out2 <- myfun(df2, "Th2vsTh1/2")
out3 <- myfun(df3, "Th1vsTh2")
out4 <- myfun(df4, "Th1vsTh0")
out5 <- myfun(df5, "Th2vsTh0")

out <- bind_rows(out3, out1, out4, out2, out5)
colnames(out) <- c("comparison", "quant. DEG", "qual. DEG")
# print the cell numbers

tbody.style = tbody_style(color = "black",
                          fill = c("grey", "white"), hjust = 0, x = 0.1)
mytable <- ggtexttable(out, rows = NULL,theme = ttheme(
  colnames.style = colnames_style(color = c("black", "black", myred), fill = "white"),
  tbody.style = tbody.style
))
mytable <- mytable %>%
  tab_add_hline(at.row = c(1, 2), row.side = "top", linewidth = 3, linetype = 1) %>%
  tab_add_hline(at.row = c(6), row.side = "bottom", linewidth = 3, linetype = 1)

ggsave("figures/fig3/deg_summary_table.png", bg = "white")