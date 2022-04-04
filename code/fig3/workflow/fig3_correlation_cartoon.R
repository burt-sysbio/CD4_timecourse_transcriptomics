# example plots that correlation works

require(readr)
require(dplyr)
require(tidyr)
require(ggplot2)
require(cowplot)
require(svglite)
require(stringr)
require(maSigPro)
require(RColorBrewer)

colors2 <- brewer.pal(n=3, name = "Set1")
myblue <- colors2[2]
myred <- colors2[1]
colors <- c(myblue, myred)

options(ggplot2.discrete.colour = colors)
options(ggplot2.discrete.fill = colors)

# plot timecourse for high correl and low correl genes
plot_timecourse <- function(df_tidy, savename, textsize, labelsize,
                            width = 4, height = 2){
  # plot time course for high/low correlations
  p1 <- ggplot(df_tidy, aes(time, value))
  p2 <- p1 + 
    geom_line(aes(colour = celltype), size = 0.2, show.legend = F)+
    geom_point(aes(colour = celltype), show.legend = F, size = 0.5) +
    theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title = element_text(size = labelsize),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          plot.margin=unit(c(0.2,0.5,0.1,0.1),"cm")) +
    ylab("expr. norm.") +
    xlab("time (h)")
  
  if(grepl("cartoons", savename)){
    p2 <- p2 + 
      scale_x_continuous(name = "", breaks = c(0,60,120)) +
      scale_y_continuous(name = "", breaks = c(-2,0,2)) + 
      theme_bw(base_size = 10) +
      theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    
    width = 1.2
    height = 1.0
    
    savename2 <- str_sub(savename, 1, -4)
    savename2 <- paste0(savename2, "svg")
    ggsave(filename = savename2, width = width, height = height, bg = "white")
  } else {
    p2 <- p2 + theme_bw(base_size = textsize) 
    
  }
  
  ggsave(filename = savename, width = width, height = height, bg = "white")
  return(p2)
}

# plot Th1 vs Th2 with time as color code
plot_correlation <- function(df_tidy, savename, textsize, labelsize,
                             width = 4, height = 2){
  df_wide <- pivot_wider(df_tidy, names_from = celltype, values_from =  value)
  p2 <-ggplot(df_wide, aes(Th1, Th2))
  p3 <- p2 + 
    geom_smooth(method = "lm", color= "grey", se = F, linetype = "dashed", size = 0.2) +
    xlab("expr. Th1") +
    ylab("expr. Th2") +
    theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title = element_text(size = labelsize),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          plot.margin=unit(c(0.2,0.5,0.1,0.1),"cm"))

  if(grepl("cartoons", savename)){
    
    p3 <- p3 + 
      geom_point(size = 0.5, colour = "black") +
      theme_bw(base_size = 10) +
      theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    
      width = 1.2
      height = 1.0
      savename2 <- str_sub(savename, 1, -4)
      savename2 <- paste0(savename2, "svg")
      ggsave(filename = savename2, width = width, height = height, bg = "white", device = "svg")
  } else {
    p3 <- p3 + 
      geom_point(aes(color = time))+
      theme_bw(base_size = textsize) +
      scale_colour_gradient(low = "grey", high = "black")
  }
    
  ggsave(filename = savename, width = width, height = height, bg = "white")
  
  return(p3)
}

# subdirectory to save figures to
savedir <- "figures/fig3/correlation_examples/"
# load timecourse data for th1 and th2
df_all <- read.csv("data/peine_logfc_tidy.csv")
df_all <- df_all %>% filter(celltype %in% c("Th1", "Th2"))

df_cor <- read.csv("data/correlation/correlation_df.csv")
df_cor <- df_cor %>% select(gene, cor_th1_th2)
df_all <- left_join(df_all, df_cor)


# get correlation add to tidy data
df_cor <- df_cor %>% pivot_longer(-gene, values_to = "value", names_to = "correlation")

# get DEGS from masig run
res <- readRDS("data/masigpro_output/tfit_fullrun_kinbackground_a0.01")

rsq <- 0.6
sigs <- get.siggenes(res, rsq = rsq, vars = "groups")
sigs <- sigs$summary$Th2vsTh1
candidates_quant_deg <- sigs[sigs!= " "]


# get qualitative DEGS
degs <- read_csv("genesets_output/quantDEGS/degs_Th1_Th2_cor0.5.csv")
df_qual_deg <- df_all %>% filter(gene %in% degs$gene)
df_quant_deg <- df_all %>% filter(gene %in% candidates_quant_deg)


# I only know my candidates and its only a cartoon so it dont matter
# just reassign instead of filter

# pick candidate
n <- 20
candidates_qual_deg <- df_qual_deg %>% select(gene, cor_th1_th2) %>% unique() %>% 
  arrange(cor_th1_th2) %>% slice_head(n=n)

candidates_quant_deg <- df_quant_deg %>% select(gene, cor_th1_th2) %>% unique() %>% 
  arrange(desc(cor_th1_th2)) %>% slice_head(n=n) 


candidates_qual_deg <- candidates_qual_deg$gene
candidates_quant_deg <- candidates_quant_deg$gene

#candidates_qual_deg <- candidates_qual_deg[1:10,1]
textsize = 12
labelsize = 11
width1 = 2.3
width2 = 3.1
height = 1.8

pipeline <- function(candidate, df, deg_type, textsize = 12, labelsize = 11, width =2.3, width2 = 3.1, height = 1.8){
  stopifnot(candidate %in% df$gene)
  df <- df[df$gene == candidate,]
  # output names
  savedir <- paste0("figures/fig3/correlation_examples/", deg_type, "/")
  savename1 <- paste0(savedir, candidate, "_timecourse.png")
  savename2 <- paste0(savedir, candidate, "_correlation.png")

  p1 <- plot_timecourse(df, savename = savename1, textsize = textsize, 
                        labelsize = labelsize, width = width1, height = height)
  
  p2 <- plot_correlation(df, savename = savename2, textsize = textsize, 
                         labelsize = labelsize, width = width2, height = height)
}

#lapply(candidates_qual_deg, pipeline, df_qual_deg, "qualitative_DEG")
#lapply(candidates_quant_deg, pipeline, df_quant_deg, "quantitative_DEG")
lapply(candidates_quant_deg, pipeline, df_quant_deg, "cartoons_quant_deg")
lapply(candidates_qual_deg, pipeline, df_qual_deg, "cartoons_qual_deg")


