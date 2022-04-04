require(ggplot2)
require(RColorBrewer)

update_geom_defaults("line", list(size = 1.0))


display.brewer.pal(n = 3, name = 'Set1')
colors2 <- brewer.pal(n=3, name = "Set1")
myblue <- colors2[2]
myred <- colors2[1]
colors <- c("darkgrey", myblue, "darkmagenta", myred)

set_colors <- function(x){
  if(x==T){
    options(ggplot2.discrete.colour = colors)
    options(ggplot2.discrete.fill = colors)    
  } else {
    options(ggplot2.discrete.colour = c("blue", "red3", "purple"))
    options(ggplot2.discrete.fill = c("blue", "red3", "purple"))
  }
}


mytheme <- function(base_size = 15){
  theme_bw(base_size) %+replace%
    theme(panel.grid.minor = element_blank(),
          legend.title = element_blank())
}

set_colors(x=T)

heatmap_colors <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100)
heatmap_colors2 <- colorRampPalette(colors = c("blue", "white", "red"))(100)