#### Notes ####

#### Libraries ####
library(ggplot2)
library(viridis)

#### Sourced functions ####
#source("../functions/annotate_boxplots.R")

#### Plotting functions ####
plot_triangle_correlation_matrix <- function(corr.matrix.diag, title="") {
  ggplot(data=corr.matrix.diag, 
         aes(x=Condition1, y=Condition2, fill=Spearman)) +
    geom_tile(color="white") +
    geom_text(aes(x=Condition1, y=Condition2, label=round(Spearman,2)), color = "white", size=4) +
    scale_fill_viridis(discrete=FALSE, option="viridis") +
    labs(title=title,
         x="Condition 1",
         y="Condition 2") +
    theme_bw() + 
    theme(plot.title=element_text(size=18),
          axis.title=element_text(size=14),
          axis.text=element_text(size=12),
          axis.text.x=element_text(angle=45, vjust=1, hjust=1),
          legend.title=element_text(size=14),
          legend.text=element_text(size=12),
          legend.position=c(0.1,0.8),
          panel.grid.major = element_blank()) +
    coord_fixed()  
} 
