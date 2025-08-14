#### Notes ####

#### Libraries ####
library(ggplot2)
library(ggpubr)

#### Sourced functions ####

#### Functions ####

# Annotate boxplot with number of observations per grouping
annotate_n <- function(x) {
  return(c(y = -0.01, label = length(x))) 
}
# Annotate boxplot with median per grouping
annotate_median <- function(x) {
  return(c(y = median(x)*1.20, label = round(median(x),2))) 
}

# Strength of purifying selection
## Boxplot with three variability categories
boxplot_omega0 <- function(selectome.df, compare.means=FALSE, title="") {
  if (compare.means==TRUE) {
    comparisons <-
      list(c("low","mid"), c("mid","high"), c("low","high"))
    
    ggplot(selectome.df, aes(x=EV, y=omega0, fill=EV)) +
      stat_boxplot(geom='errorbar') +
      geom_boxplot(notch=TRUE) +
      #ylim(c(0,0.60)) +
      ylim(c(-0.01,0.50)) +
      stat_compare_means(comparisons=comparisons, 
                         method="wilcox.test", 
                         label="p.signif",
                         #label.y=c(0.425,0.475,0.525)) +
                         label.y=c(0.325,0.375,0.425),
                         size=7) +
      stat_summary(fun.data= annotate_n, geom="text", fun=median, 
                   position= position_dodge(width = 0.75), size=7) +      
      stat_summary(fun.data= annotate_median, geom="text", fun=median, 
                   position= position_dodge(width = 0.75), size=7) +
      labs(title=title,
           x="Variability rank",
           y="ω0 (dN/dS < 1)",
           caption=paste0("produced on ", Sys.time())) +
      guides(fill=guide_legend(title="Variability rank")) +
      theme_bw() +
      theme(plot.title=element_text(size=24, face="italic"),
            plot.subtitle=element_text(size=20),
            axis.title=element_text(size=20),
            axis.text=element_text(size=20),
            legend.position="none")
  } else {
    ggplot(selectome.df, aes(x=EV, y=omega0, fill=EV)) +
      stat_boxplot(geom='errorbar') +
      geom_boxplot(notch=TRUE) +
      #ylim(c(0,0.60)) +
      ylim(c(-0.01,0.50)) +
      stat_summary(fun.data= annotate_n, geom="text", fun=median, 
                   position= position_dodge(width = 0.75), size=7) +
      stat_summary(fun.data= annotate_median, geom="text", fun=median, 
                   position= position_dodge(width = 0.75), size=7) +
      labs(title=title,
           x="Variability rank",
           y="ω0 (dN/dS < 1)",
           caption=paste0("produced on ", Sys.time())) +
      guides(fill=guide_legend(title="Variability rank")) +
      theme_bw() +
      theme(plot.title=element_text(size=24, face="italic"),
            plot.subtitle=element_text(size=20),
            axis.title=element_text(size=20),
            axis.text=element_text(size=20),
            legend.position="none")
  }
}

## Boxplot with variability rank bins
boxplot_omega0_varbins <- function(selectome.df, compare.extremes=FALSE, title="") {
  if (compare.extremes==TRUE) {
    # Test for a significant difference between the two extremes (lowly and highly variable)
    lv.hv.wilcox <- wilcox.test(
      omega0 ~ VarRankCategory,
      data = subset(selectome.df, 
                    VarRankCategory %in% c("low", "high")),
      exact = FALSE,
    )
    
    # Segments are drawn assuming 10 bins
    ggplot(selectome.df,
           aes(x=as.factor(VarRankBin), y=omega0, fill=as.factor(VarRankBin))) +
      stat_boxplot(geom ='errorbar') +
      geom_boxplot(notch=TRUE) +
      stat_summary(fun.data= annotate_n, geom="text", fun=median, 
                   position= position_dodge(width = 0.75), size=7) +      
      stat_summary(fun.data=annotate_median, geom="text", fun=median, 
                   position=position_dodge(width=0.75), size=6) +
      ylim(c(-0.01,0.50)) +
      # Long horizontal segment
      annotate("segment",
               x=1.5, xend=9.5,
               y=0.35, yend=0.35,
               size=0.6) +
      # Vertical segments
      annotate("segment",
               x=1.5, xend=1.5,
               y=0.35, yend=0.33,
               size = 0.6) +
      annotate("segment",
               x=9.5, xend=9.5,
               y=0.35, yend=0.33,
               size=0.6) +
      # Short horizontal segments
      annotate("segment",
               x=1.0, xend=2.0,
               y=0.33, yend=0.33,
               size=0.6) +
      annotate("segment",
               x=9.0, xend=10.0,
               y=0.33, yend=0.33,
               size=0.6) +
      # Results of Wilcoxon test
      annotate("text",
               x=6,y=0.35*1.05,
               label=bquote("Wilcoxon " * italic(p) == .(formatC(lv.hv.wilcox$p.value, format = "g", digits = 3))),
               size=6) +
      # Overall Spearman correlation
      stat_cor(data=selectome.bgd.groupby.gene.ev.median,
               aes(x=MedianVar, y=omega0),
               method="spearman",
               inherit.aes=FALSE,
               label.x=1,   
               label.y=0.4,
               size=6) +
      labs(title=title,
           x="Median variability rank bin",
           y="ω0 (dN/dS < 1)",
           caption=paste0("produced on ", Sys.time())) +
      theme_bw() +
      theme(plot.title=element_text(size=28, face="italic"),
            plot.subtitle=element_text(size=24),
            axis.title=element_text(size=24),
            axis.text=element_text(size=18),
            strip.text=element_text(size=24),
            legend.position="none")
    
  } else {
    ggplot(selectome.df,
           aes(x=as.factor(VarRankBin), y=omega0, fill=as.factor(VarRankBin))) +
      stat_boxplot(geom ='errorbar') +
      geom_boxplot(notch=TRUE) +
      stat_summary(fun.data= annotate_n, geom="text", fun=median, 
                   position= position_dodge(width = 0.75), size=7) +      
      stat_summary(fun.data=annotate_median, geom="text", fun=median, 
                   position=position_dodge(width=0.75), size=6) +
      ylim(c(-0.01,0.50)) +
      # Overall Spearman correlation
      stat_cor(data=selectome.bgd.groupby.gene.ev.median,
               aes(x=MedianVar, y=omega0),
               method="spearman",
               inherit.aes=FALSE,
               label.x=1,   
               label.y=0.4,
               size=6) +
      labs(title=title,
           x="Median variability rank bin",
           y="ω0 (dN/dS < 1)",
           caption=paste0("produced on ", Sys.time())) +      
      theme_bw() +
      theme(plot.title=element_text(size=28, face="italic"),
            plot.subtitle=element_text(size=24),
            axis.title=element_text(size=24),
            axis.text=element_text(size=18),
            strip.text=element_text(size=24),
            legend.position="none")    
  }
}

## Scatterplot of purifying selection as a function of median variability across conditions
scatterplot_omega0_variability <- function(selectome.df, title="") {
  ggplot(selectome.df, aes(x=MedianVar, y=omega0)) + 
    geom_point(color="#440154", size=0.2, alpha=0.6) +
    stat_density_2d(aes(fill=..level..), geom="polygon") +
    scale_fill_continuous(type="viridis") +
    stat_smooth(method="loess",
                span=0.60,
                formula=y~x,
                se=FALSE) +
    stat_cor(method="spearman", size=8) +
    labs(title=title,
         x="Variability rank",
         y="ω0 (dN/dS < 1)",
         caption=paste0("produced on ", Sys.time())) +
    theme_bw() +
    theme(plot.title=element_text(size=24, face="italic"),
          plot.subtitle=element_text(size=20),
          plot.caption=element_text(size=12),
          axis.title=element_text(size=20),
          axis.text=element_text(size=20),
          legend.position="none")
}