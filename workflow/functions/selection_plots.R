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