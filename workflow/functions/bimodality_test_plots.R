#### Notes ####

#### Libraries ####
library(ggplot2)
library(ggpubr)

#### Sourced functions ####
source("../functions/annotate_boxplots.R")

#### Functions ####
boxplot_bimodality_by_variability <- function(bimodality.test.bind, title="", subtitle="", facet=TRUE) {
  if (facet==TRUE) {
    ggplot(bimodality.test.bind,
           aes(x=as.factor(VarRankBin), y=BI, fill=as.factor(VarRankBin))) +
      stat_boxplot(geom ='errorbar') +
      geom_boxplot(notch=TRUE) +
      stat_summary(fun.data=annotate_median, geom="text", fun=median, 
                   position=position_dodge(width=0.75), size=6) +
      labs(title=title,
           subtitle=subtitle,
           x="Variability rank bin",
           y="Bimodality index",
           fill="Variability rank bin") +
      facet_wrap(vars(Condition), nrow=3) +
      theme_bw() +
      theme(plot.title=element_text(size=28, face="italic"),
            plot.subtitle=element_text(size=24),
            axis.title=element_text(size=24),
            axis.text=element_text(size=18),
            strip.text=element_text(size=24),
            legend.position="none")
    
  } else {
    ggplot(bimodality.test.bind,
           aes(x=as.factor(VarRankBin), y=BI, fill=as.factor(VarRankBin))) +
      stat_boxplot(geom ='errorbar') +
      geom_boxplot(notch=TRUE) +
      stat_summary(fun.data=annotate_median, geom="text", fun=median, 
                   position=position_dodge(width=0.75), size=6) +
      labs(title=title,
           subtitle=subtitle,
           x="Variability rank bin",
           y="Bimodality index",
           fill="Variability rank bin") +
      theme_bw() +
      theme(plot.title=element_text(size=28, face="italic"),
            plot.subtitle=element_text(size=24),
            axis.title=element_text(size=24),
            axis.text=element_text(size=18),
            strip.text=element_text(size=24),
            legend.position="none")       
  }
}

boxplot_bimodality_by_variability_v2 <- function(bimodality.test.bind, title="", subtitle="") {
  ggplot(bimodality.test.bind,
         aes(x=as.factor(VarRankBin), y=BI)) +
    stat_boxplot(geom ='errorbar') +
    geom_boxplot(notch=TRUE, fill="#00BFC4") +
    stat_summary(fun.data = annotate_median, geom = "text", fun = median, 
                 position = position_dodge(width = 0.75), size=7) +
    labs(title=title,
         subtitle=subtitle,
         x="Variability rank bin",
         y="Bimodality index",
         fill="Variability rank bin") +
    ylim(c(0,25)) +
    theme_bw() +
    theme(plot.title=element_text(size=28, face="italic"),
          plot.subtitle=element_text(size=24),
          axis.title=element_text(size=28),
          axis.text=element_text(size=24),
          legend.position="none")  
}

boxplot_bimodality_by_variability_comparison <- function(bimodality.test.bind, title="") {
  ggplot(bimodality.test.bind,
         aes(x=as.factor(VarRankBin), y=log2(BI), fill=as.factor(Dataset))) +
    stat_boxplot(geom ='errorbar') +
    geom_boxplot(notch=TRUE) +
    stat_compare_means(aes(label=..p.signif..),
                       method="wilcox.test",
                       size=8) +
    #stat_summary(fun.data=annotate_median, geom="text", fun=median, 
    #             position=position_dodge(width=0.75), size=6) +
    labs(title=title,
         x="Variability rank bin",
         y="log(BI)",
         fill="") +
    theme_bw() +
    theme(plot.title=element_text(size=28, face="italic"),
          axis.title=element_text(size=28),
          axis.text=element_text(size=24),
          legend.text=element_text(size=24),
          legend.key.size = unit(2,'cm'),
          legend.position="bottom") 
}
  
plot_density_zscore <- function(zscore.matrix, bimodality.test, title="") {
  par(mfrow=c(3,4), mai = c(0.8,0.6,0.2,0.2))
  
  b <- 0
  while (b < 1.0) {
    i <- 0.1
    
    zscore.vector <-
      as.vector(zscore.matrix[bimodality.test$GeneID[(bimodality.test$VarRankBin==round(b,1))],])
    
    plot(density(zscore.vector),
         main="",
         #cex.lab=1.5, cex.axis=1.5,
         ylab=paste("rank = ", format(b, nsmall=1)),
         xlab=paste("N =", length(zscore.vector)/ncol(zscore.matrix), "x", ncol(zscore.matrix)))
    
    b <- b + i
  }
  
  mtext(title, side=3, line=-1, outer=TRUE)
}

plot_density_zscore_by_bi <- function(zscore.matrix, bimodality.test, title="") {
  par(mfrow=c(3,4), mai = c(0.8,0.6,0.2,0.2))
  
  b <- 0
  while (b < 6.0) {
    i <- 0.5
    b <- b + i
    
    print(format(b,nsmall=1))
    zscore.vector.tmp <- 
      as.vector(zscore.matrix[bimodality.test$GeneID[
        (bimodality.test$BIBin > round(b-i,1) & bimodality.test$BIBin <= round(b,1))],])
    
    plot(density(zscore.vector.tmp),
         main="",
         ylab=paste(format(b-i, nsmall=1), "< BI", expression("\u2264"), format(b, nsmall=1)),
         xlab=paste("N =", length(zscore.vector.tmp)/ncol(zscore.matrix), "x", ncol(zscore.matrix)))
    
  }  
}