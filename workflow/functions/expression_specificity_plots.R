#### Notes ####

#### Libraries ####
library(dplyr)
library(ggplot2)

#### Sourced functions ####

#### Functions ####

# Plot tau distribution
plot_tau_distribution <- function(tau, prob, tau.low=NULL, tau.high=NULL, xlim=c(0.0,1.0), ylim=c(0.0,2.4)) {
  if (is.null(tau.low) & is.null(tau.high)) {
    quant <- quantile(tau, probs=prob)
    
    tau.low <- quant[2]
    tau.high <- quant[length(quant)-1]
  }
  
  den <- density(tau)
  
  plot(den, main="Expression specificity (tau) distribution", 
       cex.main=1.5, cex.axis=1.5, cex.lab=1.5, 
       xlim=xlim, ylim=ylim)
  polygon(c(den$x[den$x >= tau.high ], tau.high),
          c(den$y[den$x >= tau.high ], 0),
          col = "#F8766D", border=1)
  polygon(c(den$x[den$x <= tau.low ], tau.low),
          c(den$y[den$x <= tau.low ], 0),
          col = "#00BFC4", border=1)  
}

# Plot number of genes with top expression in organ/tissue
barplot_top_tissue <- function(expr.rank.first, apply.cutoff=FALSE, min.tau, title="") {
  if (apply.cutoff==FALSE) {
    top.tissue.count <- expr.rank.first %>%
      dplyr::count(PrimaryTissue) %>%
      dplyr::arrange(desc(n))
    
    ggplot(top.tissue.count, aes(x=reorder(PrimaryTissue, -n), y=n)) +
      geom_bar(stat="identity", fill="#F8766D") +
      geom_text(aes(label=n), vjust=-0.3, size=3.5) +
      labs(title=title,
           x="Primary Organ",
           y="Genes with top expression in organ") +
      theme_bw()
    
  } else {
    top.tissue.count <- expr.rank.first %>%
      dplyr::filter(tau > min.tau) %>%
      dplyr::count(PrimaryTissue) %>%
      dplyr::arrange(desc(n))
    
    ggplot(top.tissue.count, aes(x=reorder(PrimaryTissue, -n), y=n)) +
      geom_bar(stat="identity", fill="#F8766D") +
      geom_text(aes(label=n), vjust=-0.3, size=3.5) +
      labs(title=title,
           subtitle=paste("tau >", min.tau),
           x="Primary Organ",
           y="Genes with top expression in organ") +
      theme_bw()
  }
}