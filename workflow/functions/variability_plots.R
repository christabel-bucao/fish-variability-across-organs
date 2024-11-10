#### Notes ####

#### Libraries ####
library(ggplot2)
library(ggpubr)
library(cowplot)

#### Sourced functions ####

#### Functions ####
plot_sd_vs_mean <- function(summary.stats, metadata, log2=FALSE) {
  Tissue <- as.factor(metadata$Tissue)
  
  if ("Sex" %in% colnames(metadata)) {
    Sex <- as.factor(metadata$Sex)
    gg <- vector(mode="list", length(levels(Tissue))*length(levels(Sex)))
    names(gg) <- levels(interaction(Tissue, Sex, sep="_", lex.order=TRUE))
    
    for (t in levels(Tissue)) {
      for (s in levels(Sex)) {
        if (nrow(summary.stats[[t]][[s]]>0)) {
          if (log2==FALSE) {
            gg[[paste(t,s,sep="_")]] <-
              ggplot(summary.stats[[t]][[s]], aes(x=Mean, y=SD)) + 
              geom_point(color="#440154", size=0.2, alpha=0.6) +
              stat_density_2d(aes(fill=..level..), geom="polygon") +
              scale_fill_continuous(type = "viridis") +
              stat_smooth(method="loess",
                          span=0.60,
                          formula=y~x,
                          se=FALSE) +
              stat_cor(method="pearson") +
              labs(title=gsub("_"," ",t), 
                   subtitle=gsub("_"," ", s)) +
              theme_bw()
          } else {
            gg[[paste(t,s,sep="_")]] <-
              ggplot(summary.stats[[t]][[s]], aes(x=Mean, y=log2(SD))) + 
              geom_point(color="#440154", size=0.2, alpha=0.6) +
              stat_density_2d(aes(fill=..level..), geom="polygon") +
              scale_fill_continuous(type = "viridis") +
              stat_smooth(method="loess",
                          span=0.60,
                          formula=y~x,
                          se=FALSE) +
              stat_cor(method="pearson") +
              labs(title=gsub("_"," ",t), 
                   subtitle=gsub("_"," ", s)) +
              theme_bw()          
          }
        } else {
          gg[[paste(t,s,sep="_")]] <- NULL
        }
      }
    }
  } else {
    gg <- vector(mode="list", length(levels(Tissue)))
    names(gg) <- levels(Tissue)
    
    for (t in levels(Tissue)) {
      if (log2==FALSE) {
        gg[[t]] <-
          ggplot(summary.stats[[t]], aes(x=Mean, y=SD)) + 
          geom_point(color="#440154", size=0.2, alpha=0.6) +
          stat_density_2d(aes(fill=..level..), geom="polygon") +
          scale_fill_continuous(type = "viridis") +
          stat_smooth(method="loess",
                      span=0.60,
                      formula=y~x,
                      se=FALSE) +
          stat_cor(method="pearson") +
          labs(title=gsub("_"," ",t)) +
          theme_bw()
      } else {
        gg[[t]] <-
          ggplot(summary.stats[[t]], aes(x=Mean, y=log2(SD))) + 
          geom_point(color="#440154", size=0.2, alpha=0.6) +
          stat_density_2d(aes(fill=..level..), geom="polygon") +
          scale_fill_continuous(type = "viridis") +
          stat_smooth(method="loess",
                      span=0.60,
                      formula=y~x,
                      se=FALSE) +
          stat_cor(method="pearson") +
          labs(title=gsub("_"," ",t)) +
          theme_bw()          
      }
    }
  }
  
  for (i in seq(1,length(gg),4)) {
    if (i+3 > length(gg)) { print(plot_grid(plotlist=gg[i:length(gg)], nrow=2, ncol=2)) } 
    else { print(plot_grid(plotlist=gg[i:(i+3)])) }
  }  
}

plot_cv_vs_mean <- function(summary.stats, metadata, log2=FALSE) {
  Tissue <- as.factor(metadata$Tissue)
  
  if ("Sex" %in% colnames(metadata)) {
    Sex <- as.factor(metadata$Sex)
    gg <- vector(mode="list", length(levels(Tissue))*length(levels(Sex)))
    names(gg) <- levels(interaction(Tissue, Sex, sep="_", lex.order=TRUE))
    
    for (t in levels(Tissue)) {
      for (s in levels(Sex)) {
        if (nrow(summary.stats[[t]][[s]]>0)) {
          if (log2==FALSE) {
            gg[[paste(t,s,sep="_")]] <-
              ggplot(summary.stats[[t]][[s]], aes(x=Mean, y=CV2)) + 
              geom_point(color="#440154", size=0.2, alpha=0.6) +
              stat_density_2d(aes(fill=..level..), geom="polygon") +
              scale_fill_continuous(type = "viridis") +
              stat_smooth(method="loess",
                          span=0.60,
                          formula=y~x,
                          se=FALSE) +
              stat_cor(method="pearson") +
              labs(title=gsub("_"," ",t), 
                   subtitle=gsub("_"," ", s)) +
              theme_bw()
          } else {
            gg[[paste(t,s,sep="_")]] <-
              ggplot(summary.stats[[t]][[s]], aes(x=Mean, y=log2(CV2))) + 
              geom_point(color="#440154", size=0.2, alpha=0.6) +
              stat_density_2d(aes(fill=..level..), geom="polygon") +
              scale_fill_continuous(type = "viridis") +
              stat_smooth(method="loess",
                          span=0.60,
                          formula=y~x,
                          se=FALSE) +
              stat_cor(method="pearson") +
              labs(title=gsub("_"," ",t), 
                   subtitle=gsub("_"," ", s)) +
              theme_bw()          
          }
        } else {
          gg[[paste(t,s,sep="_")]] <- NULL
        }
      }
    }
  } else {
    gg <- vector(mode="list", length(levels(Tissue)))
    names(gg) <- levels(Tissue)
    
    for (t in levels(Tissue)) {
      if (log2==FALSE) {
        gg[[t]] <-
          ggplot(summary.stats[[t]], aes(x=Mean, y=CV2)) + 
          geom_point(color="#440154", size=0.2, alpha=0.6) +
          stat_density_2d(aes(fill=..level..), geom="polygon") +
          scale_fill_continuous(type = "viridis") +
          stat_smooth(method="loess",
                      span=0.60,
                      formula=y~x,
                      se=FALSE) +
          stat_cor(method="pearson") +
          labs(title=gsub("_"," ",t)) +
          theme_bw()
      } else {
        gg[[t]] <-
          ggplot(summary.stats[[t]], aes(x=Mean, y=log2(CV2))) + 
          geom_point(color="#440154", size=0.2, alpha=0.6) +
          stat_density_2d(aes(fill=..level..), geom="polygon") +
          scale_fill_continuous(type = "viridis") +
          stat_smooth(method="loess",
                      span=0.60,
                      formula=y~x,
                      se=FALSE) +
          stat_cor(method="pearson") +
          labs(title=gsub("_"," ",t)) +
          theme_bw()          
      }
    }
  }
  
  for (i in seq(1,length(gg),4)) {
    if (i+3 > length(gg)) { print(plot_grid(plotlist=gg[i:length(gg)], nrow=2, ncol=2)) } 
    else { print(plot_grid(plotlist=gg[i:(i+3)])) }
  }  
}

plot_jackknife_adj_sd_vs_mean <- function(jack.summary, metadata) {
  Tissue <- as.factor(metadata$Tissue)
  xlab <- "Mean of jackknifed means"
  ylab <- "Mean of jackknifed adjusted SD"
  
  if ("Sex" %in% colnames(metadata)) {
    Sex <- as.factor(metadata$Sex)
    gg <- vector(mode="list", length(levels(Tissue))*length(levels(Sex)))
    names(gg) <- levels(interaction(Tissue, Sex, sep="_", lex.order=TRUE))
    
    for (t in levels(Tissue)) {
      for (s in levels(Sex)) {
        if (!is.null(jack.summary[[t]][[s]])) {
          gg[[paste(t,s,sep="_")]] <-
            ggplot(jack.summary[[t]][[s]], aes(x=Mean_Mean, y=Mean_AdjSD)) + 
            geom_point(color="#440154", size=0.2, alpha=0.6) +
            stat_density_2d(aes(fill=..level..), geom="polygon") +
            scale_fill_continuous(type = "viridis") +
            stat_smooth(method="loess",
                        span=0.60,
                        formula=y~x,
                        se=FALSE) +
            stat_cor(method="pearson") +
            labs(title=gsub("_"," ",t), 
                 subtitle=gsub("_"," ", s),
                 x=xlab, y=ylab) +
            theme_bw()
        } else {
          gg[[paste(t,s,sep="_")]] <- NULL
        }
      }
    }
  } else {
    gg <- vector(mode="list", length(levels(Tissue)))
    names(gg) <- levels(Tissue)
    
    for (t in levels(Tissue)) {
      gg[[t]] <-
        ggplot(jack.summary[[t]], aes(x=Mean_Mean, y=Mean_AdjSD)) + 
        geom_point(color="#440154", size=0.2, alpha=0.6) +
        stat_density_2d(aes(fill=..level..), geom="polygon") +
        scale_fill_continuous(type = "viridis") +
        stat_smooth(method="loess",
                    span=0.60,
                    formula=y~x,
                    se=FALSE) +
        stat_cor(method="pearson") +
        labs(title=gsub("_"," ",t),
             x=xlab, y=ylab) +
        theme_bw()
    }         
  }
  
  for (i in seq(1,length(gg),4)) {
    if (i+3 > length(gg)) { print(plot_grid(plotlist=gg[i:length(gg)], nrow=2, ncol=2)) } 
    else { print(plot_grid(plotlist=gg[i:(i+3)])) }
  }  
}

plot_jackknife_log2sd_vs_mean <- function(jack.summary, metadata) {
  Tissue <- as.factor(metadata$Tissue)
  xlab <- "Mean of jackknifed means"
  ylab <- "Mean of jackknifed log2(SD)"
  
  if ("Sex" %in% colnames(metadata)) {
    Sex <- as.factor(metadata$Sex)
    gg <- vector(mode="list", length(levels(Tissue))*length(levels(Sex)))
    names(gg) <- levels(interaction(Tissue, Sex, sep="_", lex.order=TRUE))
    
    for (t in levels(Tissue)) {
      for (s in levels(Sex)) {
        if (!is.null(jack.summary[[t]][[s]])) {
          gg[[paste(t,s,sep="_")]] <-
            ggplot(jack.summary[[t]][[s]], aes(x=Mean_Mean, y=Mean_Log2SD)) + 
            geom_point(color="#440154", size=0.2, alpha=0.6) +
            stat_density_2d(aes(fill=..level..), geom="polygon") +
            scale_fill_continuous(type = "viridis") +
            stat_smooth(method="loess",
                        span=0.60,
                        formula=y~x,
                        se=FALSE) +
            stat_cor(method="pearson") +
            labs(title=gsub("_"," ",t), 
                 subtitle=gsub("_"," ", s),
                 x=xlab, y=ylab) +
            theme_bw()
        } else {
          gg[[paste(t,s,sep="_")]] <- NULL
        }
      }
    }
  } else {
    gg <- vector(mode="list", length(levels(Tissue)))
    names(gg) <- levels(Tissue)
    
    for (t in levels(Tissue)) {
      gg[[t]] <-
        ggplot(jack.summary[[t]], aes(x=Mean_Mean, y=Mean_Log2SD)) + 
        geom_point(color="#440154", size=0.2, alpha=0.6) +
        stat_density_2d(aes(fill=..level..), geom="polygon") +
        scale_fill_continuous(type = "viridis") +
        stat_smooth(method="loess",
                    span=0.60,
                    formula=y~x,
                    se=FALSE) +
        stat_cor(method="pearson") +
        labs(title=gsub("_"," ",t),
             x=xlab, y=ylab) +
        theme_bw()
    }         
  }
  
  for (i in seq(1,length(gg),4)) {
    if (i+3 > length(gg)) { print(plot_grid(plotlist=gg[i:length(gg)], nrow=2, ncol=2)) } 
    else { print(plot_grid(plotlist=gg[i:(i+3)])) }
  }  
}

plot_jackknife_log2cv_vs_mean <- function(jack.summary, metadata) {
  Tissue <- as.factor(metadata$Tissue)
  xlab <- "Mean of jackknifed means"
  ylab <- "Mean of jackknifed log2(CV^2)"
  
  if ("Sex" %in% colnames(metadata)) {
    Sex <- as.factor(metadata$Sex)
    gg <- vector(mode="list", length(levels(Tissue))*length(levels(Sex)))
    names(gg) <- levels(interaction(Tissue, Sex, sep="_", lex.order=TRUE))
    
    for (t in levels(Tissue)) {
      for (s in levels(Sex)) {
        if (!is.null(jack.summary[[t]][[s]])) {
          gg[[paste(t,s,sep="_")]] <-
            ggplot(jack.summary[[t]][[s]], aes(x=Mean_Mean, y=Mean_Log2CV)) + 
            geom_point(color="#440154", size=0.2, alpha=0.6) +
            stat_density_2d(aes(fill=..level..), geom="polygon") +
            scale_fill_continuous(type = "viridis") +
            stat_smooth(method="loess",
                        span=0.60,
                        formula=y~x,
                        se=FALSE) +
            stat_cor(method="pearson") +
            labs(title=gsub("_"," ",t), 
                 subtitle=gsub("_"," ", s),
                 x=xlab, y=ylab) +
            theme_bw()
        } else {
          gg[[paste(t,s,sep="_")]] <- NULL
        }
      }
    }
  } else {
    gg <- vector(mode="list", length(levels(Tissue)))
    names(gg) <- levels(Tissue)
    
    for (t in levels(Tissue)) {
      gg[[t]] <-
        ggplot(jack.summary[[t]], aes(x=Mean_Mean, y=Mean_Log2CV)) + 
        geom_point(color="#440154", size=0.2, alpha=0.6) +
        stat_density_2d(aes(fill=..level..), geom="polygon") +
        scale_fill_continuous(type = "viridis") +
        stat_smooth(method="loess",
                    span=0.60,
                    formula=y~x,
                    se=FALSE) +
        stat_cor(method="pearson") +
        labs(title=gsub("_"," ",t),
             x=xlab, y=ylab) +
        theme_bw()
    }         
  }
  
  for (i in seq(1,length(gg),4)) {
    if (i+3 > length(gg)) { print(plot_grid(plotlist=gg[i:length(gg)], nrow=2, ncol=2)) } 
    else { print(plot_grid(plotlist=gg[i:(i+3)])) }
  }  
}

plot_jackknife_resid_var_vs_mean <- function(jack.summary, metadata, method=c("sd","cv")) {
  Tissue <- as.factor(metadata$Tissue)
  xlab <- "Mean of jackknifed means"
  if (method=="sd") {
    ylab <- "Mean of residual log2(SD)"    
  } else if (method=="cv") {
    ylab <- "Mean of residual log2(CV^2)"
  }

  
  if ("Sex" %in% colnames(metadata)) {
    Sex <- as.factor(metadata$Sex)
    gg <- vector(mode="list", length(levels(Tissue))*length(levels(Sex)))
    names(gg) <- levels(interaction(Tissue, Sex, sep="_", lex.order=TRUE))
    
    for (t in levels(Tissue)) {
      for (s in levels(Sex)) {
        if (!is.null(jack.summary[[t]][[s]])) {
          if (method=="sd") {
            gg[[paste(t,s,sep="_")]] <-
              ggplot(jack.summary[[t]][[s]], aes(x=Mean_Mean, y=Mean_Resid_Log2SD)) + 
              geom_point(color="#440154", size=0.2, alpha=0.6) +
              stat_density_2d(aes(fill=..level..), geom="polygon") +
              scale_fill_continuous(type = "viridis") +
              stat_smooth(method="loess",
                          span=0.60,
                          formula=y~x,
                          se=FALSE) +
              stat_cor(method="pearson") +
              labs(title=gsub("_"," ",t), 
                   subtitle=gsub("_"," ", s),
                   x=xlab, y=ylab) +
              theme_bw()
          } else if (method=="cv") {
            gg[[paste(t,s,sep="_")]] <-
              ggplot(jack.summary[[t]][[s]], aes(x=Mean_Mean, y=Mean_Resid_Log2CV)) + 
              geom_point(color="#440154", size=0.2, alpha=0.6) +
              stat_density_2d(aes(fill=..level..), geom="polygon") +
              scale_fill_continuous(type = "viridis") +
              stat_smooth(method="loess",
                          span=0.60,
                          formula=y~x,
                          se=FALSE) +
              stat_cor(method="pearson") +
              labs(title=gsub("_"," ",t), 
                   subtitle=gsub("_"," ", s),
                   x=xlab, y=ylab) +
              theme_bw()            
          }

        } else {
          gg[[paste(t,s,sep="_")]] <- NULL
        }
      }
    }
  } else {
    gg <- vector(mode="list", length(levels(Tissue)))
    names(gg) <- levels(Tissue)
    
    for (t in levels(Tissue)) {
      if (method=="sd") {
        gg[[t]] <-
          ggplot(jack.summary[[t]], aes(x=Mean_Mean, y=Mean_Resid_Log2SD)) + 
          geom_point(color="#440154", size=0.2, alpha=0.6) +
          stat_density_2d(aes(fill=..level..), geom="polygon") +
          scale_fill_continuous(type = "viridis") +
          stat_smooth(method="loess",
                      span=0.60,
                      formula=y~x,
                      se=FALSE) +
          stat_cor(method="pearson") +
          labs(title=gsub("_"," ",t),
               x=xlab, y=ylab) +
          theme_bw()        
      } else if (method=="cv") {
        gg[[t]] <-
          ggplot(jack.summary[[t]], aes(x=Mean_Mean, y=Mean_Resid_Log2CV)) + 
          geom_point(color="#440154", size=0.2, alpha=0.6) +
          stat_density_2d(aes(fill=..level..), geom="polygon") +
          scale_fill_continuous(type = "viridis") +
          stat_smooth(method="loess",
                      span=0.60,
                      formula=y~x,
                      se=FALSE) +
          stat_cor(method="pearson") +
          labs(title=gsub("_"," ",t),
               x=xlab, y=ylab) +
          theme_bw()        
        }
      }         
    }
  
  for (i in seq(1,length(gg),4)) {
    if (i+3 > length(gg)) { print(plot_grid(plotlist=gg[i:length(gg)], nrow=2, ncol=2)) } 
    else { print(plot_grid(plotlist=gg[i:(i+3)])) }
  }  
}

plot_jackknife_local_var_rank_vs_mean <- function(jack.summary, metadata, method=c("sd","cv")) {
  Tissue <- as.factor(metadata$Tissue)
  xlab <- "Mean of jackknifed means"
  ylab <- "Mean of local variability rank"
  
  if ("Sex" %in% colnames(metadata)) {
    Sex <- as.factor(metadata$Sex)
    gg <- vector(mode="list", length(levels(Tissue))*length(levels(Sex)))
    names(gg) <- levels(interaction(Tissue, Sex, sep="_", lex.order=TRUE))
    
    for (t in levels(Tissue)) {
      for (s in levels(Sex)) {
        if (!is.null(jack.summary[[t]][[s]])) {
          if (method=="sd") {
            gg[[paste(t,s,sep="_")]] <-
              ggplot(jack.summary[[t]][[s]], aes(x=Mean_Mean, y=Mean_Local_Rank_Log2SD)) + 
              geom_point(color="#440154", size=0.2, alpha=0.6) +
              stat_density_2d(aes(fill=..level..), geom="polygon") +
              scale_fill_continuous(type = "viridis") +
              stat_smooth(method="loess",
                          span=0.60,
                          formula=y~x,
                          se=FALSE) +
              stat_cor(method="pearson") +
              labs(title=gsub("_"," ",t), 
                   subtitle=gsub("_"," ", s),
                   x=xlab, y=ylab) +
              theme_bw()            
          } else if (method=="cv") {
            gg[[paste(t,s,sep="_")]] <-
              ggplot(jack.summary[[t]][[s]], aes(x=Mean_Mean, y=Mean_Local_Rank_Log2CV)) + 
              geom_point(color="#440154", size=0.2, alpha=0.6) +
              stat_density_2d(aes(fill=..level..), geom="polygon") +
              scale_fill_continuous(type = "viridis") +
              stat_smooth(method="loess",
                          span=0.60,
                          formula=y~x,
                          se=FALSE) +
              stat_cor(method="pearson") +
              labs(title=gsub("_"," ",t), 
                   subtitle=gsub("_"," ", s),
                   x=xlab, y=ylab) +
              theme_bw()            
          }
        } else {
          gg[[paste(t,s,sep="_")]] <- NULL
        }
      }
    }
  } else {
    gg <- vector(mode="list", length(levels(Tissue)))
    names(gg) <- levels(Tissue)
    
    for (t in levels(Tissue)) {
      if (method=="sd") {
        gg[[t]] <-
          ggplot(jack.summary[[t]], aes(x=Mean_Mean, y=Mean_Local_Rank_Log2SD)) + 
          geom_point(color="#440154", size=0.2, alpha=0.6) +
          stat_density_2d(aes(fill=..level..), geom="polygon") +
          scale_fill_continuous(type = "viridis") +
          stat_smooth(method="loess",
                      span=0.60,
                      formula=y~x,
                      se=FALSE) +
          stat_cor(method="pearson") +
          labs(title=gsub("_"," ",t),
               x=xlab, y=ylab) +
          theme_bw()        
      } else if (method=="cv") {
        gg[[t]] <-
          ggplot(jack.summary[[t]], aes(x=Mean_Mean, y=Mean_Local_Rank_Log2CV)) + 
          geom_point(color="#440154", size=0.2, alpha=0.6) +
          stat_density_2d(aes(fill=..level..), geom="polygon") +
          scale_fill_continuous(type = "viridis") +
          stat_smooth(method="loess",
                      span=0.60,
                      formula=y~x,
                      se=FALSE) +
          stat_cor(method="pearson") +
          labs(title=gsub("_"," ",t),
               x=xlab, y=ylab) +
          theme_bw()        
      }
    }         
  }
  
  for (i in seq(1,length(gg),4)) {
    if (i+3 > length(gg)) { print(plot_grid(plotlist=gg[i:length(gg)], nrow=2, ncol=2)) } 
    else { print(plot_grid(plotlist=gg[i:(i+3)])) }
  }  
}

plot_jackknife_global_var_rank_vs_mean <- function(jack.summary, metadata, method=c("sd","cv")) {
  Tissue <- as.factor(metadata$Tissue)
  xlab <- "Mean of jackknifed means"
  ylab <- "Global variability rank"
  
  if ("Sex" %in% colnames(metadata)) {
    Sex <- as.factor(metadata$Sex)
    gg <- vector(mode="list", length(levels(Tissue))*length(levels(Sex)))
    names(gg) <- levels(interaction(Tissue, Sex, sep="_", lex.order=TRUE))
    
    for (t in levels(Tissue)) {
      for (s in levels(Sex)) {
        if (!is.null(jack.summary[[t]][[s]])) {
          if (method=="sd") {
            gg[[paste(t,s,sep="_")]] <-
              ggplot(jack.summary[[t]][[s]], aes(x=Mean_Mean, y=Global_Rank_Log2SD)) + 
              geom_point(color="#440154", size=0.2, alpha=0.6) +
              stat_density_2d(aes(fill=..level..), geom="polygon") +
              scale_fill_continuous(type = "viridis") +
              stat_smooth(method="loess",
                          span=0.60,
                          formula=y~x,
                          se=FALSE) +
              stat_cor(method="pearson") +
              labs(title=gsub("_"," ",t), 
                   subtitle=gsub("_"," ", s),
                   x=xlab, y=ylab) +
              theme_bw()            
          } else if (method=="cv") {
            gg[[paste(t,s,sep="_")]] <-
              ggplot(jack.summary[[t]][[s]], aes(x=Mean_Mean, y=Global_Rank_Log2CV)) + 
              geom_point(color="#440154", size=0.2, alpha=0.6) +
              stat_density_2d(aes(fill=..level..), geom="polygon") +
              scale_fill_continuous(type = "viridis") +
              stat_smooth(method="loess",
                          span=0.60,
                          formula=y~x,
                          se=FALSE) +
              stat_cor(method="pearson") +
              labs(title=gsub("_"," ",t), 
                   subtitle=gsub("_"," ", s),
                   x=xlab, y=ylab) +
              theme_bw()            
          }
        } else {
          gg[[paste(t,s,sep="_")]] <- NULL
        }
      }
    }
  } else {
    gg <- vector(mode="list", length(levels(Tissue)))
    names(gg) <- levels(Tissue)
    
    for (t in levels(Tissue)) {
      if (method=="sd") {
        gg[[t]] <-
          ggplot(jack.summary[[t]], aes(x=Mean_Mean, y=Global_Rank_Log2SD)) + 
          geom_point(color="#440154", size=0.2, alpha=0.6) +
          stat_density_2d(aes(fill=..level..), geom="polygon") +
          scale_fill_continuous(type = "viridis") +
          stat_smooth(method="loess",
                      span=0.60,
                      formula=y~x,
                      se=FALSE) +
          stat_cor(method="pearson") +
          labs(title=gsub("_"," ",t),
               x=xlab, y=ylab) +
          theme_bw()        
      } else if (method=="cv") {
        gg[[t]] <-
          ggplot(jack.summary[[t]], aes(x=Mean_Mean, y=Global_Rank_Log2CV)) + 
          geom_point(color="#440154", size=0.2, alpha=0.6) +
          stat_density_2d(aes(fill=..level..), geom="polygon") +
          scale_fill_continuous(type = "viridis") +
          stat_smooth(method="loess",
                      span=0.60,
                      formula=y~x,
                      se=FALSE) +
          stat_cor(method="pearson") +
          labs(title=gsub("_"," ",t),
               x=xlab, y=ylab) +
          theme_bw()        
      }
    }         
  }
  
  for (i in seq(1,length(gg),4)) {
    if (i+3 > length(gg)) { print(plot_grid(plotlist=gg[i:length(gg)], nrow=2, ncol=2)) } 
    else { print(plot_grid(plotlist=gg[i:(i+3)])) }
  }  
}

# Annotate boxplot with number of observations per grouping
give_n <- function(x) {
  return(c(y = median(x)*1.10, label = length(x))) 
}

boxplot_var_rank_lncrna_vs_coding <- function(jack.bind, sex=TRUE, method=c("sd","cv")) {
  if (sex==TRUE) { 
    if (method=="sd") {
      ggplot(jack.bind,
             aes(x=Sex, y=Mean_Local_Rank_Log2SD, fill=Biotype)) +
        geom_boxplot(notch=TRUE) +
        stat_summary(fun.data = give_n, geom = "text", fun = median, 
                     position = position_dodge(width = 0.75)) +
        stat_compare_means(aes(label=..p.signif..), method="wilcox.test") +
        geom_hline(yintercept=0.5, linetype=9) +
        labs(caption=paste0("produced on ", Sys.time()),
             y="Variability rank") +
        facet_wrap(vars(Tissue)) +
        theme_bw()        
    } else if (method=="cv") {
      ggplot(jack.bind,
             aes(x=Sex, y=Mean_Local_Rank_Log2CV, fill=Biotype)) +
        geom_boxplot(notch=TRUE) +
        stat_summary(fun.data = give_n, geom = "text", fun = median, 
                     position = position_dodge(width = 0.75)) +
        stat_compare_means(aes(label=..p.signif..), method="wilcox.test") +
        geom_hline(yintercept=0.5, linetype=9) +
        labs(caption=paste0("produced on ", Sys.time()),
             y="Variability rank") +
        facet_wrap(vars(Tissue)) +
        theme_bw()        
    }

  } else {
    if (method=="sd") {
      ggplot(jack.bind,
             aes(x=Biotype, y=Mean_Local_Rank_Log2SD, fill=Biotype)) +
        geom_boxplot(notch=TRUE) +
        stat_summary(fun.data = give_n, geom = "text", fun = median, 
                     position = position_dodge(width = 0.75)) +
        stat_compare_means(aes(label=..p.signif..), method="wilcox.test") +
        geom_hline(yintercept=0.5, linetype=9) +
        labs(caption=paste0("produced on ", Sys.time()),
             y="Variability rank") +
        facet_wrap(vars(Tissue)) +
        theme_bw()        
    } else if (method=="cv") {
      ggplot(jack.bind,
             aes(x=Biotype, y=Mean_Local_Rank_Log2CV, fill=Biotype)) +
        geom_boxplot(notch=TRUE) +
        stat_summary(fun.data = give_n, geom = "text", fun = median, 
                     position = position_dodge(width = 0.75)) +
        stat_compare_means(aes(label=..p.signif..), method="wilcox.test") +
        geom_hline(yintercept=0.5, linetype=9) +
        labs(caption=paste0("produced on ", Sys.time()),
             y="Variability rank") +
        facet_wrap(vars(Tissue)) +
        theme_bw()           
    }
  }
}