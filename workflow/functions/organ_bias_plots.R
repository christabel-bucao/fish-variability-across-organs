#### Notes ####

#### Libraries ####
library(ggplot2)
library(viridis)
library(ggpubr)

#### Sourced functions ####
source("../functions/annotate_boxplots.R")

#### Functions ####

set_facet_vars <- function(df) {
  if ("Sex" %in% colnames(df)) {
    facet.vars <- vars(Tissue, Sex)
  } else { facet.vars <- vars(Tissue) }
  
  return(facet.vars)
}

density_plot_variability_vs_mean <- function(df, title="") {
  facet.vars <- set_facet_vars(df)
  df$Tissue <- gsub("_", " ", df$Tissue)
  
  ggplot(df,
         aes(x=Mean_Mean, y=Mean_Local_Rank_Log2CV)) +
    geom_point(color="#440154", size=0.2, alpha=0.6) +
    stat_density_2d(aes(fill=..level..), geom="polygon") +
    scale_fill_continuous(type="viridis") +
    stat_smooth(method="loess",
                span=0.60,
                formula=y~x,
                se=FALSE) +
    stat_cor(method="pearson", size=5, label.y=1.1) +
    scale_y_continuous(breaks=seq(0.0,1.15,0.20)) +
    labs(title=title, 
         x="Jackknifed mean expression", 
         y="Expression variability rank") +
    facet_wrap(facet.vars, nrow=2) +
    theme_bw() +
    theme(plot.title=element_text(size=20, face="italic"),
          axis.title=element_text(size=18),
          axis.text=element_text(size=14),
          strip.text=element_text(size=14),
          panel.spacing.x=unit(1.1,"lines"),
          legend.position="none") 
}

density_plot_mean_vs_tau <- function(df, title="") {
  facet.vars <- set_facet_vars(df)
  df$Tissue <- gsub("_", " ", df$Tissue)
  
  ggplot(df,
         aes(x=tau, y=Mean_Mean)) +
    geom_point(color="#440154", size=0.2, alpha=0.6) +
    stat_density_2d(aes(fill=..level..), geom="polygon") +
    scale_fill_continuous(type="viridis") +
    stat_smooth(method="loess",
                span=0.60,
                formula=y~x,
                se=FALSE) +
    stat_cor(method="pearson", size=5) +
    labs(title=title, y="Jackknifed mean expression") +
    facet_wrap(facet.vars, nrow=2) +
    theme_bw() +
    theme(plot.title=element_text(size=20, face="italic"),
          axis.title=element_text(size=18),
          axis.text=element_text(size=10),
          strip.text=element_text(size=14),
          panel.spacing.x=unit(1.1,"lines"),
          legend.position="none")   
}

density_plot_resid_var_vs_tau <- function(df, title="") {
  facet.vars <- set_facet_vars(df)
  df$Tissue <- gsub("_", " ", df$Tissue)
  
  ggplot(df,
         aes(x=tau, y=Mean_Resid_Log2CV)) +
    geom_point(color="#440154", size=0.2, alpha=0.6) +
    stat_density_2d(aes(fill=..level..), geom="polygon") +
    scale_fill_continuous(type = "viridis") +
    stat_smooth(method="loess",
                span=0.60,
                formula=y~x,
                se=FALSE) +
    stat_cor(method="pearson", size=5) +
    labs(title=title, y="Residual expression variation") +
    facet_wrap(facet.vars, nrow=2) +
    theme_bw() +
    theme(plot.title=element_text(size=20, face="italic"),
          axis.title=element_text(size=18),
          axis.text=element_text(size=10),
          strip.text=element_text(size=14),
          panel.spacing.x=unit(1.1,"lines"),
          legend.position="none")
}

density_plot_variability_vs_tau <- function(df, title="") {
  facet.vars <- set_facet_vars(df)
  df$Tissue <- gsub("_", " ", df$Tissue)
  
  ggplot(df,
         aes(x=tau, y=Mean_Local_Rank_Log2CV)) +
    geom_point(color="#440154", size=0.2, alpha=0.6) +
    stat_density_2d(aes(fill=..level..), geom="polygon") +
    scale_fill_continuous(type = "viridis") +
    stat_smooth(method="loess",
                span=0.60,
                formula=y~x,
                se=FALSE) +
    stat_cor(method="pearson", size=5) +
    labs(title=title, y="Expression variability rank") +
    facet_wrap(facet.vars, nrow=2) +
    theme_bw() +
    theme(plot.title=element_text(size=20, face="italic"),
          axis.title=element_text(size=18),
          axis.text=element_text(size=10),
          strip.text=element_text(size=14),
          panel.spacing.x=unit(1.1,"lines"),
          legend.position="none")
}

boxplot_mean_by_tau_quantiles <- function(df, title="") {
  facet.vars <- set_facet_vars(df)
  df$Tissue <- gsub("_", " ", df$Tissue)
  
  ggplot(df,
         aes(x=Quantile, y=Mean_Mean, fill=Quantile)) +
    stat_boxplot(geom ='errorbar') +
    geom_boxplot(notch=TRUE) +
    stat_summary(fun.data=annotate_n, geom="text", fun=median, 
                 position = position_dodge(width=0.75), size=4.5) +
    labs(title=params$species.name,
         x="Organ specificity (tau) quantile",
         y="Jackknifed mean expression") +
    facet_wrap(facet.vars, nrow=2) +
    theme_bw() +
    theme(plot.title=element_text(size=20, face="italic"),
          axis.title=element_text(size=18),
          axis.text=element_text(size=14),
          strip.text=element_text(size=14),
          panel.spacing.x=unit(1.1,"lines"),
          legend.position="none")
}

boxplot_resid_var_by_tau_quantiles <- function(df, title="") {
  facet.vars <- set_facet_vars(df)
  df$Tissue <- gsub("_", " ", df$Tissue)
  
  ggplot(df,
         aes(x=Quantile, y=Mean_Resid_Log2CV, fill=Quantile)) +
    stat_boxplot(geom ='errorbar') +
    geom_boxplot(notch=TRUE) +
    stat_summary(fun.data=annotate_n, geom="text", fun=median, 
                 position = position_dodge(width=0.75), size=4.5) +
    labs(title=params$species.name,
         x="Organ specificity (tau) quantile",
         y="Residual expression variation") +
    facet_wrap(facet.vars, nrow=2) +
    theme_bw() +
    theme(plot.title=element_text(size=20, face="italic"),
          axis.title=element_text(size=18),
          axis.text=element_text(size=14),
          strip.text=element_text(size=14),
          panel.spacing.x=unit(1.1,"lines"),
          legend.position="none")  
}

boxplot_variability_by_tau_quantiles <- function(df, title="") {
  facet.vars <- set_facet_vars(df)
  df$Tissue <- gsub("_", " ", df$Tissue)
  
  ggplot(df,
         aes(x=Quantile, y=Mean_Local_Rank_Log2CV, fill=Quantile)) +
    stat_boxplot(geom ='errorbar') +
    geom_boxplot(notch=TRUE) +
    stat_summary(fun.data=annotate_n, geom="text", fun=median, 
                 position = position_dodge(width=0.75), size=4.5) +
    labs(title=params$species.name,
         x="Organ specificity (tau) quantile",
         y="Expression variability rank") +
    facet_wrap(facet.vars, nrow=2) +
    theme_bw() +
    theme(plot.title=element_text(size=20, face="italic"),
          axis.title=element_text(size=18),
          axis.text=element_text(size=14),
          strip.text=element_text(size=14),
          panel.spacing.x=unit(1.1,"lines"),
          legend.position="none")    
}

boxplot_mean_by_organ_bias <- function(df, species="", tissue="", sex=NA) {
  tissue <- gsub("_", " ", tissue)
  comparisons <- list(c("broad","focal"),c(tissue,"other"))
  
  df$Tissue <- gsub("_", " ", df$Tissue)
  df$TissueBias <- gsub("_", " ", df$TissueBias)
  
  # Reorder factor levels
  df$TissueBias <- factor(df$TissueBias, levels=c("broad","focal",tissue,"other"))   
  
  ggplot(df,
         aes(x=TissueBias, y=Mean_Mean, fill=TissueBias)) +
    stat_boxplot(geom='errorbar') +
    geom_boxplot(notch=TRUE) +
    stat_summary(fun.data=annotate_n, geom="text", fun=median,
                 position=position_dodge(width=0.75), size=4.5) +
    geom_hline(yintercept=median(df$Mean_Mean),
               linetype=3, size=0.5) +
    stat_compare_means(aes(label=..p.signif..),
                       test="wilcox.test",
                       comparisons=comparisons) +
    facet_wrap(vars(Tissue), nrow=2) +
    labs(title=paste(species, " (", sex, ")", sep=""),
         subtitle=paste(tissue, "biased genes", sep="-"),
         x="Organ bias",
         y="Jackknifed mean expression",
         fill="Organ bias") +
    theme_bw() +
    theme(plot.title=element_text(size=20),
          plot.subtitle=element_text(size=18),
          axis.title=element_text(size=18),
          axis.text=element_text(size=14),
          strip.text=element_text(size=14),
          panel.spacing.x=unit(1.1,"lines"),
          legend.position="none")
}

boxplot_resid_var_by_organ_bias <- function(df, species="", tissue="", sex=NA) {
  tissue <- gsub("_", " ", tissue)
  comparisons <- list(c("broad","focal"),c(tissue,"other"))
  
  df$Tissue <- gsub("_", " ", df$Tissue)
  df$TissueBias <- gsub("_", " ", df$TissueBias)
  
  # Reorder factor levels
  df$TissueBias <- factor(df$TissueBias, levels=c("broad","focal",tissue,"other"))  
  
  ggplot(df,
         aes(x=TissueBias, y=Mean_Resid_Log2CV, fill=TissueBias)) +
    stat_boxplot(geom='errorbar') +
    geom_boxplot(notch=TRUE) +
    stat_summary(fun.data=annotate_n, geom="text", fun=median,
                 position=position_dodge(width=0.75), size=4.5) +
    #geom_hline(yintercept=median(df$Mean_Local_Rank_Log2CV), linetype=3, size=0.5) +
    geom_hline(yintercept=0.0, linetype=3, size=0.5) +
    stat_compare_means(aes(label=..p.signif..),
                       test="wilcox.test",
                       comparisons=comparisons) +
    facet_wrap(vars(Tissue), nrow=2) +
    labs(title=paste(species, " (", sex, ")", sep=""),
         subtitle=paste(tissue, "biased genes", sep="-"),
         x="Organ bias",
         y="Residual expression variation",
         fill="Organ bias") +
    theme_bw() +
    theme(plot.title=element_text(size=20),
          plot.subtitle=element_text(size=18),
          axis.title=element_text(size=18),
          axis.text=element_text(size=14),
          strip.text=element_text(size=14),
          panel.spacing.x=unit(1.1,"lines"),
          legend.position="none")  
}

boxplot_variability_by_organ_bias <- function(df, species="", tissue="", sex=NA) {
  tissue <- gsub("_", " ", tissue)
  comparisons <- list(c("broad","focal"),c(tissue,"other"))
  
  df$Tissue <- gsub("_", " ", df$Tissue)
  df$TissueBias <- gsub("_", " ", df$TissueBias)
  
  # Reorder factor levels
  df$TissueBias <- factor(df$TissueBias, levels=c("broad","focal",tissue,"other"))  
  
  ggplot(df,
         aes(x=TissueBias, y=Mean_Local_Rank_Log2CV, fill=TissueBias)) +
    stat_boxplot(geom='errorbar') +
    geom_boxplot(notch=TRUE) +
    stat_summary(fun.data=annotate_n, geom="text", fun=median,
                 position=position_dodge(width=0.75), size=4.5) +
    #geom_hline(yintercept=median(df$Mean_Local_Rank_Log2CV), linetype=3, size=0.5) +
    geom_hline(yintercept=0.5, linetype=3, size=0.5) +
    stat_compare_means(aes(label=..p.signif..),
                       test="wilcox.test",
                       comparisons=comparisons) +
    facet_wrap(vars(Tissue), nrow=2) +
    labs(title=paste(species, " (", sex, ")", sep=""),
         subtitle=paste(tissue, "biased genes", sep="-"),
         x="Organ bias",
         y="Expression variability rank",
         fill="Organ bias") +
    theme_bw() +
    theme(plot.title=element_text(size=20),
          plot.subtitle=element_text(size=18),
          axis.title=element_text(size=18),
          axis.text=element_text(size=14),
          strip.text=element_text(size=14),
          panel.spacing.x=unit(1.1,"lines"),
          legend.position="none")
}

boxplot_variability_example <- function(df, species="", tissue="", ylab="Expression variability rank", color) {
  colors <- c("#999999","#999999",color,color,color)
  
  facet.vars <- set_facet_vars(df)
  tissue <- gsub("_", " ", tissue)
  comparisons <- list(c("broad","focal"),c(tissue,"other"))
  
  df$Tissue <- gsub("_", " ", df$Tissue)
  
  # Reorder factor levels
  df$TissueBias <- factor(df$TissueBias, levels=c("broad","focal",tissue,"other"))  
  
  ggplot(df,
         aes(x=TissueBias, y=Mean_Local_Rank_Log2CV, fill=TissueBias)) +
    stat_boxplot(geom='errorbar') +
    geom_boxplot(notch=TRUE, alpha=0.6) +
    scale_fill_manual(values=colors) +
    stat_summary(fun.data=annotate_n, geom="text", fun=median,
                 position=position_dodge(width=0.75), size=4.5) +
    geom_hline(yintercept=0.5, linetype=3, size=0.5) +
    stat_compare_means(aes(label=..p.signif..),
                       test="wilcox.test",
                       comparisons=comparisons) +
    facet_wrap(facet.vars, nrow=2) +
    labs(title=paste(tissue, "biased genes", sep="-"),
         subtitle=species,
         x="Organ bias",
         y=ylab,
         fill="Organ bias") +
    theme_bw() +
    theme(plot.title=element_text(size=20),
          plot.subtitle=element_text(size=18, face="italic"),
          axis.title=element_text(size=18),
          axis.text=element_text(size=14),
          strip.text=element_text(size=14),
          panel.spacing.x=unit(1.1,"lines"),
          legend.position="none")   
}

heatmap_effect_size_pairwise_comparisons <- function(wilcox.df, species, value=c("mean","log2cv","residual_variation","variability"), order.levels) {
  wilcox.df$Bias <- gsub("_", " ", wilcox.df$Bias)
  wilcox.df$Tissue <- gsub("_", " ", wilcox.df$Tissue)
  
  # Reorder factor levels
  wilcox.df$Bias <- factor(wilcox.df$Bias, levels=gsub("_", " ", order.levels))
  wilcox.df$Tissue <- factor(wilcox.df$Tissue, levels=gsub("_", " ", order.levels))
  
  if (value=="mean") {
    title <- "Expression levels"
  } else if (value=="log2cv") { 
    title <- "Log-coefficient of variation"
  } else if (value=="residual_variation") {
    title <- "Residual variation"      
  }else if (value=="variability") {
    title <- "Variability ranks"      
  }
  
  custom_theme <- theme(plot.title=element_text(size=20),
                        plot.subtitle=element_text(size=18, face="italic"),
                        axis.title=element_text(size=18),
                        axis.text=element_text(size=14),
                        strip.text=element_text(size=14),
                        legend.title=element_text(size=14),
                        legend.text=element_text(size=14),
                        panel.grid.major = element_blank())    

  if ("Sex" %in% colnames(wilcox.df)) {
    ggplot(wilcox.df, 
           aes(x=Bias, y=Tissue, fill=GlassDelta)) + 
      geom_tile() +
      scale_fill_viridis(discrete=FALSE, option="inferno", na.value="gray") +
      facet_wrap(vars(Sex), nrow=2) +
      labs(title=title,
           subtitle=species,
           x="Organ-biased genes",
           y="Focal organ",
           fill="Effect size") +
      theme_bw() +
      custom_theme

  } else {
    ggplot(wilcox.df, 
           aes(x=Bias, y=Tissue, fill=GlassDelta)) + 
      geom_tile() +
      scale_fill_viridis(discrete=FALSE, option="inferno", na.value="gray") +
      labs(title=title,
           subtitle=species,
           x="Organ-biased genes",
           y="Focal organ",
           fill="Effect size") +
      theme_bw() +
      custom_theme
  }
}

heatmap_effect_size_pairwise_comparisons_blue2red <- function(wilcox.df, species, value=c("mean","log2cv","residual_variation","variability"), adj.p=0.05, order.levels) {
  wilcox.df$Bias <- gsub("_", " ", wilcox.df$Bias)
  wilcox.df$Tissue <- gsub("_", " ", wilcox.df$Tissue)
  
  # Reorder factor levels
  wilcox.df$Bias <- factor(wilcox.df$Bias, levels=gsub("_", " ", order.levels))
  wilcox.df$Tissue <- factor(wilcox.df$Tissue, levels=gsub("_", " ", order.levels))
  
  if (value=="mean") {
    title <- "Expression levels"
  } else if (value=="log2cv") { 
    title <- "Log-coefficient of variation"
  } else if (value=="residual_variation") {
    title <- "Residual variation"      
  }else if (value=="variability") {
    title <- "Variability ranks"      
  }
  
  custom_theme <- theme(plot.title=element_text(size=20),
                        plot.subtitle=element_text(size=18, face="italic"),
                        axis.title=element_text(size=18),
                        axis.text=element_text(size=14),
                        strip.text=element_text(size=14),
                        legend.title=element_text(size=14),
                        legend.text=element_text(size=14),
                        panel.grid.major = element_blank())    
  
  if ("Sex" %in% colnames(wilcox.df)) {
    ggplot(wilcox.df, 
           aes(x=Bias, y=Tissue, fill=GlassDelta)) + 
      geom_tile() +
      scale_fill_gradient2(midpoint = 0, low="#00BFC4", mid="white", high="#F8766D", na.value="gray") +
      geom_text(data=subset(wilcox.df[wilcox.df$p.adj < adj.p,], p.signif!="ns"),
                size=6, aes(label=p.signif)) +
      facet_wrap(vars(Sex), nrow=2) +
      labs(title=title,
           subtitle=species,
           x="Organ-biased genes",
           y="Focal organ",
           fill="Effect size") +
      theme_bw() +
      custom_theme
    
  } else {
    ggplot(wilcox.df, 
           aes(x=Bias, y=Tissue, fill=GlassDelta)) + 
      geom_tile() +
      scale_fill_gradient2(midpoint = 0, low="#00BFC4", mid="white", high="#F8766D", na.value="gray") +
      geom_text(data=subset(wilcox.df[wilcox.df$p.adj < adj.p,], p.signif!="ns"),
                size=6, aes(label=p.signif)) +
      labs(title=title,
           subtitle=species,
           x="Organ-biased genes",
           y="Focal organ",
           fill="Effect size") +
      theme_bw() +
      custom_theme
  }
}

heatmap_effect_size_combined_comparisons <- function(wilcox.df, cross.species=TRUE, species, value=c("mean","log2cv","residual_variation","variability")) {
  if (value=="mean") {
    title <- "Expression levels"
  } else if (value=="log2cv") { 
    title <- "Log-coefficient of variation"
  } else if (value=="residual_variation") {
    title <- "Residual variation"      
  } else if (value=="variability") {
    title <- "Variability ranks"      
  }
  
  if (cross.species==TRUE) {
    ggplot(wilcox.df, 
           aes(x=Bias, y=Tissue, fill=GlassDelta)) + 
      geom_tile() +
      scale_fill_viridis(discrete=FALSE, option="inferno", na.value="gray") +
      facet_wrap(vars(SpeciesName, Sex), ncol=2, scales="free_x") +
      labs(title=title, 
           x="Organ-biased genes",
           y="Focal organ",
           fill="Effect size") +
      theme_bw() +
      theme(plot.title=element_text(size=24),
            axis.title=element_text(size=18),
            axis.text=element_text(size=14),
            strip.text=element_text(size=18),
            legend.title=element_text(size=18),
            legend.text=element_text(size=12, angle=30),
            panel.grid.major=element_blank(),
            panel.spacing.x=unit(1.1,"lines"),
            panel.spacing.y=unit(1.1,"lines"),
            legend.position="bottom")     
  } else {
    if ("Sex" %in% colnames(wilcox.df)) {
      facet.vars <- vars(Run, Sex)
    } else { facet.vars <- vars(Run) }
    
    ggplot(wilcox.df, 
           aes(x=Bias, y=Tissue, fill=GlassDelta)) + 
      geom_tile() +
      scale_fill_viridis(discrete=FALSE, option="inferno", na.value="gray") +
      facet_wrap(facet.vars, ncol=length(facet.vars), scales="free_x") +
      labs(title=title,
           subtitle=species,
           x="Organ-biased genes",
           y="Focal organ",
           fill="Effect size") +
      theme_bw() +
      theme(plot.title=element_text(size=24),
            plot.subtitle=element_text(size=18, face="italic"),
            axis.title=element_text(size=18),
            axis.text=element_text(size=14),
            strip.text=element_text(size=18),
            legend.title=element_text(size=18),
            legend.text=element_text(size=12, angle=30),
            panel.grid.major=element_blank(),
            panel.spacing.x=unit(1.1,"lines"),
            panel.spacing.y=unit(1.1,"lines"),
            legend.position="bottom")      
  }

}

heatmap_effect_size_combined_comparisons_blue2red <- function(wilcox.df, cross.species=TRUE, species, value=c("mean","log2cv","residual_variation","variability"), adj.p=0.05) {
  if (value=="mean") {
    title <- "Expression levels"
  } else if (value=="log2cv") { 
    title <- "Log-coefficient of variation"
  } else if (value=="residual_variation") {
    title <- "Residual variation"      
  } else if (value=="variability") {
    title <- "Variability ranks"      
  }
  
  if (cross.species==TRUE) {
    ggplot(wilcox.df, 
           aes(x=Bias, y=Tissue, fill=GlassDelta)) + 
      geom_tile() +
      scale_fill_gradient2(midpoint = 0, low="#00BFC4", mid="white", high="#F8766D", na.value="gray") +
      geom_text(data=subset(wilcox.df[wilcox.df$p.adj < adj.p,], p.signif!="ns"),
                size=6, aes(label=p.signif)) +
      facet_wrap(vars(SpeciesName, Sex), ncol=2, scales="free_x") +
      labs(title=title, 
           x="Organ-biased genes",
           y="Focal organ",
           fill="Effect size") +
      theme_bw() +
      theme(plot.title=element_text(size=24),
            axis.title=element_text(size=18),
            axis.text=element_text(size=14),
            strip.text=element_text(size=18),
            legend.title=element_text(size=18),
            legend.text=element_text(size=12, angle=30),
            panel.grid.major=element_blank(),
            panel.spacing.x=unit(1.1,"lines"),
            panel.spacing.y=unit(1.1,"lines"),
            legend.position="bottom")     
  } else {
    if ("Sex" %in% colnames(wilcox.df)) {
      facet.vars <- vars(Run, Sex)
    } else { facet.vars <- vars(Run) }
    
    ggplot(wilcox.df, 
           aes(x=Bias, y=Tissue, fill=GlassDelta)) + 
      geom_tile() +
      scale_fill_gradient2(midpoint = 0, low="#00BFC4", mid="white", high="#F8766D", na.value="gray") +
      geom_text(data=subset(wilcox.df[wilcox.df$p.adj < adj.p,], p.signif!="ns"),
                size=6, aes(label=p.signif)) +
      facet_wrap(facet.vars, ncol=length(facet.vars), scales="free_x") +
      labs(title=title,
           subtitle=species,
           x="Organ-biased genes",
           y="Focal organ",
           fill="Effect size") +
      theme_bw() +
      theme(plot.title=element_text(size=24),
            plot.subtitle=element_text(size=18, face="italic"),
            axis.title=element_text(size=18),
            axis.text=element_text(size=14),
            strip.text=element_text(size=18),
            legend.title=element_text(size=18),
            legend.text=element_text(size=12, angle=30),
            panel.grid.major=element_blank(),
            panel.spacing.x=unit(1.1,"lines"),
            panel.spacing.y=unit(1.1,"lines"),
            legend.position="bottom")     
  }

}