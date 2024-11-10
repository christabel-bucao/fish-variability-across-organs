#### Notes ####

#### Libraries ####
library(dplyr)
library(tidyr)
library(ggpubr)

#### Sourced functions ####

#### Functions ####

# Classify genes by organ bias based on tau cutoff
classify_genes_by_organ_bias <- function(df, tau.cutoff=0.30) {
  df.with.bias <- df
  
  # Broadly expressed, focal-biased, t-biased, and other-biased genes
  df.with.bias$TissueBias <-
    apply(df, MARGIN=1, function(r) {
      if (r[["tau"]] <= tau.cutoff) { "broad" }
      else if (r[["Tissue"]]==r[["PrimaryTissue"]]) { "focal" }
      else if (r[["PrimaryTissue"]]==t) { r[["PrimaryTissue"]] }
      else { "other" }
    })
  
  return(df.with.bias)
}

# Compute effect size (Glass's delta)
compute_glass_delta_effect_size <- function(comparisons.df, value=c("mean","log2cv","residual_variation","variability")) {
  if ("Sex" %in% colnames(comparisons.df)) {
    grouping <- c("Bias", "Tissue", "Sex")  
  } else { grouping <- c("Bias", "Tissue") }
  
  if (value=="mean") {
    value.comparisons.summary <- comparisons.df %>%
      dplyr::group_by(across(c(all_of(grouping),"TissueBias"))) %>%
      dplyr::summarise(
        Mean_Value = mean(Mean_Mean),
        SD_Value = sd(Mean_Mean),
      ) %>%
      dplyr::ungroup()
    
  } else if (value=="log2cv") {
    value.comparisons.summary <- comparisons.df %>%
      dplyr::group_by(across(c(all_of(grouping),"TissueBias"))) %>%
      dplyr::summarise(
        Mean_Value = mean(Mean_Log2CV),
        SD_Value = sd(Mean_Log2CV),
      ) %>%
      dplyr::ungroup()
    
  } else if (value=="residual_variation") {
    value.comparisons.summary <- comparisons.df %>%
      dplyr::group_by(across(c(all_of(grouping),"TissueBias"))) %>%
      dplyr::summarise(
        Mean_Value = mean(Mean_Resid_Log2CV),
        SD_Value = sd(Mean_Resid_Log2CV),
      ) %>%
      dplyr::ungroup()
    
  } else if (value=="variability") {
    value.comparisons.summary <- comparisons.df %>%
      dplyr::group_by(across(c(all_of(grouping),"TissueBias"))) %>%
      dplyr::summarise(
        Mean_Value = mean(Mean_Local_Rank_Log2CV),
        SD_Value = sd(Mean_Local_Rank_Log2CV),
      ) %>%
      dplyr::ungroup()    
  }
  
  # Separate foreground gene set and background ('other')
  mean.value.foreground <- value.comparisons.summary %>%
    dplyr::filter(TissueBias!="other")
  mean.value.background <- value.comparisons.summary %>%
    dplyr::filter(TissueBias=="other")
  
  # Difference in group means
  diff.mean.value <-
    dplyr::inner_join(
      mean.value.foreground[,c(grouping, "Mean_Value", "SD_Value")],
      mean.value.background[,c(grouping, "Mean_Value", "SD_Value")],
      by=grouping,
      suffix=c("_Foreground","_Other"))
  
  diff.mean.value$DiffMeans_Value <-
    diff.mean.value$Mean_Value_Foreground - diff.mean.value$Mean_Value_Other
  
  # Effect size (Glass's delta)
  ## Defined as the difference of the group means divided by the standard deviation of the reference group ('other')
  effect.size.value <- diff.mean.value
  effect.size.value$GlassDelta <-
    diff.mean.value$DiffMeans_Value / diff.mean.value$SD_Value_Other
  
  # Rename columns
  if (value=="mean") {
    colnames(effect.size.value) <- gsub("Value","MeanExpr",colnames(effect.size.value))
    
  } else if (value=="log2cv") { 
    colnames(effect.size.value) <- gsub("Value","MeanLog2CV",colnames(effect.size.value))
    
  } else if (value=="residual_variation") {
    colnames(effect.size.value) <- gsub("Value","MeanResidVar",colnames(effect.size.value))
    
  } else if (value=="variability") {
    colnames(effect.size.value) <- gsub("Value","MeanVarRank",colnames(effect.size.value))
  }
  
  return(effect.size.value)
}

# Pairwise Wilcoxon test
compute_pairwise_comparisons <- function(comparisons.df, effect.size.df, value=c("mean","log2cv","residual_variation","variability"), apply.pvalue.cutoff=TRUE, adj.p=0.05, order.levels) {
  if ("Sex" %in% colnames(comparisons.df)) {
    grouping <- c("Bias", "Tissue", "Sex")  
  } else { grouping <- c("Bias", "Tissue") }
  
  if (value=="mean") {
    value.wilcox.df <- 
      compare_means(formula=Mean_Mean~TissueBias,
                    data=comparisons.df,
                    group.by=grouping,
                    method="wilcox.test",
                    p.adjust.method="BH",
                    ref.group="other")
    
  } else if (value=="log2cv") {
    value.wilcox.df <- 
      compare_means(formula=Mean_Log2CV~TissueBias,
                    data=comparisons.df,
                    group.by=grouping,
                    method="wilcox.test",
                    p.adjust.method="BH",
                    ref.group="other")    
    
  } else if (value=="residual_variation") {
    value.wilcox.df <- 
      compare_means(formula=Mean_Resid_Log2CV~TissueBias,
                    data=comparisons.df,
                    group.by=grouping,
                    method="wilcox.test",
                    p.adjust.method="BH",
                    ref.group="other")      

  }
  else if (value=="variability") {
    value.wilcox.df <- 
      compare_means(formula=Mean_Local_Rank_Log2CV~TissueBias,
                    data=comparisons.df,
                    group.by=grouping,
                    method="wilcox.test",
                    p.adjust.method="BH",
                    ref.group="other")
  }
  
  # Add effect size
  effect.size.value.wilcox.df <- value.wilcox.df %>%
    dplyr::inner_join(effect.size.df, by=grouping)
  
  # Add missing comparisons (for heatmap visualization)
  if ("Sex" %in% colnames(comparisons.df)) {
    effect.size.value.wilcox.complete <- effect.size.value.wilcox.df %>%
      tidyr::complete(Bias, Tissue, Sex)
  } else {
    effect.size.value.wilcox.complete <- effect.size.value.wilcox.df %>%
      tidyr::complete(Bias, Tissue)    
  }
  
  # Retain comparisons that pass the adjusted p-value cutoff
  if (apply.pvalue.cutoff==TRUE) {
    effect.size.value.wilcox.complete.signif <- effect.size.value.wilcox.complete %>%
      dplyr::filter(is.na(effect.size.value.wilcox.complete$p.adj) | 
                      effect.size.value.wilcox.complete$p.adj < adj.p)      
    
    # Reorder levels (for heatmap visualization)
    effect.size.value.wilcox.complete.signif$Bias <- 
      factor(effect.size.value.wilcox.complete.signif$Bias, levels=order.levels)
    
    return(effect.size.value.wilcox.complete.signif)
    
  } else {
    # Reorder levels (for heatmap visualization)
    effect.size.value.wilcox.complete$Bias <- 
      factor(effect.size.value.wilcox.complete$Bias, levels=order.levels)
    
    return(effect.size.value.wilcox.complete)
  }
}

bind_pairwise_comparisons <- function(comparisons.list, apply.pvalue.cutoff=TRUE, adj.p=0.05) {
  # No adjusted p-value cutoff
  wilcox.complete.bind <-
    dplyr::bind_rows(list("LOC"=comparisons.list[["LOC"]]$wilcox.complete,
                          "ELU"=comparisons.list[["ELU"]]$wilcox.complete,
                          "DRE"=comparisons.list[["DRE"]]$wilcox.complete),
                     .id="Species") %>%
    tidyr::complete(Species, Bias, Tissue, Sex) %>% # Add missing comparisons
    dplyr::filter(!(Species=="LOC" & !is.na(Sex))) %>% # Remove some comparisons
    dplyr::filter(!(Species=="ELU" & is.na(Sex))) %>%
    dplyr::filter(!(Species=="DRE" & is.na(Sex)))
  
  # Add species names
  wilcox.complete.bind$SpeciesName <-
    vapply(wilcox.complete.bind$Species, FUN.VALUE=character(1), function(sp) {
      if (sp=="LOC") { "Spotted gar" }
      else if (sp=="ELU") { "Northern pike" }
      else if (sp=="DRE") { "Zebrafish" }
    })
  wilcox.complete.bind$SpeciesName <-
    factor(wilcox.complete.bind$SpeciesName,
           levels=c("Zebrafish", "Northern pike", "Spotted gar"))
  
  # Rearrange columns
  wilcox.complete.bind <- 
    wilcox.complete.bind[,c("SpeciesName",colnames(wilcox.complete.bind)[1:ncol(wilcox.complete.bind)-1])]
  
  # Unknown sex for spotted gar
  wilcox.complete.bind$Sex[is.na(wilcox.complete.bind$Sex) & wilcox.complete.bind$Species=="LOC"] <- "unknown"
  
  # Remove underscores (for heatmap visualization)
  wilcox.complete.bind$Bias <- gsub("_", " ", wilcox.complete.bind$Bias)
  wilcox.complete.bind$Tissue <- gsub("_", " ", wilcox.complete.bind$Tissue)
  
  # Reorder levels (for heatmap visualization)
  levels.order <- 
    c("brain","eye","gills","pectoral fin", "intestine","liver","heart","muscle","ovary","testis","gonads")
  
  wilcox.complete.bind$Bias <- factor(wilcox.complete.bind$Bias, levels=levels.order)
  wilcox.complete.bind$Tissue <- factor(wilcox.complete.bind$Tissue, levels=levels.order)
  
  if (apply.pvalue.cutoff==TRUE) {
    # Retain comparisons that pass the adjusted p-value cutoff
    wilcox.signif.bind <- wilcox.complete.bind %>%
      dplyr::filter(is.na(wilcox.complete.bind$p.adj) | wilcox.complete.bind$p.adj < adj.p)
    
    # Reorder levels (for heatmap visualization)
    wilcox.signif.bind$Bias <- factor(wilcox.signif.bind$Bias, levels=levels.order)
    wilcox.signif.bind$Tissue <- factor(wilcox.signif.bind$Tissue, levels=levels.order)
    
    return(wilcox.signif.bind)
  } else { return(wilcox.complete.bind) }

}

bind_pairwise_comparisons_observed_vs_random <- function(comparisons.list, apply.pvalue.cutoff=TRUE, adj.p=0.05) {
  # No adjusted p-value cutoff
  wilcox.complete.bind <-
    dplyr::bind_rows(list(comparisons.list[[1]]$wilcox.complete,
                          comparisons.list[[2]]$wilcox.complete,
                          comparisons.list[[3]]$wilcox.complete),
                     .id="Run")
  
  wilcox.complete.bind$Run <-
    vapply(wilcox.complete.bind$Run, FUN.VALUE=character(1), function(r) {
      if (r==1) { names(comparisons.list)[1] }
      else if (r==2) { names(comparisons.list)[2] }
      else if (r==3) { names(comparisons.list)[3] }
    })
    
  if ("Sex" %in% colnames(wilcox.complete.bind)) {
    wilcox.complete.bind <- wilcox.complete.bind %>%
      tidyr::complete(Run, Bias, Tissue, Sex) # Add missing comparisons      
  } else {
    wilcox.complete.bind <- wilcox.complete.bind %>%
      tidyr::complete(Run, Bias, Tissue)
  }

  # Remove underscores (for heatmap visualization)
  wilcox.complete.bind$Bias <- gsub("_", " ", wilcox.complete.bind$Bias)
  wilcox.complete.bind$Tissue <- gsub("_", " ", wilcox.complete.bind$Tissue)
  
  # Reorder levels (for heatmap visualization)
  levels.order <- 
    c("brain","eye","gills","pectoral fin", "intestine","liver","heart","muscle","ovary","testis","gonads")
  
  ## Exclude missing conditions from factor levels
  levels.exclude <- setdiff(levels.order, levels(as.factor(wilcox.complete.bind$Bias)))
  
  wilcox.complete.bind$Bias <-
    factor(wilcox.complete.bind$Bias, levels=levels.order, exclude=levels.exclude)
  
  wilcox.complete.bind$Tissue <-
    factor(wilcox.complete.bind$Tissue, levels=levels.order, exclude=levels.exclude)
  
  if (apply.pvalue.cutoff==TRUE) {
    # Retain comparisons that pass the adjusted p-value cutoff
    wilcox.signif.bind <- wilcox.complete.bind %>%
      dplyr::filter(is.na(wilcox.complete.bind$p.adj) | wilcox.complete.bind$p.adj < adj.p)
    
    # Reorder levels (for heatmap visualization)
    wilcox.signif.bind$Bias <- factor(wilcox.signif.bind$Bias, levels=levels.order)
    wilcox.signif.bind$Tissue <- factor(wilcox.signif.bind$Tissue, levels=levels.order)
    
    return(wilcox.signif.bind)
  } else { return(wilcox.complete.bind) }
  
}