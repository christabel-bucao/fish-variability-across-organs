#### Notes ####
# Functions to compute expression variability
# See 02_Variability_Jackknife.Rmd for sample usage

#### Libraries ####
library(dplyr)
library(matrixStats)
library(slider)

#### Functions ####

#### Compute genewise summary statistics ####

# Input: expr.matrix = expression matrix (rows - genes, columns - samples)
# Output: stats.df = data frame with summary statistics per gene
compute_gene_summary_stats <- function(expr.matrix) {
  # Change Inf to NA (if any)
  expr.matrix[is.infinite(expr.matrix)] <- NA
  
  means <- rowMeans(expr.matrix, na.rm=TRUE)
  medians <- rowMedians(expr.matrix, useNames=TRUE, na.rm=TRUE)
  vars <- rowVars(expr.matrix, useNames=TRUE, na.rm=TRUE)
  sds <- rowSds(expr.matrix, useNames=TRUE, na.rm=TRUE)
  cv.sq <- (sds/means)^2
  
  stats.df <- 
    data.frame("Mean"=means, "Median"=medians, 
               "Var"=vars, "SD"=sds, "CV2"=cv.sq)
  stats.df$Rank_Mean <- percent_rank(stats.df$Mean)
  
  return(stats.df)
}

#### Functions to compute expression variability ####

# Adjusted SD - Defined as the ratio of observed to predicted SD
# Metric adapted from Liu et al. (2020) (doi: 10.1186/s12915-020-00842-z), but using LOESS regression
# Input: expr.matrix = expression matrix
# Output: model.results = list with input (mean, SD) and output values (adj SD, log2 adj SD)
adjusted_sd <- function(expr.matrix) {
  means <- rowMeans(expr.matrix)
  sds <- rowSds(expr.matrix, useNames=TRUE)
  
  model.df <- data.frame("Mean"=means, "SD"=sds)
  
  model.fit <- loess(SD~Mean, model.df, span=0.6)
  adj.sd <- model.df$SD/model.fit$fitted
  names(adj.sd) <- rownames(model.df)
  
  model.results <- 
    list("mean"=means, "sd"=sds, 
         "adj.sd"=adj.sd, "log2.adj.sd"=log2(adj.sd))
  
  return(model.results)  
}

# Defined as the residual of log2 of the (squared) coefficient of variation (CV)
# Metric adapted from Faure et al. (2017) (doi: 10.1016/j.cels.2017.10.003)
residual_log2cv <- function(expr.matrix) {
  means <- rowMeans(expr.matrix)
  sds <- rowSds(expr.matrix, useNames=TRUE)
  cv.sq <- (sds/means)^2
  log2.cv.sq <- log2(cv.sq)
  
  model.df <- data.frame("Mean"=means, "Log2CV"=log2.cv.sq)
  model.fit <- loess(Log2CV~Mean, model.df, span=0.6)
  
  model.results <- 
    list("mean"=means, "log2.cv"=log2.cv.sq, "resid.log2.cv"=model.fit$residuals)
  
  return(model.results)
}

# Defined as the residual of log2 SD
residual_log2sd <- function(expr.matrix) {
  means <- rowMeans(expr.matrix)
  log2.sd <-log2(rowSds(expr.matrix, useNames=TRUE))
  
  model.df <- data.frame("Mean"=means, "Log2SD"=log2.sd)
  model.fit <- loess(Log2SD~Mean, model.df, span=0.6)
  
  model.results <- 
    list("mean"=means, "log2.sd"=log2.sd, "resid.log2.sd"=model.fit$residuals)
  
  return(model.results)
}

#### Compute local variability percentile rank using sliding windows ####

# This controls for dependency between mean expression and variance of variability
# Adapted from Faure et al. (2017) (doi: 10.1016/j.cels.2017.10.003)

# Input: 
## jack.mean = vector of mean expression per gene
## jack.ev = vector of computed variability per gene
## win.size = integer length of sliding window
# Output:
## local.rank = numeric vector of local variability rank (within [0, 1]) per gene
compute_local_ev_rank <- function(jack.mean, jack.ev, win.size=100) {
  ### Sort genes by expression level
  sorted.jack.mean <- sort(jack.mean)
  sorted.jack.ev <- jack.ev[names(sorted.jack.mean)]
  
  # Set windows
  windows <- 
    slider::slide(sorted.jack.ev, ~.x, 
                  .before=(win.size/2), .after=(win.size/2), .complete=TRUE)
  
  # Compute local variability rank for each window
  local.rank <- lapply(windows, dplyr::percent_rank)
  local.rank <- unlist(lapply(local.rank, function(i) { i[((win.size/2)+1)] }))
  names(local.rank) <- names(sorted.jack.ev)
  local.rank <- local.rank[names(jack.ev)]
  
  return(local.rank)
}

#### Functions for jackknife analysis ####

# For each tissue or tissue-sex condition with n replicates
# (1) Consider each subset of n-1 replicates
# (1) Compute variability and local variability rank of genes for each subset
# (3) Use the average of the results to have a jackknifed estimate of variability

# Input:
# Input:
## expr.matrix = expression matrix
## min.percentile = numeric between [0,1]; removes bottom x% of genes by expression level
## max.percentile = numeric between [0,1]; removes top (1-x)% of genes 
## (e.g. max.percentile=0.95 removes top 5% of genes by expression)
## win.size = integer length of sliding window to be used for compute_local_ev_rank


# Output:
## jack = includes list with matrices of computed results per (n-1) subset
## and summary data frame with average of jackknifed results

# Adjusted SD
jackknife_adj_sd <- function(expr.matrix, min.percentile=0.00, max.percentile=0.95, win.size=100) {
  # Remove lowly expressed genes
  filtered.matrix <- expr.matrix[rowMedians(expr.matrix)>1.0 &
                                   rowMeans(expr.matrix)>1.0,]
  
  jack <- list()
  jack$mean <- matrix(nrow=nrow(filtered.matrix), ncol=ncol(filtered.matrix),
                      dimnames=list(rownames(filtered.matrix)))
  
  jack$sd <- matrix(nrow=nrow(filtered.matrix), ncol=ncol(filtered.matrix),
                    dimnames=list(rownames(filtered.matrix)))
  jack$adj.sd <- matrix(nrow=nrow(filtered.matrix), ncol=ncol(filtered.matrix),
                        dimnames=list(rownames(filtered.matrix)))
  jack$rank.adj.sd <- matrix(nrow=nrow(filtered.matrix), ncol=ncol(filtered.matrix),
                              dimnames=list(rownames(filtered.matrix)))
  
  # Calculate statistics per subset of n-1 samples
  for (i in 1:ncol(filtered.matrix)) {
    print(paste(i, "of", ncol(filtered.matrix)))
    
    # Variability metric
    adj.sd.results <- adjusted_sd(filtered.matrix[,-i])
    
    jack$mean[,i] <- adj.sd.results$mean
    jack$sd[,i] <- adj.sd.results$sd
    jack$adj.sd[,i] <- adj.sd.results$adj.sd
    
    # Local variability rank within a sliding window
    jack$rank.adj.sd[,i] <- 
      compute_local_ev_rank(jack$mean[,i], jack$adj.sd[,i], win.size)
  }
  
  # Summarize jackknife results
  jack$summary <- 
    data.frame("Mean_Mean"=rowMeans(jack$mean), # Jackknifed mean
               "Median_Mean"=rowMedians(jack$mean, useNames=TRUE),
               "Mean_SD"=rowMeans(jack$sd),
               "Median_SD"=rowMedians(jack$sd, useNames=TRUE),
               "Mean_AdjSD"=rowMeans(jack$adj.sd),
               "Median_AdjSD"=rowMedians(jack$adj.sd, useNames=TRUE),
               row.names=rownames(filtered.matrix))
  
  # Remove top and bottom x% of genes by expression
  jack$summary$Rank_Mean <- dplyr::percent_rank(jack$summary$Mean_Mean)
  jack$summary <- jack$summary[jack$summary$Rank_Mean<max.percentile &
                                 jack$summary$Rank_Mean>min.percentile,]
  
  # Local and global variability rank
  mean.rank.adj.sd <- rowMeans(jack$rank.adj.sd, na.rm=TRUE)
  median.rank.adj.sd <- rowMedians(jack$rank.adj.sd, na.rm=TRUE, useNames=TRUE)
  
  jack$summary$Mean_Local_Rank_AdjSD <- mean.rank.adj.sd[rownames(jack$summary)]
  jack$summary$Median_Local_Rank_AdjSD <- median.rank.adj.sd[rownames(jack$summary)]
  
  # Global rank - no sliding window applied
  jack$summary$Global_Rank_AdjSD <- dplyr::percent_rank(jack$summary$Median_AdjSD)

  # Complete cases only
  jack$summary <- jack$summary[complete.cases(jack$summary),]
  
  return(jack)
}

## Residual log2 (CV)^2
jackknife_resid_log2cv <- function(expr.matrix, min.percentile=0.00, max.percentile=0.95, win.size=100) {
  # Remove lowly expressed genes
  filtered.matrix <- expr.matrix[rowMedians(expr.matrix)>1.0 &
                                   rowMeans(expr.matrix)>1.0,]
  
  jack <- list()
  jack$mean <- matrix(nrow=nrow(filtered.matrix), ncol=ncol(filtered.matrix),
                      dimnames=list(rownames(filtered.matrix)))
  jack$log2.cv <- matrix(nrow=nrow(filtered.matrix), ncol=ncol(filtered.matrix),
                         dimnames=list(rownames(filtered.matrix)))
  jack$resid.log2.cv <- matrix(nrow=nrow(filtered.matrix), ncol=ncol(filtered.matrix),
                               dimnames=list(rownames(filtered.matrix)))
  jack$rank.resid.log2.cv <- matrix(nrow=nrow(filtered.matrix), ncol=ncol(filtered.matrix),
                                     dimnames=list(rownames(filtered.matrix)))
  
  # Calculate statistics per subset of n-1 samples
  for (i in 1:ncol(filtered.matrix)) {
    print(paste(i, "of", ncol(filtered.matrix)))
    
    # Variability metric
    log2.cv.results <- residual_log2cv(filtered.matrix[,-i])
    
    jack$mean[,i] <- log2.cv.results$mean
    jack$log2.cv[,i] <- log2.cv.results$log2.cv
    jack$resid.log2.cv[,i] <- log2.cv.results$resid.log2.cv
    
    # Local variability rank within a sliding window
    jack$rank.resid.log2.cv[,i] <- 
      compute_local_ev_rank(jack$mean[,i], jack$resid.log2.cv[,i], win.size)
  }
  
  # Summarize jackknife results
  jack$summary <- 
    data.frame("Mean_Mean"=rowMeans(jack$mean), # Jackknifed mean
               "Median_Mean"=rowMedians(jack$mean, useNames=TRUE),
               "Mean_Log2CV"=rowMeans(jack$log2.cv),
               "Median_Log2CV"=rowMedians(jack$log2.cv, useNames=TRUE),
               "Mean_Resid_Log2CV"=rowMeans(jack$resid.log2.cv),
               "Median_Resid_Log2CV"=rowMedians(jack$resid.log2.cv, useNames=TRUE),
               row.names=rownames(filtered.matrix))
  
  # Remove top and bottom x% of genes by expression
  jack$summary$Rank_Mean <- dplyr::percent_rank(jack$summary$Mean_Mean)
  jack$summary <- jack$summary[jack$summary$Rank_Mean<max.percentile &
                                 jack$summary$Rank_Mean>min.percentile,]
  
  # Local and global variability rank
  mean.local.rank.log2.cv <- rowMeans(jack$rank.resid.log2.cv, na.rm=TRUE)
  median.local.rank.log2.cv <- 
    rowMedians(jack$rank.resid.log2.cv, na.rm=TRUE, useNames=TRUE)
  
  jack$summary$Mean_Local_Rank_Log2CV <- mean.local.rank.log2.cv[rownames(jack$summary)]
  jack$summary$Median_Local_Rank_Log2CV <- median.local.rank.log2.cv[rownames(jack$summary)]
  
  # Global rank - no sliding window applied
  jack$summary$Global_Rank_Log2CV <- dplyr::percent_rank(jack$summary$Median_Resid_Log2CV)
  
  # Complete cases only
  jack$summary <- jack$summary[complete.cases(jack$summary),]
  
  return(jack)  
}

# Residual log2 SD
jackknife_resid_log2sd <- function(expr.matrix, min.percentile=0.00, max.percentile=0.95, win.size=100) {
  # Remove lowly expressed genes
  filtered.matrix <- expr.matrix[rowMedians(expr.matrix)>1.0 &
                                   rowMeans(expr.matrix)>1.0,]
  
  jack <- list()
  jack$mean <- matrix(nrow=nrow(filtered.matrix), ncol=ncol(filtered.matrix),
                      dimnames=list(rownames(filtered.matrix)))
  jack$log2.sd <- matrix(nrow=nrow(filtered.matrix), ncol=ncol(filtered.matrix),
                         dimnames=list(rownames(filtered.matrix)))
  jack$resid.log2.sd <- matrix(nrow=nrow(filtered.matrix), ncol=ncol(filtered.matrix),
                               dimnames=list(rownames(filtered.matrix)))
  jack$rank.resid.log2.sd <- matrix(nrow=nrow(filtered.matrix), ncol=ncol(filtered.matrix),
                                     dimnames=list(rownames(filtered.matrix)))
  
  # Calculate statistics per subset of n-1 samples
  for (i in 1:ncol(filtered.matrix)) {
    print(paste(i, "of", ncol(filtered.matrix)))
    
    # Variability metric
    log2.sd.results <- residual_log2sd(filtered.matrix[,-i])
    
    jack$mean[,i] <- log2.sd.results$mean
    jack$log2.sd[,i] <- log2.sd.results$log2.sd
    jack$resid.log2.sd[,i] <- log2.sd.results$resid.log2.sd
    
    # Local variability rank within a sliding window
    jack$rank.resid.log2.sd[,i] <- 
      compute_local_ev_rank(jack$mean[,i], jack$resid.log2.sd[,i], win.size)
  }
  
  # Summarize jackknife results
  jack$summary <- 
    data.frame("Mean_Mean"=rowMeans(jack$mean), # Jackknifed mean
               "Median_Mean"=rowMedians(jack$mean, useNames=TRUE),
               "Mean_Log2SD"=rowMeans(jack$log2.sd),
               "Median_Log2SD"=rowMedians(jack$log2.sd, useNames=TRUE),
               "Mean_Resid_Log2SD"=rowMeans(jack$resid.log2.sd),
               "Median_Resid_Log2SD"=rowMedians(jack$resid.log2.sd, useNames=TRUE),
               row.names=rownames(filtered.matrix))
  
  # Remove top and bottom x% of genes by expression
  jack$summary$Rank_Mean <- dplyr::percent_rank(jack$summary$Mean_Mean)
  jack$summary <- jack$summary[jack$summary$Rank_Mean<max.percentile &
                                 jack$summary$Rank_Mean>min.percentile,]
  
  # Local and global variability rank
  mean.local.rank.log2.sd <- rowMeans(jack$rank.resid.log2.sd, na.rm=TRUE)
  median.local.rank.log2.sd <- 
    rowMedians(jack$rank.resid.log2.sd, na.rm=TRUE, useNames=TRUE)
  
  jack$summary$Mean_Local_Rank_Log2SD <- mean.local.rank.log2.sd[rownames(jack$summary)]
  jack$summary$Median_Local_Rank_Log2SD <- median.local.rank.log2.sd[rownames(jack$summary)]
  
  # Global rank - no sliding window applied
  jack$summary$Global_Rank_Log2SD <- dplyr::percent_rank(jack$summary$Median_Resid_Log2SD)
  
  # Complete cases only
  jack$summary <- jack$summary[complete.cases(jack$summary),]
  
  return(jack)  
}