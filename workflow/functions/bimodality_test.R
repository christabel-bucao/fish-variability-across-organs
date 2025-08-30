#### Notes ####

#### Libraries ####
library(matrixStats)

#### Sourced functions ####

#### Functions ####

# Convert expression matrix to z-score matrix
compute_zscore <- function(expr.matrix) {
  means <- rowMeans(expr.matrix, na.rm=TRUE)
  sds <- matrixStats::rowSds(expr.matrix, useNames=TRUE, na.rm=TRUE)
  
  zscore.matrix <- (expr.matrix - means) / sds
  
  return(zscore.matrix)
}