#### Notes ####

#### Libraries ####
library(matrixStats)

#### Sourced functions ####

#### Functions ####

# Convert expression matrix to z-score matrix
zscore <- function(expr.matrix) {
  means <- rowMeans(expr.matrix, na.rm=TRUE)
  sds <- rowSds(expr.matrix, useNames=TRUE, na.rm=TRUE)
  
  zscore.matrix <- (expr.matrix - means) / sds
  
  return(zscore.matrix)
}