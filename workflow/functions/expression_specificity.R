#### Notes ####
# Functions for computing expression specificity

#### Usage ####

# Step 1: Compute mean expression per tissue (if multiple samples per tissue)
## mean.expr <- mean_expr_per_tissue(expr.matrix, metadata)
# Step 2: Compute tau from mean expression
## tau <- compute_tau(mean.expr)

#### Sourced functions ####

#### Functions ####

# Compute mean expression per tissue
mean_expr_per_tissue <- function(expr.matrix, metadata) {
  Tissue <- as.factor(metadata$Tissue)
  
  mean.expr <-
    matrix(nrow=nrow(expr.matrix), ncol=length(levels(Tissue)), 
           dimnames=list(rownames(expr.matrix), levels(Tissue)))
  
  for (t in levels(Tissue)) {
    for (g in rownames(expr.matrix)) {
      mean.expr[g,t] <- 
        mean(expr.matrix[g,grepl(t,colnames(expr.matrix)), drop=FALSE])
    }
  }
  
  return(mean.expr)
}

# Scale by maximum expression
scale_by_max_expr <- function(expr.matrix) {
  max <- apply(expr.matrix, MARGIN=1, max)
  
  scaled.expr <- expr.matrix / max
  names(scaled.expr) <- rownames(expr.matrix)
  
  return(scaled.expr)
}

# Compute tau
compute_tau <- function(expr.matrix) {
  # Takes mean expression matrix
  tau <- apply(expr.matrix, MARGIN=1, function(row) {
    max <- max(row)
    exp.norm <- row/max
    numerator <- sum(1-exp.norm)
    conditions <- ncol(expr.matrix)
    
    numerator/(conditions-1)
  })
  
  names(tau) <- rownames(expr.matrix)
  return(tau)
}