#### Notes ####

#### Libraries ####
library(colorBlindness)

#### Sourced functions ####

#### Functions ####

# Permute values in each column
permute_by_cols <- function(input.matrix, seed=12345) {
  set.seed(seed)
  
  permuted.matrix <- input.matrix
  for (s in colnames(input.matrix)) { 
    permuted.matrix[,s] <- sample(input.matrix[,s], replace=FALSE) 
  }
  
  return(permuted.matrix)
}

# Permute values in each row
permute_by_rows <- function(input.matrix, seed=12345) {
  set.seed(seed)
  
  permuted.matrix <- input.matrix
  for (g in rownames(input.matrix)) {
    permuted.matrix[g,] <- sample(input.matrix[g,], replace=FALSE)
  }
  
  return(permuted.matrix)
}

# Subset expression matrix by condition
subset_expr_matrix <- function(expr.matrix, tissue=NULL, sex=NULL) {
  if (is.null(tissue) & is.null(sex)) { return(expr.matrix) }
  else if (is.null(tissue) & !is.null(sex)) { regexpr <- paste0(sex,"[0-9][0-9]","_") }
  else if (!is.null(tissue) & is.null(sex)) { regexpr <- tissue }
  else if (!is.null(tissue) & !is.null(sex)) { regexpr <- paste0(sex,"[0-9][0-9]","_",tissue)}
  
  sub.matrix <- expr.matrix[,grepl(regexpr, colnames(expr.matrix))]
  return(sub.matrix)
}

## Consider testis and ovary separately
split_gonads_df <- function(metadata) {
    metadata$Tissue[metadata$Tissue=="gonads" & metadata$Sex=="F"] <- "ovary"
    metadata$Tissue[metadata$Tissue=="gonads" & metadata$Sex=="M"] <- "testis"
    
    if ("SampleName" %in% colnames(metadata)) {
      metadata$SampleName <- paste(metadata$SexID, metadata$Tissue, sep="_")    
    }
    
    return(metadata)    
}

split_gonads_matrix <- function(expr.matrix) {
  colnames(expr.matrix)[grepl("F[0-9][0-9]_gonads",colnames(expr.matrix))] <- 
    gsub("gonads", "ovary", 
         colnames(expr.matrix[,grepl("F[0-9][0-9]_gonads",colnames(expr.matrix))]))
  colnames(expr.matrix)[grepl("M[0-9][0-9]_gonads",colnames(expr.matrix))] <- 
    gsub("gonads", "testis",
         colnames(expr.matrix[,grepl("M[0-9][0-9]_gonads",colnames(expr.matrix))]))
  
  return(expr.matrix)    
}

# Consider testis and ovary as one tissue (gonads)
join_gonads_df <- function(metadata=NULL, expr.matrix=NULL) {
  metadata$Tissue[metadata$Tissue=="ovary"] <- "gonads"
  metadata$Tissue[metadata$Tissue=="testis"] <- "gonads"
  
  if ("SampleName" %in% colnames(metadata)) {
    metadata$SampleName <- paste(metadata$SexID, metadata$Tissue, sep="_")    
  }
  
  return(metadata)
}

join_gonads_matrix <- function(expr.matrix) {
  colnames(expr.matrix)[grepl("F[0-9][0-9]_ovary",colnames(expr.matrix))] <- 
    gsub("ovary", "gonads", 
         colnames(expr.matrix[,grepl("F[0-9][0-9]_ovary",colnames(expr.matrix))]))
  colnames(expr.matrix)[grepl("M[0-9][0-9]_testis",colnames(expr.matrix))] <- 
    gsub("testis", "gonads",
         colnames(expr.matrix[,grepl("M[0-9][0-9]_testis",colnames(expr.matrix))]))
  
  return(expr.matrix)    
}

# Set fixed tissue colors for consistency across species. Uses colorblind-friendly palette.
set_tissue_colors <- function() {
  all.tissues <-
    c("blood", "brain", "eye", "gills", "gonads", "heart", "intestine", "kidney",
      "liver", "muscle", "pectoral_fin", "skin", "swim_bladder")
  
  #all.tissue.colors <- c("#D01C8B", RColorBrewer::brewer.pal(n=12, name="Paired"))
  all.tissue.colors <- colorBlindness::paletteMartin[2:14]
  names(all.tissue.colors) <- all.tissues
  
  return(all.tissue.colors)
}