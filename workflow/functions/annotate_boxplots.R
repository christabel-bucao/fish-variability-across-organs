#### Notes ####

#### Libraries ####

#### Sourced functions ####

#### Functions ####

# Annotate boxplot with number of observations per grouping
annotate_n <- function(x) {
  return(c(y = median(x)*1.10, label = length(x))) 
}

# Annotate boxplot with median value per grouping
annotate_median <- function(x) {
  return(c(y = median(x)*1.10, label = round(median(x),2))) 
}