#### Notes ####

#### Libraries ####
library(tibble)
library(dplyr)

#### Sourced functions ####

#### Functions ####

# Merge variability rank
merge_ev_rank <- function(jack.summary, metadata, by=c("local","global"), join=c("inner","full"), method=c("sd","cv")) {
  Tissue <- as.factor(metadata$Tissue)
  if ("Sex" %in% colnames(metadata)) { 
    Sex <- as.factor(metadata$Sex)
    
    for (t in levels(Tissue)) {
      for (s in levels(Sex)) {
        if (!is.null(jack.summary[[t]][[s]])) {
          colnames(jack.summary[[t]][[s]])[-1] <- 
            paste(colnames(jack.summary[[t]][[s]][-1]),t,s,sep="_") 
        } else {
          jack.summary[[t]][[s]] <- NULL
        }
      }
    }
    
    if (join=="inner") {
      rank.ev.matrix <- 
        Reduce(function(...) merge(..., by="GeneID"), unlist(jack.summary, recursive=FALSE))      
    } else if (join=="full") {
      rank.ev.matrix <- 
        Reduce(function(...) merge(..., by="GeneID", all=TRUE), unlist(jack.summary, recursive=FALSE))        
    }
    
    
  } else {
    for (t in levels(Tissue)) {
      colnames(jack.summary[[t]])[-1] <- 
        paste(colnames(jack.summary[[t]][-1]),t,sep="_")
    }
    
    if (join=="inner") {
      rank.ev.matrix <- 
        Reduce(function(...) merge(..., by="GeneID"), jack.summary)      
    } else if (join=="full") {
      rank.ev.matrix <- 
        Reduce(function(...) merge(..., by="GeneID", all=TRUE), jack.summary)      
    }
  }
  
  if (by=="local") {
    if (method=="sd") {
      rank.ev.matrix <- rank.ev.matrix[,grepl("GeneID|Mean_Local_Rank_Log2SD",colnames(rank.ev.matrix))]          
    } else if (method=="cv") {
      rank.ev.matrix <- rank.ev.matrix[,grepl("GeneID|Mean_Local_Rank_Log2CV",colnames(rank.ev.matrix))]      
    }
  } else if (by=="global") {
    if (method=="sd") {
      rank.ev.matrix <- rank.ev.matrix[,grepl("GeneID|Global_Rank_Log2SD",colnames(rank.ev.matrix))]      
    } else if (method=="cv") {
      rank.ev.matrix <- rank.ev.matrix[,grepl("GeneID|Global_Rank_Log2CV",colnames(rank.ev.matrix))]       
    }
  }
  
  rank.ev.matrix <- tibble::column_to_rownames(rank.ev.matrix, var="GeneID")
  
  return(rank.ev.matrix)
}

# Tally low, mid, and high variability genes
tally_ev <- function(ev.matrix, ev.percentile=0.20, output="all") {
  ev.tally <- vector(mode="list", length=2)
  names(ev.tally) <- c("values","summary")
  
  ev.tally$values <- ev.matrix
  
  ev.tally$values$N_Conditions <- rowSums(!is.na(ev.matrix))
  ev.tally$values$LV <- 
    apply(ev.matrix, MARGIN=1, function(r) { sum(r<=ev.percentile, na.rm=TRUE) })
  ev.tally$values$HV <- 
    apply(ev.matrix, MARGIN=1, function(r) { sum(r>=(1.0-ev.percentile), na.rm=TRUE) })
  #ev.tally$values$MV <-
  #  apply(ev.matrix, MARGIN=1, function(r) { 
  #    sum(r>=(0.5-(ev.percentile/2) & r<=(0.5+(ev.percentile/2))), na.rm=TRUE)})
  
  ev.tally$values$MV <- ev.tally$values$N_Conditions-(ev.tally$values$LV + ev.tally$values$HV)
  
  # Ratio
  ev.tally$values$LV_Ratio <- round(ev.tally$values$LV / ev.tally$values$N_Conditions, 3)
  ev.tally$values$HV_Ratio <- round(ev.tally$values$HV / ev.tally$values$N_Conditions, 3)
  ev.tally$values$MV_Ratio <- round(ev.tally$values$MV / ev.tally$values$N_Conditions, 3)
  
  # Summary
  if (output != "values") {
    tally.lv <- dplyr::count(ev.tally$values,LV)
    colnames(tally.lv) <- c("N_Conditions","LV")
    
    tally.hv <- dplyr::count(ev.tally$values,HV)
    colnames(tally.hv) <- c("N_Conditions","HV")
    
    tally.MV <- dplyr::count(ev.tally$values,MV)
    colnames(tally.MV) <- c("N_Conditions","MV")
    
    ev.tally$summary <- 
      Reduce(function(...) merge(..., all=TRUE),list(tally.lv,tally.hv,tally.MV))
    
    ev.tally$summary[is.na(ev.tally$summary)] <- 0
    
    ev.tally$summary$LV <- (ev.tally$summary$LV / nrow(ev.tally$values)) * 100
    ev.tally$summary$HV <- (ev.tally$summary$HV / nrow(ev.tally$values)) * 100
    ev.tally$summary$MV <- (ev.tally$summary$MV / nrow(ev.tally$values)) * 100
  }
  
  if (output=="all") { return(ev.tally) }
  else if (output=="df") { return(ev.tally$values) }
  else if (output=="summary") { return(ev.tally$summary) }
}

# Permute values in each column
permute_by_cols <- function(input.matrix, seed=12345) {
  set.seed(seed)
  
  permuted.matrix <- input.matrix
  for (s in colnames(input.matrix)) { 
    permuted.matrix[,s] <- sample(input.matrix[,s]) 
  }
  
  return(permuted.matrix)
}