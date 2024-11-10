#### Notes ####
# Functions to filter samples for quality control

#### Libraries ####
library(edgeR)
library(matrixStats)

#### Sourced functions ####
source("../functions/unit_tests.R")
source("../functions/helper_functions.R")

#### Functions ####

# Remove barcodes not included in metadata and map to sample names
map_barcodes_to_sample_names <- function(counts.df, metadata) {
  filtered.counts.df <- counts.df[,c("Gene_id", metadata$Barcode)]
  
  sample.names <-
    vapply(colnames(filtered.counts.df)[-1], FUN.VALUE=character(1), function(col) {
      metadata$SampleName[metadata$Barcode==col]
    })
  colnames(filtered.counts.df) <- c("GeneID", unlist(sample.names))
  
  return(filtered.counts.df)
}

# Count number of uniquely mapped reads and detected genes per sample
summarize_sequencing_results <- function(counts.matrix) {
  ## Reads per sample
  mapped.reads <- apply(counts.matrix, MARGIN=2, sum)
  
  ## Genes with count > 0
  detected.genes <- apply(counts.matrix, MARGIN=2, function(col) { length(col[col>0]) })
  
  summary <- data.frame("SampleName"=colnames(counts.matrix), 
                        "MappedReads"=mapped.reads,
                        "DetectedGenes"=detected.genes)
  
  return(summary)
}

# Remove genes with low abundance
filter_genes_by_cpm <- function(counts.matrix, metadata, min.cpm, tmm=FALSE) {
  if (tmm=="TRUE") {
    dge <- edgeR::DGEList(counts=counts.matrix, samples=metadata)
    tmm.counts <- edgeR::calcNormFactors(dge, method="TMM")
    cpm <- edgeR::cpm(tmm.counts)    
  } else {
    cpm <- edgeR::cpm(counts.matrix)
  }
  
  filtered.counts <- 
    counts.matrix[rowMeans(cpm)>min.cpm & matrixStats::rowMedians(cpm)>min.cpm,]
  return(filtered.counts)
}

# Remove genes with low abundance in all tissues
filter_genes_by_cpm_per_condition <- function(counts.matrix, metadata, min.cpm) {
  if (!samples_match_metadata(colnames(counts.matrix), metadata)) { 
    return("Error: Samples do not match metadata")
  } 
  else {
    cpm <- edgeR::cpm(counts.matrix)
    
    if ("Sex" %in% colnames(metadata)) { 
      metadata <- split_gonads_df(metadata)
      colnames(cpm) <- metadata$SampleName
      }
    
    Tissue <- as.factor(metadata$Tissue)
    mean.cpm <- matrix(nrow=nrow(cpm), ncol=length(levels(Tissue)),
                       dimnames=list(rownames(cpm), levels(Tissue)))
    median.cpm <- matrix(nrow=nrow(cpm), ncol=length(levels(Tissue)),
                       dimnames=list(rownames(cpm), levels(Tissue)))
    
    for (t in levels(Tissue)) {
      cpm.per.tissue <- cpm[,grepl(t,colnames(cpm)), drop=FALSE]
      mean.cpm[,t] <- round(rowMeans(cpm.per.tissue),4)
      median.cpm[,t] <- round(matrixStats::rowMedians(cpm.per.tissue, useNames=TRUE),4)
    }
    
    ## Keep genes with mean cpm > min.cpm & median cpm > min.cpm in at least one condition
    keep.mean <- 
      apply(mean.cpm, MARGIN=1, function(row) { any(row>min.cpm) })
    keep.median <-
      apply(median.cpm, MARGIN=1, function(row) { any(row>min.cpm) })
    filtered.counts <- counts.matrix[keep.mean & keep.median,]
    
    return(filtered.counts) 
  }
}

# Transform counts to TMM-CPM or log2 TMM-CPM
transform_counts_to_tmm_cpm <- function(counts.matrix, metadata, log2=FALSE, shift.min=FALSE) {
  dge <- edgeR::DGEList(counts=counts.matrix, samples=metadata)
  tmm.counts <- edgeR::calcNormFactors(dge, method="TMM")
  
  if (log2==FALSE) { 
    tmm.cpm <- edgeR::cpm(tmm.counts, log=FALSE)
    return(tmm.cpm)
  }
  else {
    log2.tmm.cpm <- edgeR::cpm(tmm.counts, log=TRUE)
    if (shift.min==FALSE) {
      return(log2.tmm.cpm)
    }
    else { ## Shift log2 (TMM-CPM) values by minimum
      expr.min <- min(log2.tmm.cpm) ## This is actually 0
      log2.tmm.cpm.shift <- round(log2.tmm.cpm - expr.min, 8)
      return(log2.tmm.cpm.shift)
    }
  }
}

# Keep only reduced expression matrix without any zeroes
keep_full_nonzero_expression_matrix <- function(expr.matrix) {
  row.nonzeroes <-
    apply(expr.matrix, MARGIN=1, function(r) {
      sum(r>0.0)
    })
  
  nonzero.expr.matrix <- expr.matrix[which(row.nonzeroes==ncol(expr.matrix)),]
  
  return(nonzero.expr.matrix)
}

# Remove samples with less than minimum detected genes OR less than minimum mapped reads
filter_samples_by_detected_genes_mapped_reads <- function(counts.matrix, metadata, min.genes, min.reads) {
  filtered.metadata <- 
    metadata[!(metadata$DetectedGenes<min.genes | metadata$MappedReads<min.reads),]
  
  filtered.counts <-
    counts.matrix[,colnames(counts.matrix) %in% filtered.metadata$SampleName]
  
  filtered.data <- list("metadata"=filtered.metadata, "counts"=filtered.counts)
  return(filtered.data)
}

# Run principal component analysis
run_pca <- function(expr.matrix) {
  pca <- list()
  pca$results <- prcomp(t(expr.matrix), scale.=F, retx=T) ### PCA function
  pca$variance <- (pca$results$sdev)^2
  pca$percent <- round(pca$variance/sum(pca$variance)*100, 1)  ## Percentage of variance
  pca$df <- as.data.frame(pca$results$x)  
  
  return(pca)
}

# Note samples in which the average within-tissue correlation is lower than the maximum between-tissue correlation
tag_samples_by_correlation <- function(cor.matrix, metadata, round=2) {
  Tissue <- as.factor(metadata$Tissue)
  
  mean.tissue.cor <- matrix(nrow=nrow(cor.matrix), ncol=length(levels(Tissue)),
                     dimnames=list(rownames(cor.matrix), levels(Tissue)))
  
  samples <- NULL
    
  for (id in rownames(cor.matrix)) {
    id.ind <- metadata$IndID[metadata$SampleName==id]
      
    # Reference tissue
    tissue.ref <- as.character(metadata$Tissue[metadata$SampleName==id])
      
    for (t in levels(Tissue)) {
      if (t==tissue.ref) {
        # Sample IDs for the reference tissue, excluding sample from the same individual
        tissue.ids <- 
          subset(rownames(cor.matrix), 
                 rownames(cor.matrix) %in% metadata$SampleName[metadata$Tissue==t]
                  & !(rownames(cor.matrix) %in% metadata$SampleName[metadata$IndID==id.ind]))
      }
      else {
        # Sample IDs for a specified tissue
        tissue.ids <- 
          subset(rownames(cor.matrix), 
                 rownames(cor.matrix) %in% metadata$SampleName[metadata$Tissue==t])
      }
      
      mean.tissue.cor[id,t] <- mean(cor.matrix[id,tissue.ids])
    }
      
    # Note samples in which the average within-tissue correlation 
    # is lower than the maximum between-tissue correlation
    if (round(mean.tissue.cor[id,tissue.ref], round) 
        != round(max(mean.tissue.cor[id,]), round)) { 
      samples <- c(samples, id) 
      } 
    }
    
    cor.analysis <- list("cor"=mean.tissue.cor, "samples"=samples)
    return(cor.analysis)
}

# Iterative removal of samples that are not well-correlated with annotated tissue
iterate_filter_samples_by_correlation <- function(counts.matrix, metadata, min.cpm, round, shift.min=TRUE) {
  if (!samples_match_metadata(colnames(counts.matrix), metadata)) { 
    return("Error: Samples do not match metadata")
  } 
  else {
    if ("Sex" %in% colnames(metadata)) {
      orig.colnames <- colnames(counts.matrix)
      metadata <- split_gonads_df(metadata)
      colnames(counts.matrix) <- metadata$SampleName
    }    
    
    iter=TRUE
    while (iter==TRUE) {
      
      print("Filtering and normalizing counts...")
      counts.gf <- 
        filter_genes_by_cpm_per_condition(counts.matrix, metadata, min.cpm)
      log2.tmm.cpm.gf <- transform_counts_to_tmm_cpm(counts.gf, metadata, log2=TRUE, shift.min)
      
      cor.matrix <- cor(log2.tmm.cpm.gf, method="pearson")
      cor.test <- 
        tag_samples_by_correlation(cor.matrix, metadata, round)
      
      if (length(cor.test$samples)>0) {
        print("Removing samples:")
        print(cor.test$samples)
        
        metadata <- 
          metadata[!(metadata$SampleName %in% cor.test$samples),]
        counts.matrix <- 
          counts.matrix[,colnames(counts.matrix) %in% metadata$SampleName]
      }
      else { 
        print("No samples to remove")
        iter=FALSE 
      }
    }
    
    if ("Sex" %in% colnames(metadata)) { 
      metadata <- join_gonads_df(metadata)
      colnames(counts.matrix) <- metadata$SampleName
      colnames(counts.gf) <- metadata$SampleName
      colnames(log2.tmm.cpm.gf) <- metadata$SampleName
      colnames(cor.matrix) <- metadata$SampleName
      rownames(cor.matrix) <- metadata$SampleName
      conditions <- dplyr::count(metadata, Tissue, Sex)
    }
    else {
      conditions <- dplyr::count(metadata, Tissue)
    }
    
    filtered.samples <- 
      list("tissue.cor"=cor.test$cor, "conditions"=conditions,
           "metadata"=metadata, "cor.matrix"=as.matrix(cor.matrix),
           "counts"=as.matrix(counts.matrix), "counts.gf"=as.matrix(counts.gf), 
           "log2.tmm.cpm.gf"=as.matrix(log2.tmm.cpm.gf))
    
    return(filtered.samples)
  }
}

# Remove tissue-sex conditions with few replicates
filter_conditions_by_min_replicates <- function(counts.matrix, metadata, min.reps) {
  if ("Sex" %in% colnames(metadata)) {
    conditions <- dplyr::count(metadata, Tissue, Sex)
    conditions$regex <- paste0(conditions$Sex,"[0-9][0-9]_", conditions$Tissue)
    conditions.to.remove <- conditions$regex[conditions$n<min.reps]
    
    if (length(conditions.to.remove)==0) {
      print("No conditions to remove")
      filtered.data <- list("conditions"=conditions,
                            "metadata"=metadata, 
                            "counts"=counts.matrix)
      return(filtered.data)
    }
    else {
      print("Removing conditions:")
      print(conditions.to.remove)
      
      filtered.metadata <- metadata
      for (rgx in conditions.to.remove) {
        filtered.metadata <- filtered.metadata[!(grepl(rgx, filtered.metadata$SampleName)),]
      }
      filtered.conditions <- dplyr::count(filtered.metadata, Tissue, Sex)      
    }
  }
  else { 
    conditions <- dplyr::count(metadata, Tissue)
    conditions.to.remove <- conditions$Tissue[conditions$n<min.reps]
    print(conditions.to.remove)
    
    if (length(conditions.to.remove)==0) {
      print("No conditions to remove")
      filtered.data <- list("conditions"=conditions,
                            "metadata"=metadata, 
                            "counts"=counts.matrix)
      return(filtered.data)
    }
    else {
      print("Removing conditions:")
      print(conditions.to.remove)
      
      filtered.metadata <- metadata[!(metadata$Tissue %in% conditions.to.remove),]
      filtered.conditions <- dplyr::count(filtered.metadata, Tissue)      
    }
  }
  filtered.counts <- counts.matrix[,filtered.metadata$SampleName]
  filtered.data <- list("conditions"=filtered.conditions,
                        "metadata"=filtered.metadata, 
                        "counts"=filtered.counts)
  return(filtered.data)
}