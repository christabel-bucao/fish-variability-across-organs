#### Notes ####

#### Libraries ####
library(ggplot2)
library(gplots)
library(viridis)
library(colorBlindness)

#### Sourced functions ####
source("../functions/external/heatmap.3.R")
source("../functions/unit_tests.R")
#source("../functions/annotate_boxplots.R")

#### Plotting functions ####
plot_detected_genes_vs_mapped_reads <- function(metadata, min.genes, min.reads, set.scale=TRUE) {
  Tissue <- as.factor(metadata$Tissue)
  tissue.colors <- set_tissue_colors()[levels(Tissue)]
  
  # set.scale is for aesthetic purposes only
  if (set.scale==TRUE) {
    ggplot(metadata) +
      aes(x=MappedReads/1e6, y=DetectedGenes, color=Tissue, shape=ProjectSeqBatch) +
      geom_point() +
      scale_color_discrete(type=tissue.colors) +
      geom_hline(yintercept=min.genes, linetype="longdash", colour="#2C467A", size=0.5) +
      geom_vline(xintercept=min.reads/1e6, linetype="longdash", colour="#2C467A", size=0.5) +
      labs(x="Uniquely mapped reads (millions)", 
           y="Detected genes (counts > 0)",
           title=paste0(metadata$Species[1], " (n = ", nrow(metadata), ")"),
           subtitle=paste0("x = ", min.reads/1e6, ", y = ", min.genes),
           caption=paste0("produced on ", Sys.time())) +
      scale_x_continuous(breaks=seq(0, 5.5, by=1.0), limits=c(0,5.5)) +
      scale_y_continuous(breaks=seq(0, 22500, by=5000), limits=c(0, 22500)) +
      theme_bw()    
  } else {
    ggplot(metadata) +
      aes(x=MappedReads/1e6, y=DetectedGenes, color=Tissue, shape=ProjectSeqBatch) +
      geom_point() +
      scale_color_discrete(type=tissue.colors) +
      geom_hline(yintercept=min.genes, linetype="longdash", colour="#2C467A", size=0.5) +
      geom_vline(xintercept=min.reads/1e6, linetype="longdash", colour="#2C467A", size=0.5) +
      labs(x="Uniquely mapped reads (millions)", 
           y="Detected genes (counts > 0)",
           title=paste0(metadata$Species[1], " (n = ", nrow(metadata), ")"),
           subtitle=paste0("x = ", min.reads/1e6, ", y = ", min.genes),
           caption=paste0("produced on ", Sys.time())) +
      theme_bw()        
  }

}

plot_pca_by_batch <- function(pca, metadata) {
  if (!samples_match_metadata(rownames(pca$df), metadata)) { 
    return("Error: Samples do not match metadata")
  }
  else {
    Batch <- as.factor(metadata$Batch)
    ggplot(pca$df) +
      aes(x=PC1, y=PC2, color=Batch) +
      viridis::scale_color_viridis(discrete=TRUE) +
      geom_point(size=1.0) +
      xlab(paste0("PC1 (", pca$percent[1], "%", ")")) + 
      ylab(paste0("PC2 (", pca$percent[2], "%", ")")) +
      labs(title="PC1 - PC2", 
           subtitle=paste0("log2 (TMM-CPM) ", "(n = ", ncol(pca$df), ")"),
           caption=paste0("produced on ", Sys.time())) +
      xlim(c(-200,200)) +
      ylim(c(-150,150)) +
      coord_fixed() +
      theme_bw() 
  }  
}

plot_pca_by_mapped_reads <- function(pca, metadata, min.reads) {
  if (!samples_match_metadata(rownames(pca$df), metadata)) { 
    return("Error: Samples do not match metadata")
  }
  else {
    Mapped_reads <- metadata$MappedReads
    ggplot(pca$df) +
      aes(x=PC1, y=PC2, color=(Mapped_reads/1e6)>(min.reads/1e6)) +
      geom_point(size=1.0) +
      xlab(paste0("PC1 (", pca$percent[1], "%", ")")) + 
      ylab(paste0("PC2 (", pca$percent[2], "%", ")")) +
      labs(title="PC1 - PC2", 
           subtitle=paste0("log2 (TMM-CPM) ", "(n = ", ncol(pca$df), ")"),
           color=paste0("Mapped reads > ", min.reads/1e6, "M"),
           caption=paste0("produced on ", Sys.time())) +
      xlim(c(-200,200)) +
      ylim(c(-150,150)) +
      coord_fixed() +
      theme_bw()     
  }
}

plot_pca_by_mapped_reads_gradient <- function(pca, metadata) {
  if (!samples_match_metadata(rownames(pca$df), metadata)) { 
    return("Error: Samples do not match metadata")
  }
  else {
    Mapped_reads <- metadata$MappedReads
    ggplot(pca$df) +
      aes(x=PC1, y=PC2, color=Mapped_reads/1e6) +
      #scale_color_gradientn(colors=gplots::colorpanel(5, "pink","slategray","darkslateblue")) +
      scale_color_gradientn(colors=viridis::viridis(10)) +
      geom_point(size=1.0) +
      xlab(paste0("PC1 (", pca$percent[1], "%", ")")) + 
      ylab(paste0("PC2 (", pca$percent[2], "%", ")")) +
      labs(title="PC1 - PC2", 
           subtitle=paste0("log2 (TMM-CPM) ", "(n = ", ncol(pca$df), ")"),
           color="Mapped reads (millions)",
           caption=paste0("produced on ", Sys.time())) +
      xlim(c(-200,200)) +
      ylim(c(-150,150)) +
      coord_fixed() +
      theme_bw()
  }
}

plot_pca_by_detected_genes <- function(pca, metadata, min.genes) {
  if (!samples_match_metadata(rownames(pca$df), metadata)) { 
    return("Error: Samples do not match metadata")
  }
  else {
    Detected_genes <- metadata$DetectedGenes
    ggplot(pca$df) +
      aes(x=PC1, y=PC2, color=Detected_genes>min.genes) +
      geom_point(size=1.0) +
      xlab(paste0("PC1 (", pca$percent[1], "%", ")")) + 
      ylab(paste0("PC2 (", pca$percent[2], "%", ")")) +
      labs(title="PC1 - PC2", 
           subtitle=paste0("log2 (TMM-CPM) ", "(n = ", ncol(pca$df), ")"),
           color=paste0("Detected genes > ", min.genes),
           caption=paste0("produced on ", Sys.time())) +
      xlim(c(-200,200)) +
      ylim(c(-150,150)) +
      coord_fixed() +
      theme_bw()
  }
}

plot_pca_by_detected_genes_gradient <- function(pca, metadata) {
  if (!samples_match_metadata(rownames(pca$df), metadata)) { 
    return("Error: Samples do not match metadata")
  }
  else {
    Detected_genes <- metadata$DetectedGenes
    ggplot(pca$df) +
      aes(x=PC1, y=PC2, color=Detected_genes) +
      #scale_color_gradientn(colors=gplots::colorpanel(5,"pink","slategray","darkslateblue")) +
      scale_color_gradientn(colors=viridis::viridis(10)) +
      geom_point(size=1.0) +
      xlab(paste0("PC1 (", pca$percent[1], "%", ")")) + 
      ylab(paste0("PC2 (", pca$percent[2], "%", ")")) +
      labs(title="PC1 - PC2", 
           subtitle=paste0("log2 (TMM-CPM) ", "(n = ", ncol(pca$df), ")"),
           color="Detected genes",
           caption=paste0("produced on ", Sys.time())) +
      xlim(c(-200,200)) +
      ylim(c(-150,150)) +
      coord_fixed() +
      theme_bw()    
  }
}

plot_pca_by_condition <- function(pca, metadata) {
  if (!samples_match_metadata(rownames(pca$df), metadata)) { 
    return("Error: Samples do not match metadata")
  }
  else {
    Tissue <- as.factor(metadata$Tissue)
    tissue.colors <- set_tissue_colors()[levels(Tissue)]
    
    ## Check if sex is annotated
    if ("Sex" %in% colnames(metadata)) {
      Sex <- as.factor(metadata$Sex)
      
      ggplot(pca$df) +
        aes(x=PC1, y=PC2, color=Tissue, shape=Sex) +
        geom_point(size=1.0) +
        scale_color_discrete(type=tissue.colors) +
        xlab(paste0("PC1 (", pca$percent[1], "%", ")")) + 
        ylab(paste0("PC2 (", pca$percent[2], "%", ")")) +
        labs(title="PC1 - PC2", 
             subtitle=paste0("log2 (TMM-CPM) ", "(n = ", ncol(pca$df), ")"),
             caption=paste0("produced on ", Sys.time())) +
        xlim(c(-200,200)) +
        ylim(c(-150,150)) +
        coord_fixed() +
        theme_bw() 
    }
    else {
      ggplot(pca$df) +
        aes(x=PC1, y=PC2, color=Tissue) +
        geom_point(size=1.0) +
        scale_color_discrete(type=tissue.colors) +
        xlab(paste0("PC1 (", pca$percent[1], "%", ")")) + 
        ylab(paste0("PC2 (", pca$percent[2], "%", ")")) +
        labs(title="PC1 - PC2", 
             subtitle=paste0("log2 (TMM-CPM) ", "(n = ", ncol(pca$df), ")"),
             caption=paste0("produced on ", Sys.time())) +
        xlim(c(-200,200)) +
        ylim(c(-150,150)) +
        coord_fixed() +
        theme_bw() 
    }    
  }
}

plot_pca_by_condition_pc34 <- function(pca, metadata) {
  if (!samples_match_metadata(rownames(pca$df), metadata)) { 
    return("Error: Samples do not match metadata")
  }
  else {
    Tissue <- as.factor(metadata$Tissue)
    tissue.colors <- set_tissue_colors()[levels(Tissue)]
    
    ## Check if sex is annotated
    if ("Sex" %in% colnames(metadata)) {
      Sex <- as.factor(metadata$Sex)
      
      ggplot(pca$df) +
        aes(x=PC3, y=PC4, color=Tissue, shape=Sex) +
        geom_point(size=1.0) +
        scale_color_discrete(type=tissue.colors) +
        xlab(paste0("PC3 (", pca$percent[3], "%", ")")) + 
        ylab(paste0("PC4 (", pca$percent[4], "%", ")")) +
        labs(title="PC3 - PC4", 
             subtitle=paste0("log2 (TMM-CPM) ", "(n = ", ncol(pca$df), ")"),
             caption=paste0("produced on ", Sys.time())) +
        xlim(c(-200,200)) +
        ylim(c(-150,150)) +
        coord_fixed() +
        theme_bw() 
    }
    else {
      ggplot(pca$df) +
        aes(x=PC3, y=PC4, color=Tissue) +
        geom_point(size=1.0) +
        scale_color_discrete(type=tissue.colors) +
        xlab(paste0("PC3 (", pca$percent[3], "%", ")")) + 
        ylab(paste0("PC4 (", pca$percent[4], "%", ")")) +
        labs(title="PC3 - PC4", 
             subtitle=paste0("log2 (TMM-CPM) ", "(n = ", ncol(pca$df), ")"),
             caption=paste0("produced on ", Sys.time())) +
        xlim(c(-200,200)) +
        ylim(c(-150,150)) +
        coord_fixed() +
        theme_bw() 
    }    
  }
}

boxplot_detected_genes <- function(metadata) {
  Tissue <- as.factor(metadata$Tissue)
  tissue.colors <- set_tissue_colors()[levels(Tissue)]
  
  Species <- levels(as.factor(metadata$Species))
  
  custom_theme <- theme(plot.title=element_text(size=24, face="italic"),
                        axis.title=element_text(size=22),
                        axis.text=element_text(size=20),
                        axis.text.x=element_text(size=20, angle=30),
                        strip.text=element_text(size=20),
                        panel.spacing.x=unit(1.1,"lines"),
                        legend.position="none")
  
  min.genes <- min(metadata$DetectedGenes)
  
  annotate_n <- function(x) {
    return(c(y = min.genes*0.95, label = length(x))) 
  }
  
  if ("Sex" %in% colnames(metadata)) {
    ggplot(metadata,
           aes(x=Tissue, y=DetectedGenes, fill=Tissue)) +
      stat_boxplot(geom='errorbar') +
      geom_boxplot(alpha=0.6, outlier.color=NA) +
      geom_jitter(size=2, alpha=0.9) +
      scale_fill_manual(values=tissue.colors) +
      stat_summary(fun.data=annotate_n, geom="text", fun=median,
                   position=position_dodge(width=0.75), size=8) +
      facet_wrap(vars(Sex), ncol=2) +
      labs(title=gsub("_"," ", Species),
           x="Organ",
           y="Detected genes (counts > 0)",
           fill="Organ") +
      theme_bw() +
      custom_theme
    
  } else {
    ggplot(metadata,
           aes(x=Tissue, y=DetectedGenes, fill=Tissue)) +
      stat_boxplot(geom='errorbar') +
      geom_boxplot(alpha=0.6, outlier.color=NA) +
      geom_jitter(size=2, alpha=0.9) +
      scale_fill_manual(values=tissue.colors) +
      stat_summary(fun.data=annotate_n, geom="text", fun=median,
                   position=position_dodge(width=0.75), size=8) +
      labs(title=gsub("_"," ", Species),
           x="Organ",
           y="Detected genes (counts > 0)",
           fill="Organ") +
      theme_bw() +
      custom_theme    
  }
  
}

plot_correlation_heatmap <- function(cor.matrix, metadata, clust=FALSE) {
  if (!samples_match_metadata(colnames(cor.matrix), metadata)) { 
    return("Error: Samples do not match metadata")
  }  
  else {
    Tissue <- as.factor(metadata$Tissue)
    tissue.colors <- set_tissue_colors()[levels(Tissue)]
    
    organ <- 
      vapply(colnames(cor.matrix), FUN.VALUE=character(1), function(id) {
        tissue.colors[metadata$Tissue[metadata$SampleName==id]]
      })
    
    ## Check if sex is annotated
    if ("Sex" %in% colnames(metadata)) {
      #sex.colors <- c("F"="#EF8A62","M"="#67A9CF")
      sex.colors <- c("F"="#BB00BB","M"="#00BB00")
      sex <- 
        vapply(colnames(cor.matrix), FUN.VALUE=character(1), function(id) {
          sex.colors[metadata$Sex[metadata$SampleName==id]]
        })
      
      col.labels <- cbind(organ, sex)
      row.labels <- t(col.labels)
      row.labels <- rbind(row.labels[2,], row.labels[1,])
    }
    else { 
      col.labels <- cbind(organ)
      row.labels <- rbind(t(col.labels))
    }
    
    if (clust==FALSE) {
      heatmap.3(cor.matrix, 
                #scale = "row"
                Rowv=FALSE, 
                Colv=FALSE,
                dendrogram="none",
                col=gplots::colorpanel(100,"#F0F0F0","#BDBDBD","#636363"),
                breaks=seq(0,1,length.out=101),
                labCol=NA,
                labRow=NA,
                trace="none",
                density.info="none",
                cexRow=1, cexCol=1,
                ColSideColors=col.labels,
                RowSideColors=row.labels,
                lwid=c(1.0,4.0),
                lhei=c(1.0,4.0),
                main=paste0("log2 (TMM-CPM) Pearson correlation (n = ", nrow(metadata), ")"))
    }
    else if (clust==TRUE) {
      clust <- hclust(as.dist(1-cor.matrix), method="complete") 
      
      heatmap.3(cor.matrix, 
                #scale = "row"
                Rowv=as.dendrogram(clust), 
                Colv=as.dendrogram(clust),
                dendrogram="column",
                col=gplots::colorpanel(100,"#F0F0F0","#BDBDBD","#636363"),
                breaks=seq(0,1,length.out=101),
                labCol=NA,
                labRow=NA,
                trace="none",
                density.info="none",
                cexRow=1, cexCol=1,
                ColSideColors=col.labels,
                RowSideColors=row.labels,
                lwid=c(1.0,4.0),
                lhei=c(1.0,4.0),
                main=paste0("log2 (TMM-CPM) Pearson correlation (n = ", nrow(metadata), ")"))
    }
    
    legend("left", title="Organs", legend=levels(Tissue), 
           fill=tissue.colors, 
           cex=1.0, bty="n")
    if ("Sex" %in% colnames(metadata)) {
      legend("bottomleft", title="Sex", legend=names(sex.colors), 
             fill=sex.colors, 
             cex=1.0, bty="n")
    }    
  }
}