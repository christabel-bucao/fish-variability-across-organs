#### Notes ####

#### Libraries ####
library(GOstats)

#### Sourced functions ####

#### Functions ####

# Retrieve low, mid, and high variability genes
retrieve_variability_gene_lists <- function(jack.summary, ev.percentile=0.20, col="Mean_Local_Rank_Log2CV") {
  bgd.genes <- jack.summary$GeneID
  lv.genes <- jack.summary$GeneID[jack.summary[[col]]<=ev.percentile]
  hv.genes <- jack.summary$GeneID[jack.summary[[col]]>=(1.0-ev.percentile)]
  mv.genes <- 
    jack.summary$GeneID[jack.summary[[col]]>=(0.5-(ev.percentile/2)) &
                          jack.summary[[col]]<=(0.5+(ev.percentile/2))]
  
  ev.genes <- list("LV"=lv.genes, "MV"=mv.genes,"HV"=hv.genes, "Bgd"=bgd.genes)
  return(ev.genes)
}

# Combine GO results across conditions
combine_go_results <- function(go.summary, sex=TRUE) {
  go.combined <- 
    lapply(c("LV"="LV","MV"="MV","HV"="HV"), function(ev) {
      lapply(c("over"="over","under"="under"), function(test) { NULL })
    })
  
  if (sex==TRUE) {
    for (t in levels(Tissue)) {
      for (s in levels(Sex)) {
        for (ev in c("LV","MV","HV")) {
          for (test in c("over","under")) {
            if (class(go.summary[[t]][[s]][[ev]][[test]])=="data.frame") {
              go.combined[[ev]][[test]] <-
                rbind(go.combined[[ev]][[test]], go.summary[[t]][[s]][[ev]][[test]])              
            }
          }
        }
      }
    }  
  } else {
    for (t in levels(Tissue)) {
      for (ev in c("LV","MV","HV")) {
        for (test in c("over","under")) {
          if (class(go.summary[[t]][[ev]][[test]])=="data.frame") {
            go.combined[[ev]][[test]] <-
              rbind(go.combined[[ev]][[test]], go.summary[[t]][[ev]][[test]])
          }
        }
      }
    }   
  }
  
  go.combined <- 
    lapply(c("LV"="LV","MV"="MV","HV"="HV"), function(ev) {
      lapply(c("over"="over","under"="under"), function(test) {
        dplyr::group_by(go.combined[[ev]][[test]], GOBPID, Term) %>%
          dplyr::summarise(
            MeanPvalue=mean(Pvalue),
            NConditions=n()
          )
      })
    })
  
  return(go.combined)
}

# Export GO results for external plotting
export_go_results <- function(go.combined, condition.cutoff, output.dir, prefix) {
  terms.dir <- file.path(output.dir, "terms")
  if (!dir.exists(terms.dir)) { dir.create(terms.dir, recursive=TRUE) }
  
  go.export <- go.combined
  
  for (ev in c("LV","MV","HV")) {
    for (test in c("over","under")) {
      go.export[[ev]][[test]] <- 
        go.export[[ev]][[test]][go.export[[ev]][[test]]$NConditions>=condition.cutoff, 
                                c("GOBPID","MeanPvalue", "NConditions")]
    }
  }
  
  for (ev in c("LV","MV","HV")) {
    for (test in c("over","under")) {
      write.table(
        go.export[[ev]][[test]], 
        file.path(terms.dir, paste(paste(prefix,tolower(ev),test,sep="_"),"txt",sep=".")),
        quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t"
      )
    }
  }  
}

run_gostats_pipeline <- function(foreground, background, annotation.db, pvalue.cutoff, go.ontology="BP") {
  # Initialize
  gostats.analysis <- list()
  
  gostats.analysis$genes <- list()
  gostats.analysis$genes$foreground <- foreground
  gostats.analysis$genes$background <- background
  
  gostats.analysis$over <- list()
  gostats.analysis$under <- list()
  
  # Set parameters
  # Enrichment
  print("Setting parameters")
  gostats.analysis$over$params <-
    new("GOHyperGParams",
        geneIds=foreground,
        universeGeneIds=background,
        annotation=annotation.db,
        ontology=go.ontology,
        pvalueCutoff=pvalue.cutoff,
        conditional=TRUE,
        testDirection="over")
  
  # Depletion
  gostats.analysis$under$params <-
    new("GOHyperGParams",
        geneIds=foreground,
        universeGeneIds=background,
        annotation=annotation.db,
        ontology=go.ontology,
        pvalueCutoff=pvalue.cutoff,
        conditional=TRUE,
        testDirection="under")
  
  # Run conditional hypergeometric test
  print("Running conditional hypergeometric test")
  gostats.analysis$over$results <- hyperGTest(gostats.analysis$over$params)
  gostats.analysis$under$results <- hyperGTest(gostats.analysis$under$params)
  
  # Summarize results
  print("Summarizing results")
  gostats.analysis$over$summary <- summary(gostats.analysis$over$results)
  gostats.analysis$under$summary <- summary(gostats.analysis$under$results)
  
  return(gostats.analysis)
}