#### Notes on revisions ####


#### Set species ####
species.name <- c("Lepisosteus oculatus","Esox lucius","Danio rerio")
names(species.name) <- c("LOC","ELU","DRE")


#### Step 1: (Skipped) Filter samples ####

#### Step 2: Estimate gene expression variability - With jackknife ####

##### Set input and parameters #####
step02.variability.jackknife <- list()
step02.variability.jackknife$input <- "analysis/02_Variability_Jackknife.Rmd"
step02.variability.jackknife$output.dir <- "../results/02_Variability_Jackknife"

# Create results directory
if (!dir.exists(step02.variability.jackknife$output.dir)) {
  dir.create(step02.variability.jackknife$output.dir, recursive=TRUE, showWarnings=FALSE)
}

# List parameters
step02.variability.jackknife$params <- list(
  list(species="LOC", species.name=species.name[["LOC"]], min.percentile=0.00, max.percentile=0.95, win.size=100, qc1.only=TRUE, all.nonzero.matrix=FALSE, combat=FALSE, run.example=FALSE),
  list(species="ELU", species.name=species.name[["ELU"]], min.percentile=0.00, max.percentile=0.95, win.size=100, qc1.only=TRUE, all.nonzero.matrix=FALSE, combat=FALSE, run.example=FALSE),
  list(species="DRE", species.name=species.name[["DRE"]], min.percentile=0.00, max.percentile=0.95, win.size=100, qc1.only=TRUE, all.nonzero.matrix=FALSE, combat=FALSE, run.example=FALSE),
  list(species="LOC", species.name=species.name[["LOC"]], min.percentile=0.00, max.percentile=0.95, win.size=100, qc1.only=TRUE, all.nonzero.matrix=FALSE, combat=TRUE, run.example=FALSE),
  list(species="ELU", species.name=species.name[["ELU"]], min.percentile=0.00, max.percentile=0.95, win.size=100, qc1.only=TRUE, all.nonzero.matrix=FALSE, combat=TRUE, run.example=FALSE),
  list(species="DRE", species.name=species.name[["DRE"]], min.percentile=0.00, max.percentile=0.95, win.size=100, qc1.only=TRUE, all.nonzero.matrix=FALSE, combat=TRUE, run.example=FALSE)
)
# If qc1.only=TRUE, we use samples that have been filtered by (1) sequencing quality but not (2) within-sample correlation

##### Run analysis script #####
for (p in (1:length(step02.variability.jackknife$params))) {
  print(paste("Species:", step02.variability.jackknife$params[[p]]$species))
  print(paste("Run jackknife?:", !step02.variability.jackknife$params[[p]]$run.example))
  
  if (step02.variability.jackknife$params[[p]]$combat==FALSE) {
    output.dir <- file.path(step02.variability.jackknife$output.dir, "no_combat", "qc1_only")
  } else {
    output.dir <- file.path(step02.variability.jackknife$output.dir, "combat", "qc1_only")
  }
  
  if (!dir.exists(output.dir)) { dir.create(output.dir, recursive=TRUE, showWarnings=FALSE) }
  
  if (step02.variability.jackknife$params[[p]]$all.nonzero.matrix==FALSE) {
    if (step02.variability.jackknife$params[[p]]$run.example==TRUE) {
      rmarkdown::render(
        input=step02.variability.jackknife$input,
        params=step02.variability.jackknife$params[[p]],
        output_file=file.path(
          "..", output.dir,
          paste(step02.variability.jackknife$params[[p]]$species, "example.html", sep="_")))
    } else {
      rmarkdown::render(
        input=step02.variability.jackknife$input,
        params=step02.variability.jackknife$params[[p]],
        output_file=file.path(
          "..", output.dir,
          paste(step02.variability.jackknife$params[[p]]$species, ".html", sep="")))
    }       
    
  } else {
    if (step02.variability.jackknife$params[[p]]$run.example==TRUE) {
      rmarkdown::render(
        input=step02.variability.jackknife$input,
        params=step02.variability.jackknife$params[[p]],
        output_file=file.path(
          "..", output.dir,
          paste(step02.variability.jackknife$params[[p]]$species, "nonzero", "example.html", sep="_")))
    } else {
      rmarkdown::render(
        input=step02.variability.jackknife$input,
        params=step02.variability.jackknife$params[[p]],
        output_file=file.path(
          "..", output.dir,
          paste(step02.variability.jackknife$params[[p]]$species, "nonzero.html", sep="_")))
    }       
  }
}

# Clear workspace
keep <- ls()[grepl("step[0-9]*|species.name",ls())]
rm(list=setdiff(ls(), keep))
gc()


#### Step 2b: Estimate gene expression variability - Without jackknife ####

##### Set input and parameters #####
step02b.variability.nojack <- list()
step02b.variability.nojack$input <- "analysis/02b_Variability_NoJack.Rmd"
step02b.variability.nojack$output.dir <- "../results/02b_Variability_NoJack"

# Create results directory
if (!dir.exists(step02b.variability.nojack$output.dir)) {
  dir.create(step02b.variability.nojack$output.dir, recursive=TRUE, showWarnings=FALSE)
}

# List parameters
step02b.variability.nojack$params <- list(
  list(species="LOC", species.name=species.name[["LOC"]], min.percentile=0.00, max.percentile=0.95, win.size=100, qc1.only=TRUE, all.nonzero.matrix=FALSE, combat=FALSE, run.example=FALSE),
  list(species="ELU", species.name=species.name[["ELU"]], min.percentile=0.00, max.percentile=0.95, win.size=100, qc1.only=TRUE, all.nonzero.matrix=FALSE, combat=FALSE, run.example=FALSE),
  list(species="DRE", species.name=species.name[["DRE"]], min.percentile=0.00, max.percentile=0.95, win.size=100, qc1.only=TRUE, all.nonzero.matrix=FALSE, combat=FALSE, run.example=FALSE),
  list(species="LOC", species.name=species.name[["LOC"]], min.percentile=0.00, max.percentile=0.95, win.size=100, qc1.only=TRUE, all.nonzero.matrix=FALSE, combat=TRUE, run.example=FALSE),
  list(species="ELU", species.name=species.name[["ELU"]], min.percentile=0.00, max.percentile=0.95, win.size=100, qc1.only=TRUE, all.nonzero.matrix=FALSE, combat=TRUE, run.example=FALSE),
  list(species="DRE", species.name=species.name[["DRE"]], min.percentile=0.00, max.percentile=0.95, win.size=100, qc1.only=TRUE, all.nonzero.matrix=FALSE, combat=TRUE, run.example=FALSE)
)
# If qc1.only=TRUE, we use samples that have been filtered by (1) sequencing quality but not (2) within-sample correlation

##### Run analysis script #####
for (p in (1:length(step02b.variability.nojack$params))) {
  print(paste("Species:", step02b.variability.nojack$params[[p]]$species))
  
  if (step02b.variability.nojack$params[[p]]$combat==FALSE) {
    output.dir <- file.path(step02b.variability.nojack$output.dir, "no_combat", "qc1_only")
  } else {
    output.dir <- file.path(step02b.variability.nojack$output.dir, "combat", "qc1_only")
  }
  
  if (!dir.exists(output.dir)) { dir.create(output.dir, recursive=TRUE, showWarnings=FALSE) }
  
  if (step02b.variability.nojack$params[[p]]$all.nonzero.matrix==FALSE) {
    if (step02b.variability.nojack$params[[p]]$run.example==TRUE) {
      rmarkdown::render(
        input=step02b.variability.nojack$input,
        params=step02b.variability.nojack$params[[p]],
        output_file=file.path(
          "..", output.dir,
          paste(step02b.variability.nojack$params[[p]]$species, "example.html", sep="_")))
    } else {
      rmarkdown::render(
        input=step02b.variability.nojack$input,
        params=step02b.variability.nojack$params[[p]],
        output_file=file.path(
          "..", output.dir,
          paste(step02b.variability.nojack$params[[p]]$species, ".html", sep="")))
    }       
    
  } else {
    if (step02b.variability.nojack$params[[p]]$run.example==TRUE) {
      rmarkdown::render(
        input=step02b.variability.nojack$input,
        params=step02b.variability.nojack$params[[p]],
        output_file=file.path(
          "..", output.dir,
          paste(step02b.variability.nojack$params[[p]]$species, "nonzero", "example.html", sep="_")))
    } else {
      rmarkdown::render(
        input=step02b.variability.nojack$input,
        params=step02b.variability.nojack$params[[p]],
        output_file=file.path(
          "..", output.dir,
          paste(step02b.variability.nojack$params[[p]]$species, "nonzero.html", sep="_")))
    }       
  }
}

# Clear workspace
keep <- ls()[grepl("step[0-9]*|species.name",ls())]
rm(list=setdiff(ls(), keep))
gc()


#### Step 3: Plot expression variability metrics ####

##### Set input and parameters #####
step03.variability.plots <- list()
step03.variability.plots$input <- "analysis/03_Variability_Plots.Rmd"
step03.variability.plots$output.dir <- "../results/03_Variability_Plots"

# Create results directory
if (!dir.exists(step03.variability.plots$output.dir)) {
  dir.create(step03.variability.plots$output.dir, recursive=TRUE, showWarnings=FALSE)
}

# List parameters
step03.variability.plots$params <- list(
  list(species="LOC", species.name=species.name[["LOC"]], min.percentile=0.00, max.percentile=0.95, qc1.only=TRUE, all.nonzero.matrix=FALSE, combat=FALSE),
  list(species="ELU", species.name=species.name[["ELU"]], min.percentile=0.00, max.percentile=0.95, qc1.only=TRUE, all.nonzero.matrix=FALSE, combat=FALSE),
  list(species="DRE", species.name=species.name[["DRE"]], min.percentile=0.00, max.percentile=0.95, qc1.only=TRUE, all.nonzero.matrix=FALSE, combat=FALSE),
  list(species="LOC", species.name=species.name[["LOC"]], min.percentile=0.00, max.percentile=0.95, qc1.only=TRUE, all.nonzero.matrix=FALSE, combat=TRUE),
  list(species="ELU", species.name=species.name[["ELU"]], min.percentile=0.00, max.percentile=0.95, qc1.only=TRUE, all.nonzero.matrix=FALSE, combat=TRUE),
  list(species="DRE", species.name=species.name[["DRE"]], min.percentile=0.00, max.percentile=0.95, qc1.only=TRUE, all.nonzero.matrix=FALSE, combat=TRUE)
)

##### Run analysis script #####
for (p in (1:length(step03.variability.plots$params))) {
  print(paste("Species:", step03.variability.plots$params[[p]]$species))
  
  if (step03.variability.plots$params[[p]]$combat==FALSE) {
    output.dir <- file.path(step03.variability.plots$output.dir, "no_combat", "qc1_only")
  } else {
    output.dir <- file.path(step03.variability.plots$output.dir, "combat", "qc1_only")
  }
  
  if (!dir.exists(output.dir)) { dir.create(output.dir, recursive=TRUE, showWarnings=FALSE) }
  
  if (step03.variability.plots$params[[p]]$all.nonzero.matrix==FALSE) {
    rmarkdown::render(
      input=step03.variability.plots$input,
      params=step03.variability.plots$params[[p]],
      output_file=file.path(
        "..", output.dir,
        paste(step03.variability.plots$params[[p]]$species, ".html", sep="")))      
    
    
  } else {
    rmarkdown::render(
      input=step03.variability.plots$input,
      params=step03.variability.plots$params[[p]],
      output_file=file.path(
        "..", output.dir,
        paste(step03.variability.plots$params[[p]]$species, "nonzero.html", sep="_")))
    
  }
}

# Clear workspace
keep <- ls()[grepl("step[0-9]*|species.name",ls())]
rm(list=setdiff(ls(), keep))
gc()


#### Step 4: Check for correlation across conditions ####

##### Set input and parameters #####
step04.correlation.conditions <- list()
step04.correlation.conditions$input <- "analysis/04_Correlation_Across_Conditions.Rmd"
step04.correlation.conditions$output.dir <- "../results/04_Correlation_Across_Conditions"

# Create results directory
if (!dir.exists(step04.correlation.conditions$output.dir)) {
  dir.create(step04.correlation.conditions$output.dir, recursive=TRUE, showWarnings=FALSE)
}

# List parameters
step04.correlation.conditions$params <- list(
  list(species="LOC", species.name=species.name[["LOC"]], qc1.only=TRUE, all.nonzero.matrix=FALSE, combat=FALSE),
  list(species="ELU", species.name=species.name[["ELU"]], qc1.only=TRUE, all.nonzero.matrix=FALSE, combat=FALSE),
  list(species="DRE", species.name=species.name[["DRE"]], qc1.only=TRUE, all.nonzero.matrix=FALSE, combat=FALSE),
  list(species="LOC", species.name=species.name[["LOC"]], qc1.only=TRUE, all.nonzero.matrix=FALSE, combat=TRUE),
  list(species="ELU", species.name=species.name[["ELU"]], qc1.only=TRUE, all.nonzero.matrix=FALSE, combat=TRUE),
  list(species="DRE", species.name=species.name[["DRE"]], qc1.only=TRUE, all.nonzero.matrix=FALSE, combat=TRUE)
)

##### Run analysis script #####
for (p in (1:length(step04.correlation.conditions$params))) {
  print(paste("Species:", step04.correlation.conditions$params[[p]]$species))
  
  if (step04.correlation.conditions$params[[p]]$combat==FALSE) {
    output.dir <- file.path(step04.correlation.conditions$output.dir, "no_combat", "qc1_only")
  } else {
    output.dir <- file.path(step04.correlation.conditions$output.dir, "combat", "qc1_only")
  }
  
  if (!dir.exists(output.dir)) { dir.create(output.dir, recursive=TRUE, showWarnings=FALSE) }
  
  if (step04.correlation.conditions$params[[p]]$all.nonzero.matrix==FALSE) {
    rmarkdown::render(
      input=step04.correlation.conditions$input,
      params=step04.correlation.conditions$params[[p]],
      output_file=file.path(
        "..", output.dir,
        paste(step04.correlation.conditions$params[[p]]$species, ".html", sep="")))      
    
  } else {
    rmarkdown::render(
      input=step04.correlation.conditions$input,
      params=step04.correlation.conditions$params[[p]],
      output_file=file.path(
        "..", output.dir,
        paste(step04.correlation.conditions$params[[p]]$species, "nonzero.html", sep="_")))      
  }
}

# Clear workspace
keep <- ls()[grepl("step[0-9]*|species.name",ls())]
rm(list=setdiff(ls(), keep))
gc()


#### Step 5: Check for bimodality ####

##### Set input and parameters #####
step05.bimodality.test <- list()
step05.bimodality.test$input <- "analysis/05_Bimodality_Test.Rmd"
step05.bimodality.test$output.dir <- "../results/05_Bimodality_Test"

# Create results directory
if (!dir.exists(step05.bimodality.test$output.dir)) {
  dir.create(step05.bimodality.test$output.dir, recursive=TRUE, showWarnings=FALSE)
}

# List parameters
step05.bimodality.test$params <- list(
  list(species="LOC", species.name=species.name[["LOC"]], min.replicates=10, qc1.only=TRUE, combat=FALSE),
  list(species="ELU", species.name=species.name[["ELU"]], min.replicates=10, qc1.only=TRUE, combat=FALSE),
  list(species="DRE", species.name=species.name[["DRE"]], min.replicates=10, qc1.only=TRUE, combat=FALSE),
  list(species="LOC", species.name=species.name[["LOC"]], min.replicates=10, qc1.only=TRUE, combat=TRUE),
  list(species="ELU", species.name=species.name[["ELU"]], min.replicates=10, qc1.only=TRUE, combat=TRUE),
  list(species="DRE", species.name=species.name[["DRE"]], min.replicates=10, qc1.only=TRUE, combat=TRUE)
)

##### Run analysis script #####
for (p in (1:length(step05.bimodality.test$params))) {
  print(paste("Species:", step05.bimodality.test$params[[p]]$species))
  
  if (step05.bimodality.test$params[[p]]$combat==FALSE) {
    output.dir <- file.path(step05.bimodality.test$output.dir, "no_combat", "qc1_only")
  } else {
    output.dir <- file.path(step05.bimodality.test$output.dir, "combat", "qc1_only")
  }
  
  if (!dir.exists(output.dir)) { dir.create(output.dir, recursive=TRUE, showWarnings=FALSE) }
  
  rmarkdown::render(
    input=step05.bimodality.test$input,
    params=step05.bimodality.test$params[[p]],
    output_file=file.path(
      "..", output.dir,
      paste(step05.bimodality.test$params[[p]]$species, ".html", sep="")))
  
}

# Clear workspace
rm(list=setdiff(ls(), ls()[grepl("step[0-9]*|species.name",ls())]))
gc()


#### Step 6: (Skipped) Simulate strong bimodality ####


#### Step 7: Run GO enrichment ####

##### Set input and parameters #####
step07.go.enrichment <- list()
step07.go.enrichment$input <- "analysis/07_GO_Enrichment.Rmd"
step07.go.enrichment$output.dir <- "../results/07_GO_Enrichment"

# Create results directory
if (!dir.exists(step07.go.enrichment$output.dir)) {
  dir.create(step07.go.enrichment$output.dir, recursive=TRUE, showWarnings=FALSE)
}

# List parameters
step07.go.enrichment$params <- list(
  list(species="DRE", species.name=species.name[["DRE"]], ev.percentile=0.20, pvalue.cutoff=0.01, condition.cutoff=3, qc1.only=TRUE, combat=FALSE),
  list(species="DRE", species.name=species.name[["DRE"]], ev.percentile=0.20, pvalue.cutoff=0.01, condition.cutoff=3, qc1.only=TRUE, combat=TRUE)
)

##### Run analysis script #####
for (p in (1:length(step07.go.enrichment$params))) {
  print(paste("Species:", step07.go.enrichment$params[[p]]$species))
  
  if (step07.go.enrichment$params[[p]]$combat==FALSE) {
    output.dir <- file.path(step07.go.enrichment$output.dir, "no_combat", "qc1_only")
  } else {
    output.dir <- file.path(step07.go.enrichment$output.dir, "combat", "qc1_only")
  }
  
  if (!dir.exists(output.dir)) { dir.create(output.dir, recursive=TRUE, showWarnings=FALSE) }
  
  rmarkdown::render(
    input=step07.go.enrichment$input,
    params=step07.go.enrichment$params[[p]],
    output_file=file.path(
      "..", output.dir,
      paste(step07.go.enrichment$params[[p]]$species, ".html", sep="")))
  
}
# The output for this step is used to run GO-Figure! outside of R
# see run_go_figure.sh

# Clear workspace
rm(list=setdiff(ls(), ls()[grepl("step[0-9]*|species.name",ls())]))
gc()

#### Step 7b: Run cameraPR to test for highly/lowly variable GO terms ####

##### Set input and parameters #####
step07b.go.camera <- list()
step07b.go.camera$input <- "analysis/07b_GO_cameraPR.Rmd"
step07b.go.camera$output.dir <- "../results/07b_GO_cameraPR"

# Create results directory
if (!dir.exists(step07b.go.camera$output.dir)) {
  dir.create(step07b.go.camera$output.dir, recursive=TRUE, showWarnings=FALSE)
}

# List parameters
step07b.go.camera$params <- list(
  list(species="DRE", species.name=species.name[["DRE"]], combat=FALSE, qc1.only=FALSE, min.genes=20, max.genes=1000, condition.cutoff=3),
  list(species="DRE", species.name=species.name[["DRE"]], combat=TRUE, qc1.only=FALSE, min.genes=20, max.genes=1000, condition.cutoff=3),
  list(species="DRE", species.name=species.name[["DRE"]], combat=FALSE, qc1.only=TRUE, min.genes=20, max.genes=1000, condition.cutoff=3),   
  list(species="DRE", species.name=species.name[["DRE"]], combat=TRUE, qc1.only=TRUE, min.genes=20, max.genes=1000, condition.cutoff=3)  
)

##### Run analysis script #####
for (p in (1:length(step07b.go.camera$params))) {
  print(step07b.go.camera$params[[p]])
  
  if (step07b.go.camera$params[[p]]$combat==FALSE) {
    if (step07b.go.camera$params[[p]]$qc1.only==FALSE) {
      output.dir <- file.path(step07b.go.camera$output.dir, "no_combat", "full_qc")
    } else {
      output.dir <- file.path(step07b.go.camera$output.dir, "no_combat", "qc1_only")      
    }
  } else {
    if (step07b.go.camera$params[[p]]$qc1.only==FALSE) {
      output.dir <- file.path(step07b.go.camera$output.dir, "combat", "full_qc")      
    } else {
      output.dir <- file.path(step07b.go.camera$output.dir, "combat", "qc1_only")      
    }
  }
  
  if (!dir.exists(output.dir)) { dir.create(output.dir, recursive=TRUE, showWarnings=FALSE) }
  
  rmarkdown::render(
    input=step07b.go.camera$input,
    params=step07b.go.camera$params[[p]],
    output_file=file.path(
      "..", output.dir,
      paste(step07b.go.camera$params[[p]]$species, ".html", sep="")))
  
}
# The output for this step is used to run GO-Figure! outside of R
# see run_go_figure.sh

# Clear workspace
rm(list=setdiff(ls(), ls()[grepl("step[0-9]*|species.name",ls())]))
gc()


#### Step 8: Check for signals of selection ####

##### Set input and parameters #####
step08.selection <- list()
step08.selection$input <- "analysis/08_Selection.Rmd"
step08.selection$output.dir <- "../results/08_Selection"

# Create results directory
if (!dir.exists(step08.selection$output.dir)) {
  dir.create(step08.selection$output.dir, recursive=TRUE, showWarnings=FALSE)
}

# List parameters
step08.selection$params <- list(
  list(species="LOC", species.name=species.name[["LOC"]], ev.percentile=0.20, set.seed=12345, n.permutations=2000, qc1.only=TRUE, combat=FALSE),
  list(species="ELU", species.name=species.name[["ELU"]], ev.percentile=0.20, set.seed=12345, n.permutations=2000, qc1.only=TRUE, combat=FALSE),
  list(species="DRE", species.name=species.name[["DRE"]], ev.percentile=0.20, set.seed=12345, n.permutations=2000, qc1.only=TRUE, combat=FALSE),
  list(species="LOC", species.name=species.name[["LOC"]], ev.percentile=0.20, set.seed=12345, n.permutations=2000, qc1.only=TRUE, combat=TRUE),
  list(species="ELU", species.name=species.name[["ELU"]], ev.percentile=0.20, set.seed=12345, n.permutations=2000, qc1.only=TRUE, combat=TRUE),
  list(species="DRE", species.name=species.name[["DRE"]], ev.percentile=0.20, set.seed=12345, n.permutations=2000, qc1.only=TRUE, combat=TRUE)
)

##### Run analysis script #####
for (p in (1:length(step08.selection$params))) {
  print(paste("Species:", step08.selection$params[[p]]$species))
  
  if (step08.selection$params[[p]]$combat==FALSE) {
    output.dir <- file.path(step08.selection$output.dir, "no_combat", "qc1_only")
  } else {
    output.dir <- file.path(step08.selection$output.dir, "combat", "qc1_only")
  }
  
  if (!dir.exists(output.dir)) { dir.create(output.dir, recursive=TRUE, showWarnings=FALSE) }
  
  rmarkdown::render(
    input=step08.selection$input,
    params=step08.selection$params[[p]],
    output_file=file.path(
      "..", output.dir,
      paste(step08.selection$params[[p]]$species, ".html", sep="")))    
  
}

# Clear workspace
rm(list=setdiff(ls(), ls()[grepl("step[0-9]*|species.name",ls())]))
gc()


#### Step 9: (Skipped) Plot Selectome species tree ####


#### Step 10: Compute organ expression specificity ####

##### Set input and parameters #####
step10.expression.specificity <- list()
step10.expression.specificity$input <- "analysis/10_Expression_Specificity.Rmd"
step10.expression.specificity$output.dir <- "../results/10_Expression_Specificity"

# Create results directory
if (!dir.exists(step10.expression.specificity$output.dir)) {
  dir.create(step10.expression.specificity$output.dir, recursive=TRUE, showWarnings=FALSE)
}

# List parameters
step10.expression.specificity$params <- list(
  list(min.cpm=1.0, tau.cutoff=0.30, qc1.only=TRUE, all.nonzero.matrix=FALSE, randomize.mean.expr=FALSE, combat=FALSE),
  list(min.cpm=1.0, tau.cutoff=0.50, qc1.only=TRUE, all.nonzero.matrix=FALSE, randomize.mean.expr=FALSE, combat=FALSE),
  list(min.cpm=1.0, tau.cutoff=0.30, qc1.only=TRUE, all.nonzero.matrix=FALSE, randomize.mean.expr=FALSE, combat=TRUE),
  list(min.cpm=1.0, tau.cutoff=0.50, qc1.only=TRUE, all.nonzero.matrix=FALSE, randomize.mean.expr=FALSE, combat=TRUE)
)

##### Run analysis script #####
for (p in (1:length(step10.expression.specificity$params))) {
  
  if (step10.expression.specificity$params[[p]]$combat==FALSE) {
    output.dir <- file.path(step10.expression.specificity$output.dir, "no_combat", "qc1_only")
  } else {
    output.dir <- file.path(step10.expression.specificity$output.dir, "combat", "qc1_only")
  }
  
  if (!dir.exists(output.dir)) { dir.create(output.dir, recursive=TRUE, showWarnings=FALSE) }
  
  if (step10.expression.specificity$params[[p]]$all.nonzero.matrix==FALSE) {
    if (step10.expression.specificity$params[[p]]$randomize.mean.expr==FALSE) {
      rmarkdown::render(
        input=step10.expression.specificity$input,
        params=step10.expression.specificity$params[[p]],
        output_file=file.path(
          "..", output.dir,
          paste("LOC_ELU_DRE",
                paste0("tau", step10.expression.specificity$params[[p]]$tau.cutoff, ".html"), sep="_")))
      
    } else {
      rmarkdown::render(
        input=step10.expression.specificity$input,
        params=step10.expression.specificity$params[[p]],
        output_file=file.path(
          "..", output.dir,
          paste("LOC_ELU_DRE",
                paste0("tau", step10.expression.specificity$params[[p]]$tau.cutoff),
                paste0("rnd", step10.expression.specificity$params[[p]]$set.seed, ".html"), sep="_")))
    }
    
  } else {
    if (step10.expression.specificity$params[[p]]$randomize.mean.expr==FALSE) {
      rmarkdown::render(
        input=step10.expression.specificity$input,
        params=step10.expression.specificity$params[[p]],
        output_file=file.path(
          "..", output.dir,
          paste("LOC_ELU_DRE", "nonzero",
                paste0("tau", step10.expression.specificity$params[[p]]$tau.cutoff, ".html"), sep="_")))
      
    } else {
      rmarkdown::render(
        input=step10.expression.specificity$input,
        params=step10.expression.specificity$params[[p]],
        output_file=file.path(
          "..", output.dir,
          paste("LOC_ELU_DRE", "nonzero",
                paste0("tau", step10.expression.specificity$params[[p]]$tau.cutoff),
                paste0("rnd", step10.expression.specificity$params[[p]]$set.seed, ".html"), sep="_")))
      
    }
  }
}

# Clear workspace
keep <- ls()[grepl("step[0-9]*|species.name",ls())]
rm(list=setdiff(ls(), keep))
gc()


#### Step 11: Check selection on organ-biased genes ####

##### Set input and parameters #####
step11.selection.organ.bias <- list()
step11.selection.organ.bias$input <- "analysis/11_Selection_Organ_Bias.Rmd"
step11.selection.organ.bias$output.dir <- "../results/11_Selection_Organ_Bias"

# Create results directory
if (!dir.exists(step11.selection.organ.bias$output.dir)) {
  dir.create(step11.selection.organ.bias$output.dir, recursive=TRUE, showWarnings=FALSE)
}

# List parameters
step11.selection.organ.bias$params <- list(
  list(species="LOC", species.name=species.name[["LOC"]], tau.cutoff=0.30, qc1.only=TRUE, combat=FALSE),
  list(species="ELU", species.name=species.name[["ELU"]], tau.cutoff=0.30, qc1.only=TRUE, combat=FALSE),
  list(species="DRE", species.name=species.name[["DRE"]], tau.cutoff=0.30, qc1.only=TRUE, combat=FALSE),
  list(species="LOC", species.name=species.name[["LOC"]], tau.cutoff=0.50, qc1.only=TRUE, combat=FALSE),
  list(species="ELU", species.name=species.name[["ELU"]], tau.cutoff=0.50, qc1.only=TRUE, combat=FALSE),
  list(species="DRE", species.name=species.name[["DRE"]], tau.cutoff=0.50, qc1.only=TRUE, combat=FALSE),
  list(species="LOC", species.name=species.name[["LOC"]], tau.cutoff=0.30, qc1.only=TRUE, combat=TRUE),
  list(species="ELU", species.name=species.name[["ELU"]], tau.cutoff=0.30, qc1.only=TRUE, combat=TRUE),
  list(species="DRE", species.name=species.name[["DRE"]], tau.cutoff=0.30, qc1.only=TRUE, combat=TRUE),
  list(species="LOC", species.name=species.name[["LOC"]], tau.cutoff=0.50, qc1.only=TRUE, combat=TRUE),
  list(species="ELU", species.name=species.name[["ELU"]], tau.cutoff=0.50, qc1.only=TRUE, combat=TRUE),
  list(species="DRE", species.name=species.name[["DRE"]], tau.cutoff=0.50, qc1.only=TRUE, combat=TRUE)
)

##### Run analysis script #####
for (p in (1:length(step11.selection.organ.bias$params))) {
  print(paste("Species:", step11.selection.organ.bias$params[[p]]$species))
  
  if (step11.selection.organ.bias$params[[p]]$combat==FALSE) {
    output.dir <- file.path(step11.selection.organ.bias$output.dir, "no_combat", "qc1_only")
  } else {
    output.dir <- file.path(step11.selection.organ.bias$output.dir, "combat", "qc1_only")
  }
  
  if (!dir.exists(output.dir)) { dir.create(output.dir, recursive=TRUE, showWarnings=FALSE) }
  
  rmarkdown::render(
    input=step11.selection.organ.bias$input,
    params=step11.selection.organ.bias$params[[p]],
    output_file=file.path(
      "..", output.dir,
      paste(step11.selection.organ.bias$params[[p]]$species,
            paste0("tau", step11.selection.organ.bias$params[[p]]$tau.cutoff, ".html"), sep="_")))    
  
}

# Clear workspace
keep <- ls()[grepl("step[0-9]*|species.name",ls())]
rm(list=setdiff(ls(), keep))
gc()


#### Step 12: Check variability of organ-biased genes ####

##### Set input and parameters #####
step12.organ.bias <- list()
step12.organ.bias$input <- "analysis/12_Organ_Bias.Rmd"
step12.organ.bias$output.dir <- "../results/12_Organ_Bias"

# Create results directory
if (!dir.exists(step12.organ.bias$output.dir)) {
  dir.create(step12.organ.bias$output.dir, recursive=TRUE, showWarnings=FALSE)
}

# List parameters
step12.organ.bias$params <- list(
  list(species="LOC", species.name=species.name[["LOC"]], tau.cutoff=0.30, qc1.only=TRUE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=FALSE),
  list(species="ELU", species.name=species.name[["ELU"]], tau.cutoff=0.30, qc1.only=TRUE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=FALSE),
  list(species="DRE", species.name=species.name[["DRE"]], tau.cutoff=0.30, qc1.only=TRUE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=FALSE),
  list(species="LOC", species.name=species.name[["LOC"]], tau.cutoff=0.50, qc1.only=TRUE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=FALSE),
  list(species="ELU", species.name=species.name[["ELU"]], tau.cutoff=0.50, qc1.only=TRUE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=FALSE),
  list(species="DRE", species.name=species.name[["DRE"]], tau.cutoff=0.50, qc1.only=TRUE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=FALSE),
  list(species="LOC", species.name=species.name[["LOC"]], tau.cutoff=0.30, qc1.only=TRUE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=TRUE),
  list(species="ELU", species.name=species.name[["ELU"]], tau.cutoff=0.30, qc1.only=TRUE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=TRUE),
  list(species="DRE", species.name=species.name[["DRE"]], tau.cutoff=0.30, qc1.only=TRUE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=TRUE),
  list(species="LOC", species.name=species.name[["LOC"]], tau.cutoff=0.50, qc1.only=TRUE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=TRUE),
  list(species="ELU", species.name=species.name[["ELU"]], tau.cutoff=0.50, qc1.only=TRUE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=TRUE),
  list(species="DRE", species.name=species.name[["DRE"]], tau.cutoff=0.50, qc1.only=TRUE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=TRUE)
)

##### Run analysis script #####
for (p in (1:length(step12.organ.bias$params))) {
  print(paste("Species:", step12.organ.bias$params[[p]]$species))
  
  if (step12.organ.bias$params[[p]]$combat==FALSE) {
    output.dir <- file.path(step12.organ.bias$output.dir, "no_combat", "qc1_only")
  } else {
    output.dir <- file.path(step12.organ.bias$output.dir, "combat", "qc1_only")
  }
  
  if (!dir.exists(output.dir)) { dir.create(output.dir, recursive=TRUE, showWarnings=FALSE) }
  
  if (step12.organ.bias$params[[p]]$all.nonzero.matrix==FALSE) {
    if (step12.organ.bias$params[[p]]$randomized.mean.expr==FALSE) {
      rmarkdown::render(
        input=step12.organ.bias$input,
        params=step12.organ.bias$params[[p]],
        output_file=file.path(
          "..", output.dir,
          paste(step12.organ.bias$params[[p]]$species,
                paste0("tau", step12.organ.bias$params[[p]]$tau.cutoff, ".html"), sep="_")))
      
    } else {
      rmarkdown::render(
        input=step12.organ.bias$input,
        params=step12.organ.bias$params[[p]],
        output_file=file.path(
          "..", output.dir,
          paste(step12.organ.bias$params[[p]]$species,
                paste0("tau", step12.organ.bias$params[[p]]$tau.cutoff),
                paste0("rnd", step12.organ.bias$params[[p]]$seed.run, ".html"), sep="_")))
      
    }
  } else {
    if (step12.organ.bias$params[[p]]$randomized.mean.expr==FALSE) {
      rmarkdown::render(
        input=step12.organ.bias$input,
        params=step12.organ.bias$params[[p]],
        output_file=file.path(
          "..", output.dir,
          paste(step12.organ.bias$params[[p]]$species, "nonzero",
                paste0("tau", step12.organ.bias$params[[p]]$tau.cutoff, ".html"), sep="_")))
      
    } else {
      rmarkdown::render(
        input=step12.organ.bias$input,
        params=step12.organ.bias$params[[p]],
        output_file=file.path(
          "..", output.dir,
          paste(step12.organ.bias$params[[p]]$species, "nonzero",
                paste0("tau", step12.organ.bias$params[[p]]$tau.cutoff),
                paste0("rnd", step12.organ.bias$params[[p]]$seed.run, ".html"), sep="_")))
      
    }
  }
}


# Clear workspace
keep <- ls()[grepl("step[0-9]*|species.name",ls())]
rm(list=setdiff(ls(), keep))
gc()


#### Step 12b: Check variability (residual MAD) of organ-biased genes ####

##### Set input and parameters #####
step12b.organ.bias.nojack <- list()
step12b.organ.bias.nojack$input <- "analysis/12b_Organ_Bias_NoJack.Rmd"
step12b.organ.bias.nojack$output.dir <- "../results/12b_Organ_Bias_NoJack"

# Create results directory
if (!dir.exists(step12b.organ.bias.nojack$output.dir)) {
  dir.create(step12b.organ.bias.nojack$output.dir, recursive=TRUE, showWarnings=FALSE)
}

# List parameters
step12b.organ.bias.nojack$params <- list(
  list(species="LOC", species.name=species.name[["LOC"]], tau.cutoff=0.30, qc1.only=FALSE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=FALSE),
  list(species="ELU", species.name=species.name[["ELU"]], tau.cutoff=0.30, qc1.only=FALSE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=FALSE),
  list(species="DRE", species.name=species.name[["DRE"]], tau.cutoff=0.30, qc1.only=FALSE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=FALSE),
  list(species="LOC", species.name=species.name[["LOC"]], tau.cutoff=0.50, qc1.only=FALSE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=FALSE),
  list(species="ELU", species.name=species.name[["ELU"]], tau.cutoff=0.50, qc1.only=FALSE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=FALSE),
  list(species="DRE", species.name=species.name[["DRE"]], tau.cutoff=0.50, qc1.only=FALSE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=FALSE),
  list(species="LOC", species.name=species.name[["LOC"]], tau.cutoff=0.30, qc1.only=FALSE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=TRUE),
  list(species="ELU", species.name=species.name[["ELU"]], tau.cutoff=0.30, qc1.only=FALSE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=TRUE),
  list(species="DRE", species.name=species.name[["DRE"]], tau.cutoff=0.30, qc1.only=FALSE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=TRUE),
  list(species="LOC", species.name=species.name[["LOC"]], tau.cutoff=0.50, qc1.only=FALSE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=TRUE),
  list(species="ELU", species.name=species.name[["ELU"]], tau.cutoff=0.50, qc1.only=FALSE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=TRUE),
  list(species="DRE", species.name=species.name[["DRE"]], tau.cutoff=0.50, qc1.only=FALSE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=TRUE),
  list(species="LOC", species.name=species.name[["LOC"]], tau.cutoff=0.30, qc1.only=TRUE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=FALSE),
  list(species="ELU", species.name=species.name[["ELU"]], tau.cutoff=0.30, qc1.only=TRUE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=FALSE),
  list(species="DRE", species.name=species.name[["DRE"]], tau.cutoff=0.30, qc1.only=TRUE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=FALSE),
  list(species="LOC", species.name=species.name[["LOC"]], tau.cutoff=0.50, qc1.only=TRUE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=FALSE),
  list(species="ELU", species.name=species.name[["ELU"]], tau.cutoff=0.50, qc1.only=TRUE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=FALSE),
  list(species="DRE", species.name=species.name[["DRE"]], tau.cutoff=0.50, qc1.only=TRUE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=FALSE),
  list(species="LOC", species.name=species.name[["LOC"]], tau.cutoff=0.30, qc1.only=TRUE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=TRUE),
  list(species="ELU", species.name=species.name[["ELU"]], tau.cutoff=0.30, qc1.only=TRUE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=TRUE),
  list(species="DRE", species.name=species.name[["DRE"]], tau.cutoff=0.30, qc1.only=TRUE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=TRUE),
  list(species="LOC", species.name=species.name[["LOC"]], tau.cutoff=0.50, qc1.only=TRUE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=TRUE),
  list(species="ELU", species.name=species.name[["ELU"]], tau.cutoff=0.50, qc1.only=TRUE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=TRUE),
  list(species="DRE", species.name=species.name[["DRE"]], tau.cutoff=0.50, qc1.only=TRUE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=TRUE)
)

##### Run analysis script #####
for (p in (1:length(step12b.organ.bias.nojack$params))) {
  print(unlist(step12b.organ.bias.nojack$params[[p]]))
  
  if (step12b.organ.bias.nojack$params[[p]]$combat==FALSE) {
    if (step12b.organ.bias.nojack$params[[p]]$qc1.only==TRUE) {
      output.dir <- file.path(step12b.organ.bias.nojack$output.dir, "no_combat", "qc1_only")      
    } else {
      output.dir <- file.path(step12b.organ.bias.nojack$output.dir, "no_combat", "full_qc")      
    }
  } else {
    if (step12b.organ.bias.nojack$params[[p]]$qc1.only==TRUE) {
      output.dir <- file.path(step12b.organ.bias.nojack$output.dir, "combat", "qc1_only")      
    } else {
      output.dir <- file.path(step12b.organ.bias.nojack$output.dir, "combat", "full_qc")
    }
  }
  
  if (!dir.exists(output.dir)) { dir.create(output.dir, recursive=TRUE, showWarnings=FALSE) }
  
  if (step12b.organ.bias.nojack$params[[p]]$all.nonzero.matrix==FALSE) {
    rmarkdown::render(
      input=step12b.organ.bias.nojack$input,
      params=step12b.organ.bias.nojack$params[[p]],
      output_file=file.path(
        "..", output.dir,
        paste(step12b.organ.bias.nojack$params[[p]]$species,
              paste0("tau", step12b.organ.bias.nojack$params[[p]]$tau.cutoff,".html"), sep="_")))
    
    } else {
      rmarkdown::render(
        input=step12b.organ.bias.nojack$input,
        params=step12b.organ.bias.nojack$params[[p]],
        output_file=file.path(
          "..", output.dir,
          paste(step12b.organ.bias.nojack$params[[p]]$species, "nonzero",
                paste0("tau", step12b.organ.bias.nojack$params[[p]]$tau.cutoff,".html"), sep="_")))
    }
}

# Clear workspace
keep <- ls()[grepl("step[0-9]*|species.name",ls())]
rm(list=setdiff(ls(), keep))
gc()


#### Step 13: Combine results from step 12 across species ####

##### Set input and parameters #####
step13.organ.bias.combined <- list()
step13.organ.bias.combined$input <- "analysis/13_Organ_Bias_Combined.Rmd"
step13.organ.bias.combined$output.dir <- "../results/13_Organ_Bias_Combined"

# Create results directory
if (!dir.exists(step13.organ.bias.combined$output.dir)) {
  dir.create(step13.organ.bias.combined$output.dir, recursive=TRUE, showWarnings=FALSE)
}

# List parameters
step13.organ.bias.combined$params <- list(
  list(tau.cutoff=0.30, qc1.only=TRUE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=FALSE),
  list(tau.cutoff=0.50, qc1.only=TRUE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=FALSE),
  list(tau.cutoff=0.30, qc1.only=TRUE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=TRUE),
  list(tau.cutoff=0.50, qc1.only=TRUE, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE, combat=TRUE)
)

##### Run analysis script #####
for (p in (1:length(step13.organ.bias.combined$params))) {
  
  if (step13.organ.bias.combined$params[[p]]$combat==FALSE) {
    output.dir <- file.path(step13.organ.bias.combined$output.dir, "no_combat", "qc1_only")
  } else {
    output.dir <- file.path(step13.organ.bias.combined$output.dir, "combat", "qc1_only")
  }
  
  if (!dir.exists(output.dir)) { dir.create(output.dir, recursive=TRUE, showWarnings=FALSE) }
  
  if (step13.organ.bias.combined$params[[p]]$all.nonzero.matrix==FALSE) {
    if (step13.organ.bias.combined$params[[p]]$randomized.mean.expr==FALSE) {
      rmarkdown::render(
        input=step13.organ.bias.combined$input,
        params=step13.organ.bias.combined$params[[p]],
        output_file=file.path(
          "..", output.dir,
          paste("LOC_ELU_DRE",
                paste0("tau", step13.organ.bias.combined$params[[p]]$tau.cutoff, ".html"), sep="_")))
      
    } else {
      rmarkdown::render(
        input=step13.organ.bias.combined$input,
        params=step13.organ.bias.combined$params[[p]],
        output_file=file.path(
          "..", output.dir,
          paste("LOC_ELU_DRE",
                paste0("tau", step13.organ.bias.combined$params[[p]]$tau.cutoff),
                paste0("rnd", step13.organ.bias.combined$params[[p]]$seed.run, ".html"), sep="_")))
      
    }
  } else {
    if (step13.organ.bias.combined$params[[p]]$randomized.mean.expr==FALSE) {
      rmarkdown::render(
        input=step13.organ.bias.combined$input,
        params=step13.organ.bias.combined$params[[p]],
        output_file=file.path(
          "..", output.dir,
          paste("LOC_ELU_DRE", "nonzero",
                paste0("tau", step13.organ.bias.combined$params[[p]]$tau.cutoff, ".html"), sep="_")))
      
    } else {
      rmarkdown::render(
        input=step13.organ.bias.combined$input,
        params=step13.organ.bias.combined$params[[p]],
        output_file=file.path(
          "..", output.dir,
          paste("LOC_ELU_DRE", "nonzero", 
                paste0("tau", step13.organ.bias.combined$params[[p]]$tau.cutoff),
                paste0("rnd", step13.organ.bias.combined$params[[p]]$seed.run, ".html"), sep="_")))
      
    }
  }
}


# Clear workspace
keep <- ls()[grepl("step[0-9]*|species.name",ls())]
rm(list=setdiff(ls(), keep))
gc()


#### Step 14: (Skipped) Combine results from step 12 for observed and randomized data ####


#### Step 15: Run WCGNA to check for expression covariation ####

##### Set input and parameters #####
step15.wgcna <- list()
step15.wgcna$input <- "analysis/15_WGCNA.Rmd"
step15.wgcna$output.dir <- "../results/15_WGCNA"

# Create results directory
if (!dir.exists(step15.wgcna$output.dir)) {
  dir.create(step15.wgcna$output.dir, recursive=TRUE, showWarnings=FALSE)
}

# List parameters
step15.wgcna$params <- list(
  list(species="LOC", species.name=species.name[["LOC"]], tau.cutoff=0.30, min.replicates=10, qc1.only=FALSE, combat=FALSE, set.seed=12345, power=10),
  list(species="ELU", species.name=species.name[["ELU"]], tau.cutoff=0.30, min.replicates=10, qc1.only=FALSE, combat=FALSE, set.seed=12345, power=10),
  list(species="DRE", species.name=species.name[["DRE"]], tau.cutoff=0.30, min.replicates=10, qc1.only=FALSE, combat=FALSE, set.seed=12345, power=10),
  #list(species="LOC", species.name=species.name[["LOC"]], tau.cutoff=0.30, min.replicates=10, qc1.only=FALSE, combat=FALSE, set.seed=67890, power=10),
  #list(species="ELU", species.name=species.name[["ELU"]], tau.cutoff=0.30, min.replicates=10, qc1.only=FALSE, combat=FALSE, set.seed=67890, power=10),
  #list(species="DRE", species.name=species.name[["DRE"]], tau.cutoff=0.30, min.replicates=10, qc1.only=FALSE, combat=FALSE, set.seed=67890, power=10),
  list(species="LOC", species.name=species.name[["LOC"]], tau.cutoff=0.30, min.replicates=10, qc1.only=FALSE, combat=TRUE, set.seed=12345, power=10),
  list(species="ELU", species.name=species.name[["ELU"]], tau.cutoff=0.30, min.replicates=10, qc1.only=FALSE, combat=TRUE, set.seed=12345, power=10),
  list(species="DRE", species.name=species.name[["DRE"]], tau.cutoff=0.30, min.replicates=10, qc1.only=FALSE, combat=TRUE, set.seed=12345, power=10),
  #list(species="LOC", species.name=species.name[["LOC"]], tau.cutoff=0.30, min.replicates=10, qc1.only=FALSE, combat=TRUE, set.seed=67890, power=10),
  #list(species="ELU", species.name=species.name[["ELU"]], tau.cutoff=0.30, min.replicates=10, qc1.only=FALSE, combat=TRUE, set.seed=67890, power=10),
  #list(species="DRE", species.name=species.name[["DRE"]], tau.cutoff=0.30, min.replicates=10, qc1.only=FALSE, combat=TRUE, set.seed=67890, power=10),
  list(species="LOC", species.name=species.name[["LOC"]], tau.cutoff=0.50, min.replicates=10, qc1.only=FALSE, combat=FALSE, set.seed=12345, power=10),
  list(species="ELU", species.name=species.name[["ELU"]], tau.cutoff=0.50, min.replicates=10, qc1.only=FALSE, combat=FALSE, set.seed=12345, power=10),
  list(species="DRE", species.name=species.name[["DRE"]], tau.cutoff=0.50, min.replicates=10, qc1.only=FALSE, combat=FALSE, set.seed=12345, power=10),
  #list(species="LOC", species.name=species.name[["LOC"]], tau.cutoff=0.50, min.replicates=10, qc1.only=FALSE, combat=FALSE, set.seed=67890, power=10),
  #list(species="ELU", species.name=species.name[["ELU"]], tau.cutoff=0.50, min.replicates=10, qc1.only=FALSE, combat=FALSE, set.seed=67890, power=10),
  #list(species="DRE", species.name=species.name[["DRE"]], tau.cutoff=0.50, min.replicates=10, qc1.only=FALSE, combat=FALSE, set.seed=67890, power=10),
  list(species="LOC", species.name=species.name[["LOC"]], tau.cutoff=0.50, min.replicates=10, qc1.only=FALSE, combat=TRUE, set.seed=12345, power=10),
  list(species="ELU", species.name=species.name[["ELU"]], tau.cutoff=0.50, min.replicates=10, qc1.only=FALSE, combat=TRUE, set.seed=12345, power=10),
  list(species="DRE", species.name=species.name[["DRE"]], tau.cutoff=0.50, min.replicates=10, qc1.only=FALSE, combat=TRUE, set.seed=12345, power=10)#,
  #list(species="LOC", species.name=species.name[["LOC"]], tau.cutoff=0.50, min.replicates=10, qc1.only=FALSE, combat=TRUE, set.seed=67890, power=10),
  #list(species="ELU", species.name=species.name[["ELU"]], tau.cutoff=0.50, min.replicates=10, qc1.only=FALSE, combat=TRUE, set.seed=67890, power=10),
  #list(species="DRE", species.name=species.name[["DRE"]], tau.cutoff=0.50, min.replicates=10, qc1.only=FALSE, combat=TRUE, set.seed=67890, power=10)
)


##### Run analysis script #####
for (p in (1:length(step15.wgcna$params))) {
  print(unlist(step15.wgcna$params[[p]]))
  
  if (step15.wgcna$params[[p]]$combat==FALSE) {
    if (step15.wgcna$params[[p]]$qc1.only==TRUE) {
      output.dir <- file.path(step15.wgcna$output.dir, "no_combat", "qc1_only")      
    } else {
      output.dir <- file.path(step15.wgcna$output.dir, "no_combat", "full_qc")      
    }
  } else {
    if (step15.wgcna$params[[p]]$qc1.only==TRUE) {
      output.dir <- file.path(step15.wgcna$output.dir, "combat", "qc1_only")      
    } else {
      output.dir <- file.path(step15.wgcna$output.dir, "combat", "full_qc")
    }
  }
  
  if (!dir.exists(output.dir)) { dir.create(output.dir, recursive=TRUE, showWarnings=FALSE) }
  
  rmarkdown::render(
    input=step15.wgcna$input,
    params=step15.wgcna$params[[p]],
    output_file=file.path(
      "..", output.dir,
      paste(step15.wgcna$params[[p]]$species,
            paste0("tau", step15.wgcna$params[[p]]$tau.cutoff),
            paste0("rnd", step15.wgcna$params[[p]]$set.seed,".html"), sep="_")))
  
}

# Clear workspace
keep <- ls()[grepl("step[0-9]*|species.name",ls())]
rm(list=setdiff(ls(), keep))
gc()


#### Save all parameters ####
save.image(file.path("../results", "run_pipeline_revisions.Rdata"))