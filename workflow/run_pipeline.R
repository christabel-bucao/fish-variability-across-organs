#### Control package versions ####
#renv::init()
BiocManager::install(version="3.13", force=TRUE)

renv::restore() # Restore CRAN package versions specified in conda environment
BiocManager::install(c("edgeR","limma","GOstats","org.Dr.eg.db"), force=TRUE)
BiocManager::install(
  c("KEGGREST","genefilter","Category","RBGL","Rgraphviz","BiocGenerics", 
    "zlibbioc","Biostrings","Biobase","IRanges","GenomeInfoDb","AnnotationDbi",
    "GSEABase","S4Vectors", "graph","annotate","XVector","GO.db","AnnotationForge"),
  force=TRUE) # For consistency, other packages should have the same Bioconductor version
renv::snapshot(packages=c("edgeR","limma","GOstats","org.Dr.eg.db"))

BiocManager::install(c("ggtree","tidytree","treeio"), force=TRUE) # For plotting species tree
renv::snapshot(packages=c("ggtree","tidytree","treeio"))

# For any discrepancies, please refer to session info included in all HTML output files

#### Set species ####
species.name <- c("Lepisosteus oculatus","Esox lucius","Danio rerio")
names(species.name) <- c("LOC","ELU","DRE")
  
#### Step 1: Filter samples ####

##### Set input and parameters #####
step01.filter.samples <- list()
step01.filter.samples$input <- "analysis/01_Filter_Samples.Rmd"
step01.filter.samples$output.dir <- "../results/01_Filter_Samples"

# Create results directory
if (!dir.exists(step01.filter.samples$output.dir)) {
  dir.create(step01.filter.samples$output.dir, recursive=TRUE)
}

# List parameters
## For expression variability analysis, we require a minimum of 4 replicates per condition
## For organ expression specificity analysis, we require 2 replicates for spotted gar and 1 replicate for zebrafish and northern pike
step01.filter.samples$params <- list(
  list(species="LOC", species.name=species.name[["LOC"]], min.genes=6000, min.reads=0.2e6, min.replicates=4, min.cpm=1.0, normalize.by.condition=TRUE, set.scale=FALSE),
  list(species="ELU", species.name=species.name[["ELU"]], min.genes=10000, min.reads=0.3e6, min.replicates=4, min.cpm=1.0, normalize.by.condition=TRUE, set.scale=TRUE),
  list(species="DRE", species.name=species.name[["DRE"]], min.genes=10000, min.reads=0.5e6, min.replicates=4, min.cpm=1.0, normalize.by.condition=TRUE, set.scale=TRUE),
  list(species="LOC", species.name=species.name[["LOC"]], min.genes=6000, min.reads=0.2e6, min.replicates=2, min.cpm=1.0, normalize.by.condition=FALSE, set.scale=FALSE),
  list(species="ELU", species.name=species.name[["ELU"]], min.genes=10000, min.reads=0.3e6, min.replicates=1, min.cpm=1.0, normalize.by.condition=FALSE, set.scale=TRUE),
  list(species="DRE", species.name=species.name[["DRE"]], min.genes=10000, min.reads=0.5e6, min.replicates=1, min.cpm=1.0, normalize.by.condition=FALSE, set.scale=TRUE)
)

##### Run analysis script #####
for (p in (1:length(step01.filter.samples$params))) {
  print(paste("Species:", step01.filter.samples$params[[p]]$species))
  print(paste("Minimum replicates:", step01.filter.samples$params[[p]]$min.replicates))
  
  rmarkdown::render(
    input=step01.filter.samples$input,
    params=step01.filter.samples$params[[p]],
    output_file=file.path(
      "..", step01.filter.samples$output.dir,
      paste(step01.filter.samples$params[[p]]$species, step01.filter.samples$params[[p]]$min.replicates, sep="_", "reps.html")))
}

# Clear workspace
rm(list=setdiff(ls(), ls()[grepl("step[0-9]*|species.name",ls())]))
gc()


#### Step 2: Estimate gene expression variability ####

##### Set input and parameters #####
step02.variability.jackknife <- list()
step02.variability.jackknife$input <- "analysis/02_Variability_Jackknife.Rmd"
step02.variability.jackknife$output.dir <- "../results/02_Variability_Jackknife"

# Create results directory
if (!dir.exists(step02.variability.jackknife$output.dir)) {
  dir.create(step02.variability.jackknife$output.dir, recursive=TRUE)
}

# List parameters
step02.variability.jackknife$params <- list(
  list(species="DRE", species.name=species.name[["DRE"]], min.percentile=0.00, max.percentile=0.95, win.size=100, all.nonzero.matrix=FALSE, run.example=TRUE, ex.tissue="brain", ex.sex="F"), 
  list(species="LOC", species.name=species.name[["LOC"]], min.percentile=0.00, max.percentile=0.95, win.size=100, all.nonzero.matrix=FALSE, run.example=FALSE),
  list(species="ELU", species.name=species.name[["ELU"]], min.percentile=0.00, max.percentile=0.95, win.size=100, all.nonzero.matrix=FALSE, run.example=FALSE),
  list(species="DRE", species.name=species.name[["DRE"]], min.percentile=0.00, max.percentile=0.95, win.size=100, all.nonzero.matrix=FALSE, run.example=FALSE),
  list(species="DRE", species.name=species.name[["DRE"]], min.percentile=0.00, max.percentile=0.95, win.size=100, all.nonzero.matrix=TRUE, run.example=TRUE, ex.tissue="brain", ex.sex="F"),
  list(species="LOC", species.name=species.name[["LOC"]], min.percentile=0.00, max.percentile=0.95, win.size=100, all.nonzero.matrix=TRUE, run.example=FALSE),
  list(species="ELU", species.name=species.name[["ELU"]], min.percentile=0.00, max.percentile=0.95, win.size=100, all.nonzero.matrix=TRUE, run.example=FALSE),
  list(species="DRE", species.name=species.name[["DRE"]], min.percentile=0.00, max.percentile=0.95, win.size=100, all.nonzero.matrix=TRUE, run.example=FALSE)
)
# If all.nonzero.matrix=TRUE, we consider only the set of genes with nonzero counts across all samples
# The last one is for an example computation without jackknife resampling

##### Run analysis script #####
for (p in (1:length(step02.variability.jackknife$params))) {
  print(paste("Species:", step02.variability.jackknife$params[[p]]$species))
  print(paste("Run jackknife?:", !step02.variability.jackknife$params[[p]]$run.example))
  
  if (step02.variability.jackknife$params[[p]]$all.nonzero.matrix==FALSE) {
    if (step02.variability.jackknife$params[[p]]$run.example==TRUE) {
      rmarkdown::render(
        input=step02.variability.jackknife$input,
        params=step02.variability.jackknife$params[[p]],
        output_file=file.path(
          "..", step02.variability.jackknife$output.dir,
          paste(step02.variability.jackknife$params[[p]]$species, "example.html", sep="_")))        
    } else {
      rmarkdown::render(
        input=step02.variability.jackknife$input,
        params=step02.variability.jackknife$params[[p]],
        output_file=file.path(
          "..", step02.variability.jackknife$output.dir,
          paste(step02.variability.jackknife$params[[p]]$species, ".html", sep="")))    
    }    
  } else {
    if (step02.variability.jackknife$params[[p]]$run.example==TRUE) {
      rmarkdown::render(
        input=step02.variability.jackknife$input,
        params=step02.variability.jackknife$params[[p]],
        output_file=file.path(
          "..", step02.variability.jackknife$output.dir,
          paste(step02.variability.jackknife$params[[p]]$species, "nonzero", "example.html", sep="_")))        
    } else {
      rmarkdown::render(
        input=step02.variability.jackknife$input,
        params=step02.variability.jackknife$params[[p]],
        output_file=file.path(
          "..", step02.variability.jackknife$output.dir,
          paste(step02.variability.jackknife$params[[p]]$species, "nonzero.html", sep="_")))    
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
  dir.create(step03.variability.plots$output.dir, recursive=TRUE)
}

# List parameters
step03.variability.plots$params <- list(
  list(species="LOC", species.name=species.name[["LOC"]], min.percentile=0.00, max.percentile=0.95, all.nonzero.matrix=FALSE),
  list(species="ELU", species.name=species.name[["ELU"]], min.percentile=0.00, max.percentile=0.95, all.nonzero.matrix=FALSE),
  list(species="DRE", species.name=species.name[["DRE"]], min.percentile=0.00, max.percentile=0.95, all.nonzero.matrix=FALSE),
  list(species="LOC", species.name=species.name[["LOC"]], min.percentile=0.00, max.percentile=0.95, all.nonzero.matrix=TRUE),
  list(species="ELU", species.name=species.name[["ELU"]], min.percentile=0.00, max.percentile=0.95, all.nonzero.matrix=TRUE),
  list(species="DRE", species.name=species.name[["DRE"]], min.percentile=0.00, max.percentile=0.95, all.nonzero.matrix=TRUE)
)

##### Run analysis script #####
for (p in (1:length(step03.variability.plots$params))) {
  print(paste("Species:", step03.variability.plots$params[[p]]$species))
  
  if (step03.variability.plots$params[[p]]$all.nonzero.matrix==FALSE) {
    rmarkdown::render(
      input=step03.variability.plots$input,
      params=step03.variability.plots$params[[p]],
      output_file=file.path(
        "..", step03.variability.plots$output.dir,
        paste(step03.variability.plots$params[[p]]$species, ".html", sep="")))    
  } else {
    rmarkdown::render(
      input=step03.variability.plots$input,
      params=step03.variability.plots$params[[p]],
      output_file=file.path(
        "..", step03.variability.plots$output.dir,
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
  dir.create(step04.correlation.conditions$output.dir, recursive=TRUE)
}

# List parameters
step04.correlation.conditions$params <- list(
  list(species="LOC", species.name=species.name[["LOC"]], all.nonzero.matrix=FALSE),
  list(species="ELU", species.name=species.name[["ELU"]], all.nonzero.matrix=FALSE),
  list(species="DRE", species.name=species.name[["DRE"]], all.nonzero.matrix=FALSE),
  list(species="LOC", species.name=species.name[["LOC"]], all.nonzero.matrix=TRUE),
  list(species="ELU", species.name=species.name[["ELU"]], all.nonzero.matrix=TRUE),
  list(species="DRE", species.name=species.name[["DRE"]], all.nonzero.matrix=TRUE)
)

##### Run analysis script #####
for (p in (1:length(step04.correlation.conditions$params))) {
  print(paste("Species:", step04.correlation.conditions$params[[p]]$species))
  
  if (step04.correlation.conditions$params[[p]]$all.nonzero.matrix==FALSE) {
    rmarkdown::render(
      input=step04.correlation.conditions$input,
      params=step04.correlation.conditions$params[[p]],
      output_file=file.path(
        "..", step04.correlation.conditions$output.dir,
        paste(step04.correlation.conditions$params[[p]]$species, ".html", sep="")))    
  } else {
    rmarkdown::render(
      input=step04.correlation.conditions$input,
      params=step04.correlation.conditions$params[[p]],
      output_file=file.path(
        "..", step04.correlation.conditions$output.dir,
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
  dir.create(step05.bimodality.test$output.dir, recursive=TRUE)
}

# List parameters
step05.bimodality.test$params <- list(
  list(species="LOC", species.name=species.name[["LOC"]], min.replicates=10),
  list(species="ELU", species.name=species.name[["ELU"]], min.replicates=10),
  list(species="DRE", species.name=species.name[["DRE"]], min.replicates=10)
)

##### Run analysis script #####
for (p in (1:length(step05.bimodality.test$params))) {
  print(paste("Species:", step05.bimodality.test$params[[p]]$species))
  
  rmarkdown::render(
    input=step05.bimodality.test$input,
    params=step05.bimodality.test$params[[p]],
    output_file=file.path(
      "..", step05.bimodality.test$output.dir,
      paste(step05.bimodality.test$params[[p]]$species, ".html", sep="")))   
}

# Clear workspace
rm(list=setdiff(ls(), ls()[grepl("step[0-9]*|species.name",ls())]))
gc()


#### Step 6: Simulate strong bimodality ####

##### Set input and parameters #####
step06.simulated.bimodality <- list()
step06.simulated.bimodality$input <- "analysis/06_Simulated_Bimodality.Rmd"
step06.simulated.bimodality$output.dir <- "../results/06_Simulated_Bimodality"

# Create results directory
if (!dir.exists(step06.simulated.bimodality$output.dir)) {
  dir.create(step06.simulated.bimodality$output.dir, recursive=TRUE)
}

# List parameters
step06.simulated.bimodality$params <- list(
  list(species="ELU", species.name=species.name[["ELU"]], sex="F", total.replicates=10, tissue1="brain", tissue2="gonads",
       set.seed=12345, min.cpm=1.0, min.percentile=0.00, max.percentile=0.95, win.size=100)
)

##### Run analysis script #####
for (p in (1:length(step06.simulated.bimodality$params))) {
  print(paste("Species:", step06.simulated.bimodality$params[[p]]$species))
  
  rmarkdown::render(
    input=step06.simulated.bimodality$input,
    params=step06.simulated.bimodality$params[[p]],
    output_file=file.path(
      "..", step06.simulated.bimodality$output.dir,
      paste(step06.simulated.bimodality$params[[p]]$species, ".html", sep="")))   
}

# Clear workspace
rm(list=setdiff(ls(), ls()[grepl("step[0-9]*|species.name",ls())]))
gc()


#### Step 7: Run GO enrichment ####

##### Set input and parameters #####
step07.go.enrichment <- list()
step07.go.enrichment$input <- "analysis/07_GO_Enrichment.Rmd"
step07.go.enrichment$output.dir <- "../results/07_GO_Enrichment"

# Create results directory
if (!dir.exists(step07.go.enrichment$output.dir)) {
  dir.create(step07.go.enrichment$output.dir, recursive=TRUE)
}

# List parameters
step07.go.enrichment$params <- list(
  list(species="DRE", species.name=species.name[["DRE"]], ev.percentile=0.20, pvalue.cutoff=0.01, condition.cutoff=3)
)

##### Run analysis script #####
for (p in (1:length(step07.go.enrichment$params))) {
  print(paste("Species:", step07.go.enrichment$params[[p]]$species))
  
  rmarkdown::render(
    input=step07.go.enrichment$input,
    params=step07.go.enrichment$params[[p]],
    output_file=file.path(
      "..", step07.go.enrichment$output.dir,
      paste(step07.go.enrichment$params[[p]]$species, ".html", sep="")))   
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
  dir.create(step08.selection$output.dir, recursive=TRUE)
}

# List parameters
step08.selection$params <- list(
  list(species="LOC", species.name=species.name[["LOC"]], ev.percentile=0.20, set.seed=12345, n.permutations=2000),
  list(species="ELU", species.name=species.name[["ELU"]], ev.percentile=0.20, set.seed=12345, n.permutations=2000),
  list(species="DRE", species.name=species.name[["DRE"]], ev.percentile=0.20, set.seed=12345, n.permutations=2000)
)

##### Run analysis script #####
for (p in (1:length(step08.selection$params))) {
  print(paste("Species:", step08.selection$params[[p]]$species))
  
  rmarkdown::render(
    input=step08.selection$input,
    params=step08.selection$params[[p]],
    output_file=file.path(
      "..", step08.selection$output.dir,
      paste(step08.selection$params[[p]]$species, ".html", sep="")))   
}

# Clear workspace
rm(list=setdiff(ls(), ls()[grepl("step[0-9]*|species.name",ls())]))
gc()


#### Step 9: Plot Selectome species tree ####

##### Set input and parameters #####
step09.selection.species.tree <- list()
step09.selection.species.tree$input <- "analysis/09_Selection_Species_Tree.Rmd"
step09.selection.species.tree$output.dir <- "../results/09_Selection_Species_Tree"

# Create results directory
if (!dir.exists(step09.selection.species.tree$output.dir)) {
  dir.create(step09.selection.species.tree$output.dir, recursive=TRUE)
}

##### Run analysis script #####
# No added parameters
rmarkdown::render(
  input=step09.selection.species.tree$input,
  output_file=file.path(
    "..", step09.selection.species.tree$output.dir,
    paste("Selection_Species_Tree", ".html", sep="")))

# Clear workspace
rm(list=setdiff(ls(), ls()[grepl("step[0-9]*|species.name",ls())]))
gc()


#### Step 10: Compute organ expression specificity ####

##### Set input and parameters #####
step10.expression.specificity <- list()
step10.expression.specificity$input <- "analysis/10_Expression_Specificity.Rmd"
step10.expression.specificity$output.dir <- "../results/10_Expression_Specificity"

# Create results directory
if (!dir.exists(step10.expression.specificity$output.dir)) {
  dir.create(step10.expression.specificity$output.dir, recursive=TRUE)
}

# List parameters
step10.expression.specificity$params <- list(
  list(min.cpm=1.0, tau.cutoff=0.30, all.nonzero.matrix=FALSE, randomize.mean.expr=FALSE),
  list(min.cpm=1.0, tau.cutoff=0.30, all.nonzero.matrix=TRUE, randomize.mean.expr=FALSE),
  list(min.cpm=1.0, tau.cutoff=0.30, all.nonzero.matrix=FALSE, randomize.mean.expr=TRUE, set.seed=12345),
  list(min.cpm=1.0, tau.cutoff=0.30, all.nonzero.matrix=FALSE, randomize.mean.expr=TRUE, set.seed=67890)
)

##### Run analysis script #####
for (p in (1:length(step10.expression.specificity$params))) {
  if (step10.expression.specificity$params[[p]]$all.nonzero.matrix==FALSE) {
    if (step10.expression.specificity$params[[p]]$randomize.mean.expr==FALSE) {
      rmarkdown::render(
        input=step10.expression.specificity$input,
        params=step10.expression.specificity$params[[p]],
        output_file=file.path(
          "..", step10.expression.specificity$output.dir,
          paste("LOC_ELU_DRE.html", sep="")))
      
    } else {
      rmarkdown::render(
        input=step10.expression.specificity$input,
        params=step10.expression.specificity$params[[p]],
        output_file=file.path(
          "..", step10.expression.specificity$output.dir,
          paste("LOC_ELU_DRE", paste0("rnd", step10.expression.specificity$params[[p]]$set.seed, ".html"), sep="_")))       
    }
   
  } else {
    if (step10.expression.specificity$params[[p]]$randomize.mean.expr==FALSE) {
      rmarkdown::render(
        input=step10.expression.specificity$input,
        params=step10.expression.specificity$params[[p]],
        output_file=file.path(
          "..", step10.expression.specificity$output.dir,
          paste("LOC_ELU_DRE", "nonzero.html", sep="_")))
      
    } else {
      rmarkdown::render(
        input=step10.expression.specificity$input,
        params=step10.expression.specificity$params[[p]],
        output_file=file.path(
          "..", step10.expression.specificity$output.dir,
          paste("LOC_ELU_DRE", "nonzero", paste0("rnd", step10.expression.specificity$params[[p]]$set.seed, ".html"), sep="_"))) 
      
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
  dir.create(step11.selection.organ.bias$output.dir, recursive=TRUE)
}

# List parameters
step11.selection.organ.bias$params <- list(
  list(species="LOC", species.name=species.name[["LOC"]], tau.cutoff=0.30),
  list(species="ELU", species.name=species.name[["ELU"]], tau.cutoff=0.30),
  list(species="DRE", species.name=species.name[["DRE"]], tau.cutoff=0.30)
)

##### Run analysis script #####
for (p in (1:length(step11.selection.organ.bias$params))) {
  print(paste("Species:", step11.selection.organ.bias$params[[p]]$species))
  
  rmarkdown::render(
    input=step11.selection.organ.bias$input,
    params=step11.selection.organ.bias$params[[p]],
    output_file=file.path(
      "..", step11.selection.organ.bias$output.dir,
      paste(step11.selection.organ.bias$params[[p]]$species, ".html", sep="")))   
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
  dir.create(step12.organ.bias$output.dir, recursive=TRUE)
}

# List parameters
step12.organ.bias$params <- list(
  list(species="LOC", species.name=species.name[["LOC"]], tau.cutoff=0.30, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE),
  list(species="ELU", species.name=species.name[["ELU"]], tau.cutoff=0.30, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE),
  list(species="DRE", species.name=species.name[["DRE"]], tau.cutoff=0.30, all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE),
  list(species="LOC", species.name=species.name[["LOC"]], tau.cutoff=0.30, all.nonzero.matrix=TRUE, randomized.mean.expr=FALSE),
  list(species="ELU", species.name=species.name[["ELU"]], tau.cutoff=0.30, all.nonzero.matrix=TRUE, randomized.mean.expr=FALSE),
  list(species="DRE", species.name=species.name[["DRE"]], tau.cutoff=0.30, all.nonzero.matrix=TRUE, randomized.mean.expr=FALSE),
  list(species="LOC", species.name=species.name[["LOC"]], tau.cutoff=0.30, all.nonzero.matrix=FALSE, randomized.mean.expr=TRUE, seed.run=12345),
  list(species="ELU", species.name=species.name[["ELU"]], tau.cutoff=0.30, all.nonzero.matrix=FALSE, randomized.mean.expr=TRUE, seed.run=12345),
  list(species="DRE", species.name=species.name[["DRE"]], tau.cutoff=0.30, all.nonzero.matrix=FALSE, randomized.mean.expr=TRUE, seed.run=12345),
  list(species="LOC", species.name=species.name[["LOC"]], tau.cutoff=0.30, all.nonzero.matrix=FALSE, randomized.mean.expr=TRUE, seed.run=67890),
  list(species="ELU", species.name=species.name[["ELU"]], tau.cutoff=0.30, all.nonzero.matrix=FALSE, randomized.mean.expr=TRUE, seed.run=67890),
  list(species="DRE", species.name=species.name[["DRE"]], tau.cutoff=0.30, all.nonzero.matrix=FALSE, randomized.mean.expr=TRUE, seed.run=67890)
)

##### Run analysis script #####
for (p in (1:length(step12.organ.bias$params))) {
  print(paste("Species:", step12.organ.bias$params[[p]]$species))
  
  if (step12.organ.bias$params[[p]]$all.nonzero.matrix==FALSE) {
    if (step12.organ.bias$params[[p]]$randomized.mean.expr==FALSE) {
      rmarkdown::render(
        input=step12.organ.bias$input,
        params=step12.organ.bias$params[[p]],
        output_file=file.path(
          "..", step12.organ.bias$output.dir,
          paste(step12.organ.bias$params[[p]]$species, ".html", sep="")))
      
    } else {
      rmarkdown::render(
        input=step12.organ.bias$input,
        params=step12.organ.bias$params[[p]],
        output_file=file.path(
          "..", step12.organ.bias$output.dir,
          paste(step12.organ.bias$params[[p]]$species,
                paste0("rnd", step12.organ.bias$params[[p]]$seed.run, ".html"), sep="_")))

    }
  } else {
    if (step12.organ.bias$params[[p]]$randomized.mean.expr==FALSE) {
      rmarkdown::render(
        input=step12.organ.bias$input,
        params=step12.organ.bias$params[[p]],
        output_file=file.path(
          "..", step12.organ.bias$output.dir,
          paste(step12.organ.bias$params[[p]]$species, "nonzero.html", sep="_")))
      
    } else {
      rmarkdown::render(
        input=step12.organ.bias$input,
        params=step12.organ.bias$params[[p]],
        output_file=file.path(
          "..", step12.organ.bias$output.dir,
          paste(step12.organ.bias$params[[p]]$species, "nonzero", 
                paste0("rnd", step12.organ.bias$params[[p]]$seed.run, ".html"), sep="_")))
      
    }
  }
}

# Clear workspace
keep <- ls()[grepl("step[0-9]*|species.name",ls())]
rm(list=setdiff(ls(), keep))
gc()

#### Step 13: Combine results from step 8 across species ####

##### Set input and parameters #####
step13.organ.bias.combined <- list()
step13.organ.bias.combined$input <- "analysis/13_Organ_Bias_Combined.Rmd"
step13.organ.bias.combined$output.dir <- "../results/13_Organ_Bias_Combined"

# Create results directory
if (!dir.exists(step13.organ.bias.combined$output.dir)) {
  dir.create(step13.organ.bias.combined$output.dir, recursive=TRUE)
}

# List parameters
step13.organ.bias.combined$params <- list(
  list(all.nonzero.matrix=FALSE, randomized.mean.expr=FALSE),
  list(all.nonzero.matrix=TRUE, randomized.mean.expr=FALSE),
  list(all.nonzero.matrix=FALSE, randomized.mean.expr=TRUE, seed.run=12345),
  list(all.nonzero.matrix=FALSE, randomized.mean.expr=TRUE, seed.run=67890)
)

##### Run analysis script #####
for (p in (1:length(step13.organ.bias.combined$params))) {
  if (step13.organ.bias.combined$params[[p]]$all.nonzero.matrix==FALSE) {
    if (step13.organ.bias.combined$params[[p]]$randomized.mean.expr==FALSE) {
      rmarkdown::render(
        input=step13.organ.bias.combined$input,
        params=step13.organ.bias.combined$params[[p]],
        output_file=file.path(
          "..", step13.organ.bias.combined$output.dir,
          paste("LOC_ELU_DRE", ".html", sep="")))
      
    } else {
      rmarkdown::render(
        input=step13.organ.bias.combined$input,
        params=step13.organ.bias.combined$params[[p]],
        output_file=file.path(
          "..", step13.organ.bias.combined$output.dir,
          paste("LOC_ELU_DRE",
                paste0("rnd", step13.organ.bias.combined$params[[p]]$seed.run, ".html"), sep="_")))
      
    }
  } else {
    if (step13.organ.bias.combined$params[[p]]$randomized.mean.expr==FALSE) {
      rmarkdown::render(
        input=step13.organ.bias.combined$input,
        params=step13.organ.bias.combined$params[[p]],
        output_file=file.path(
          "..", step13.organ.bias.combined$output.dir,
          paste("LOC_ELU_DRE", "nonzero.html", sep="_")))
      
    } else {
      rmarkdown::render(
        input=step13.organ.bias.combined$input,
        params=step13.organ.bias.combined$params[[p]],
        output_file=file.path(
          "..", step13.organ.bias.combined$output.dir,
          paste("LOC_ELU_DRE", "nonzero", 
                paste0("rnd", step13.organ.bias.combined$params[[p]]$seed.run, ".html"), sep="_")))
      
    }
  }
}

# Clear workspace
keep <- ls()[grepl("step[0-9]*|species.name",ls())]
rm(list=setdiff(ls(), keep))
gc()

#### Step 14: Combine results from step 8 for observed and randomized data ####

##### Set input and parameters #####
step14.organ.bias.observed.vs.random <- list()
step14.organ.bias.observed.vs.random$input <- "analysis/14_Organ_Bias_Observed_vs_Random.Rmd"
step14.organ.bias.observed.vs.random$output.dir <- "../results/14_Organ_Bias_Observed_vs_Random"

# Create results directory
if (!dir.exists(step14.organ.bias.observed.vs.random$output.dir)) {
  dir.create(step14.organ.bias.observed.vs.random$output.dir, recursive=TRUE)
}

# List parameters
step14.organ.bias.observed.vs.random$params <- list(
  list(species="LOC", species.name=species.name[["LOC"]], seed.run1=12345, seed.run2=67890, plot.height=16, plot.width=10),
  list(species="ELU", species.name=species.name[["ELU"]], seed.run1=12345, seed.run2=67890, plot.height=18, plot.width=22),
  list(species="DRE", species.name=species.name[["DRE"]], seed.run1=12345, seed.run2=67890, plot.height=18, plot.width=22)
)

##### Run analysis script #####
for (p in (1:length(step14.organ.bias.observed.vs.random$params))) {
  print(paste("Species:", step14.organ.bias.observed.vs.random$params[[p]]$species))
  
  rmarkdown::render(
    input=step14.organ.bias.observed.vs.random$input,
    params=step14.organ.bias.observed.vs.random$params[[p]],
    output_file=file.path(
      "..", step14.organ.bias.observed.vs.random$output.dir,
      paste(step14.organ.bias.observed.vs.random$params[[p]]$species, ".html", sep="")))
}

# Clear workspace
rm(list=setdiff(ls(), ls()[grepl("step[0-9]*|species.name",ls())]))
gc()

#### Save all parameters ####
save.image(file.path("../results", "run_pipeline.Rdata"))