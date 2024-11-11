# Inter-individual gene expression variability implies stable regulation of brain-biased genes across organs

Christabel F. Bucao<sup>1,2*</sup>, Consolée Aletti<sup>1</sup>, Sébastien Moretti<sup>1,2</sup>, Andrew W. Thompson<sup>3,4,5</sup>, Brett L. Racicot<sup>4</sup>, Catherine A. Wilson<sup>6</sup>, Julien Bobe<sup>7</sup>, Ingo Braasch<sup>4,5</sup>, Yann Guiguen<sup>7</sup>, John H. Postlethwait<sup>6</sup>, Marc Robinson-Rechavi<sup>1,2</sup>

<sup>1</sup> Department of Ecology and Evolution, University of Lausanne, Lausanne, Switzerland
\
<sup>2</sup> SIB Swiss Institute of Bioinformatics, Lausanne, Switzerland
\
<sup>3</sup> Department of Biological Sciences, Western Michigan University, Kalamazoo, Michigan, USA
\
<sup>4</sup> ﻿Department of Integrative Biology, Michigan State University, East Lansing, Michigan, USA
\
<sup>5</sup> Ecology, Evolution, and Behavior Program, Michigan State University, East Lansing, Michigan, USA
\
<sup>6</sup> ﻿Institute of Neuroscience, University of Oregon, Eugene, Oregon, USA
\
<sup>7</sup> INRAE, LPGP, Rennes 35000, France

<sup>*</sup> Corresponding author: christabelfloi.bucao@unil.ch

## Abstract

Phenotypic variation among individuals plays a key role in evolution, since variation provides the material on which natural selection can act. One important link between genetic and phenotypic variation is gene expression. As for other phenotypes, the range of accessible expression variation is limited and biased by different evolutionary and developmental constraints. Gene expression variability broadly refers to the tendency of a gene to vary in expression (i.e., between individuals or cells) due to stochastic fluctuations or differences in genetic, epigenetic, or environmental factors, separately from the differences between e.g. organs. Variability due to biomolecular stochasticity (transcriptional ‘noise’) and cell-to-cell heterogeneity has been well-studied in isogenic populations of unicellular organisms such as bacteria and yeasts. However, for more complex organisms with multiple cells, tissues, and organs sharing the same genetic background, the interplay between inter-individual expression variability, gene and organ function, and gene regulation remains an open question. In this study, we used highly multiplexed 3’-end Bulk RNA Barcoding and sequencing (BRB-seq) to generate transcriptome profiles spanning at least nine organs in outbred individuals of three ray-finned fish species: zebrafish, Northern pike, and spotted gar. For each condition, we measured expression variation per gene independent of mean expression level. We observed that lowly variable genes are enriched in cellular housekeeping functions whereas highly variable genes are enriched in stimulus-response functions. Furthermore, genes with highly variable expression between individuals evolve under weaker purifying selection at the coding sequence level, indicating that intra-species gene expression variability predicts inter-species protein sequence divergence. Genes that are broadly expressed across organs tend to be both highly expressed and lowly variable between individuals, whereas organ-biased genes are typically highly variable within their top organ of expression. For genes with organ-biased expression profiles, we inferred differences in selective pressure on gene regulation depending on their top organ. We found that genes with peak expression in the brain have low inter-individual expression variability across non-nervous organs, suggesting stabilizing selection on regulatory evolution of brain-biased genes. Conversely, liver-biased genes have highly variable expression across organs, implying weaker regulatory constraints. These patterns show that gene regulatory mechanisms evolved differently based on constraints on the primary organ.

## Directory Structure
- `config/`: Contains YAML file indicating package versions for conda environment
  
- `data/`: Contains input data
  -  `counts/`: Contains counts and UMI-deduplicated counts. Currently under embargo and will be made available upon acceptance for publication.
  -  `gene_metadata/`: Contains gene biotype information from Ensembl
  -  `sample_metadata/`: Contains sample metadata files for each species
  -  `selectome/`: Contains selection statistics from the [Selectome](https://selectome.org/) database
    
- `results/`: Contains output files sorted by subfolders labeled after each step of the analysis pipeline. Only R notebook HTML files are available on the Git repository, please check Zenodo for R data files.
  - `run_pipeline.Rdata`: Contains all parameters used for each step of the analysis pipeline
    
- `workflow/`: Contains scripts used for the analysis pipeline
  - `analysis/`: Contains all steps of the analysis pipeline, available as .Rmd files
  - `functions/`: Contains all functions used for analysis/
  - `renv/`: Used for package management in R
  - `run_pipeline.R`: Runs all the steps under analysis/
  - `run_go_figure.sh`: Runs [GO-Figure!](https://gitlab.com/evogenlab/GO-Figure) 1.0.0 (downloaded separately)
  - `demultiplex_brbseq_fastq.sh`: Used for demultiplexing BRB-seq fastq files using [BRB-seqTools](https://github.com/DeplanckeLab/BRB-seqTools) 1.6.1 (downloaded separately) for uploading to NCBI SRA
  - `rename_fastq_files.sh`: Used for renaming demultiplexed fastq files by mapping each barcode to their corresponding sample name
  - `renv.lock`: Lockfile for managing R package versions

## Species Codes
- **LOC**: *Lepisosteus oculatus* (spotted gar)
- **ELU**: *Esox lucius* (Northern pike)
- **DRE**: *Danio rerio* (zebrafish)

