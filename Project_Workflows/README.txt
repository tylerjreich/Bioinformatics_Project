This directory contains the complete set of computational workflows developed for this project, spanning raw genomic data processing through quantitative trait locus (QTL) analysis. Together, these workflows establish a reproducible pipeline for generating genotype data from low-coverage whole-genome sequencing and mapping genetic loci associated with phenotypic traits.

The workflows are organized into two main components:
  1. GATK_Demo – Variant calling and genotype generation
  2. rQTL_Demo – QTL mapping and statistical analysis

Each workflow is fully documented and includes example data and annotated scripts for proof-of-concept validation.

Workflow Structure:

1. GATK_Demo - This workflow demonstrates how raw sequencing data are processed into high-quality genotype calls using the Genome Analysis Toolkit (GATK) and supporting tools.

  Key steps include:
    1a. Read alignment and BAM processing
    1b.  Duplicate marking and read group assignment
    1c. Variant calling in GVCF mode
    1d.  Genotyping of GVCFs into standard VCF files
    1e.  Quality control and inspection of outputs
    
    The GATK workflow is validated using publicly available human sequencing data (NA12878) and is designed to run on Linux or Windows via WSL2.

  Contents:
    - README_Project_Workflow_GATK.txt – Complete step-by-step workflow
    - Tools_Required/ – Installation instructions and resources
    - Example_Data/ – Sample sequencing data and intermediate files
    - Reference_Genomes/ – Reference genome resources and indexing guidance

2. rQTL_Demo - This workflow demonstrates quantitative trait locus (QTL) mapping using the R/qtl package. 
  It provides a fully annotated, end-to-end example of QTL analysis using an F2 intercross dataset.

  Key steps include:
  2a. Importing genotype and phenotype data
  2b. Quality control and segregation checks
  2c. Calculation of genotype probabilities
  2d. Genome-wide QTL scans using Haley–Knott regression
  2e. Permutation testing for significance
  2f. Estimation of QTL confidence intervals
  2g. Visualization of QTL effects

  The workflow is validated using the built-in listeria F2 intercross dataset provided by R/qtl and serves as a proof-of-concept for future analyses.

  Contents:
    - Annotated R/qtl workflow in multiple formats (.Rmd, .html, .pdf, .tex)
    - Figures/ – Individual plots generated during the analysis

Intended Use - These workflows were developed and validated as a proof-of-concept for analyzing F2 intercross data. Once experimental data become available, the same workflows can be applied with minimal modification by updating:
  - Reference genomes
  - Input sequencing files
  - Phenotype measurements

Together, the GATK and R/qtl workflows provide a complete and reproducible framework for linking genomic variation to phenotypic traits.

Notes:
  - All workflows are documented with extensive comments to support learning and reproducibility.
  - Example datasets are used solely for validation and demonstration purposes.
 
