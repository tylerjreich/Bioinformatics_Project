This directory contains a fully annotated R/qtl workflow demonstrating quantitative trait loci (QTL) mapping using the listeria F2 intercross dataset, which is included with the R/qtl package. The purpose of this demo is to validate the downstream analysis portion of the project and to provide a reproducible, well-documented example of QTL mapping.

The workflow covers data inspection, quality control, genotype probability calculation, genome-wide QTL scans, permutation-based significance testing, confidence interval estimation, and visualization of marker effects.

File Descriptions:

1. rQTL_Demo.Rmd - This is the primary source file for the R/qtl analysis. 
  It contains:
    Fully annotated R code
    Explanatory text describing each analytical step 
    A complete QTL mapping workflow using the listeria dataset
    Figure generation 
  All other file formats in this directory are generated from this R Markdown file.

2. rQTL_Demo.html - An HTML-rendered version of the workflow, suitable for:
  Quick viewing in a web browser
  Sharing results without requiring R or LaTeX
  Reviewing figures

3. rQTL_Demo.pdf - A PDF-rendered version of the workflow, intended for:
  Formal review and submission
  Archival documentation
  Offline viewing

  This version presents the same content as the R Markdown file in a static, publication-style format.

4. rQTL_Demo.tex - The LaTeX source file generated during PDF rendering. 
  This file is included for:
    Transparency and reproducibility
    Troubleshooting LaTeX compilation issues
    Advanced customization of document formatting (if needed)

5. Figures/ - This folder contains individual PDF files for each figure generated during the R/qtl workflow, including:
  Genome-wide LOD score plots
  Chromosome-specific LOD profiles
  Confidence interval visualizations
  Marker effect plots
  
  Separating figures into their own directory allows for easy reuse in manuscripts, presentations, or reports.

Dataset Used:
All analyses in this demo use the listeria dataset, an example F2 intercross provided with the R/qtl package (Broman et al., 2003). 
This dataset is used strictly as a proof-of-concept to verify that the workflow functions correctly for F2 populations.

Intended Use: - This R/qtl demo workflow is designed to:
  Validate the QTL mapping pipeline
  Demonstrate correct data handling and analysis in R/qtl
  Provide a reusable template for future QTL studies
  Serve as the downstream complement to the GATK variant calling workflow

  Once experimental data are available, this workflow can be applied directly by replacing:
    The input genotype CSV file
    The phenotype data
    Marker names and chromosome identifiers
  No changes to the analytical logic are required.

Relationship to GATK Workflow:
  The R/qtl workflow is intended to be used after variant calling and genotype generation using the GATK pipeline provided in the GATK_Demo/ directory. Together, these workflows form a complete pipeline from raw sequencing data to QTL discovery.

Reference:

Broman, K. W., Wu, H., Sen, Ś., & Churchill, G. A. (2003). R/qtl: QTL mapping in experimental crosses. Bioinformatics, 19(7), 889–890. https://doi.org/10.1093/bioinformatics/btg112
