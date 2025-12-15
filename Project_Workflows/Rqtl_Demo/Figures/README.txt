This directory contains individual figure files generated during the R/qtl demonstration workflow using the listeria F2 intercross dataset. These figures visualize key steps in the QTL mapping process, including data quality assessment, genome-wide scans, significance testing, confidence interval estimation, and genotype–phenotype relationships.

All figures are provided as PDF files to facilitate reuse in reports, presentations, and manuscripts.

Figure Descriptions:

1. Missing_Genotypes.pdf
  Displays the pattern of missing genotype data across individuals and markers. 
  This plot is used as a quality control check to identify markers or samples with excessive missing data that could bias QTL analyses.

2. QTL_Scan.pdf
  Shows the genome-wide LOD score profile from the Haley–Knott regression. 
  This plot summarizes evidence for QTL across all chromosomes prior to applying statistical significance thresholds.

3. QTL_Scan_Significant_Threshold.pdf
  The genome-wide QTL scan with the permutation-derived 5% significance threshold overlaid. 
  This figure highlights genomic regions where LOD scores exceed the threshold, indicating statistically significant QTL.

4. QTL_Peak_Chromosome5.pdf
  A chromosome-specific LOD profile for chromosome 5. 
  This figure illustrates the location of a significant QTL peak, the corresponding confidence interval, and the marker with the maximum LOD score.

5. QTL_Peak_Chromosome13.pdf
  A chromosome-specific LOD profile for chromosome 13, highlighting a second significant QTL. 
  The plot includes confidence interval boundaries and the peak marker associated with the strongest signal.

6. Genotype_Effect_Marker.pdf
  Displays the effect of genotype classes at a peak marker on the phenotype. 
  This plot visualizes how phenotypic values differ among genotypes, providing insight into the direction and magnitude of the QTL effect.

Notes:
  All figures were generated directly from the annotated R Markdown workflow (rQTL_Demo.Rmd).
  Marker names, chromosome numbers, and trait values are specific to the listeria dataset and are used here strictly for proof-of-concept.
  These figures demonstrate the expected outputs of the workflow and serve as templates for future analyses using experimental data.
