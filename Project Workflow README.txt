Project Workflow README

Overview:
This repository contains a workflow for processing genomic data 
using GATK (Genome Analysis Toolkit). This README provides step-by-step 
instructions for accessing the project, setting up the environment, and running analyses.

First, open your terminal (Linux/macOS) or WSL/Command Prompt (Windows with Linux environment).

1. Accessing Project Files

1a. Clone the repository (if using version control):
git clone <repository_url>
cd <repository_name> (cd C:\Bioinformatics_Project)

1b. Download data files (if not included in the repository):
Ensure all input files (e.g., FASTQ, BAM, reference genome) are stored in the data/ directory.
Maintain the same file structure to avoid breaking the workflow.

1c. Verify file integrity (optional but recommended):
CertUtil -hashfile <file_name> MD5

2. Setting Up the Environment 

2a. Load required modules (example for a cluster environment):
module load java/1.8.0
module load gatk/4.3.0.0

2b. Verify GATK installation:
java -jar C:\bioinformatics_tools\gatk\gatk-package-4.6.2.0-local.jar
You should see the version of GATK printed to confirm it’s loaded correctly.

3. Running the Workflow (Testing with Sample Data)
from https://console.cloud.google.com/storage/browser/_details/genomics-public-data/test-data/dna/wgs/hiseq2500/NA12878/H06HDADXX130110.1.ATCACGAT.200k_reads.bam?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22)))

All following commands should be run in the terminal or on the cluster, depending on where GATK is loaded.

3a. Sorting BAM files
java -jar C:\bioinformatics_tools\gatk\gatk-package-4.6.2.0-local.jar SortSam \
    -I C:\Bioinformatics_Project\Data\test-data_dna_wgs_hiseq2500_NA12878_H06HDADXX130110.1.ATCACGAT.200k_reads.bam \
    -O C:\Bioinformatics_Project\Data\test_sorted.bam \
    -SO coordinate
Description: Sorts the BAM file by genomic coordinates.
Output: test_sorted.bam

3b. Marking duplicates
java -jar C:\bioinformatics_tools\gatk\gatk-package-4.6.2.0-local.jar MarkDuplicates \
    -I C:\Bioinformatics_Project\Data\test_sorted.bam \
    -O C:\Bioinformatics_Project\Data\test_dedup.bam \
    -M C:\Bioinformatics_Project\Data\test_metrics.txt
Description: Marks duplicate reads that can arise from PCR amplification.
Outputs: test_dedup.bam and test_metrics.txt

3c. Indexing BAM files
java -jar C:\bioinformatics_tools\gatk\gatk-package-4.6.2.0-local.jar BuildBamIndex \
    -I C:\Bioinformatics_Project\Data\test_dedup.bam
Description: Creates an index file for the BAM to enable rapid random access.
Output: test_dedup.bai


