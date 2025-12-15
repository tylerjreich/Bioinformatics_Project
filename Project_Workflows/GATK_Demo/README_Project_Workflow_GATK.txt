Project Workflow README
Tyler J. Reich

Preface on this Project:

With limited background in computer science, learning to work with genomic data processing pipelines has been a significant challenge, 
but one I am eager to tackle. I am grateful that we are allowed to use generative AI in this project, as ChatGPT has been an important 
resource in helping me understand how to set up and run the GATK workflow. While ChatGPT assisted me in writing scripts and commands for 
aligning reads, calling variants, and generating VCF files, I carefully reviewed and tested each step to ensure I understand how and why 
it works. I am prepared to explain this workflow to demonstrate my understanding of the computational processes involved in genomic data analysis.

Overview:

This repository contains a workflow for processing genomic data using GATK (Genome Analysis Toolkit).
This README provides step-by-step instructions for accessing project files, setting up the environment, and running analyses.

All commands should be run in a terminal:
- Linux/macOS: native terminal
- Windows: WSL2 (Windows Subsystem for Linux 2) or Linux-compatible terminal

0. Setting Up WSL2 (Windows Only)
This workflow uses Linux-native tools (samtools, bcftools, bwa) on Windows via WSL2.

Step 0a: Enable WSL2 and Virtualization
1.  Open PowerShell as Administrator.
2. Enable required Windows features:
> dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart
> dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart
3. Restart your computer.
4. Enable CPU virtualization in BIOS:
	- Intel → Intel VT-x / Intel Virtualization Technology
	- AMD → SVM Mode / AMD-V

Step 0b: Install Ubuntu
1. Open Microsoft Store.
2. Search for Ubuntu 20.04 LTS and click Install.
3. Launch Ubuntu after installation.

Step 0c: Create Linux User
- Username: any name (e.g., tyler)
- Password: secure password

Step 0d: Update Ubuntu Packages
> sudo apt update
> sudo apt upgrade -y

Step 0e: Install Required Tools
> sudo apt install samtools bcftools bwa -y

Verify installations:
> samtools --version
> bcftools --version
> bwa

1. Accessing Project Files

1. Access respository
> cd /mnt/c/Bioinformatics_Project

1b. Download data files
Ensure all input files (FASTQ, BAM, reference genome) are in the data/ directory. Maintain the same file structure.

1c. Verify file integrity (optional)
> CertUtil -hashfile <file_name> MD5

2. Setting Up the GATK Environment

2a. Load modules (cluster example)
> module load java/1.8.0
> module load gatk/4.3.0.0

2b. Verify GATK installation
> java -jar C:\bioinformatics_tools\gatk\gatk-package-4.6.2.0-local.jar
You should see the GATK version printed.

3. Running the Workflow (Sample Data)
Sample Data:https://console.cloud.google.com/storage/browser/genomics-public-data/test-data/dna/wgs/hiseq2500/NA12878?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))
Sample Reference Genome: https://console.cloud.google.com/storage/browser/genomics-public-data/references/Homo_sapiens_assembly19_1000genomes_decoy?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))

3a. Sort BAM
> java -jar C:\bioinformatics_tools\gatk\gatk-package-4.6.2.0-local.jar SortSam \
    -I C:\Bioinformatics_Project\Data\test-data_dna_wgs_hiseq2500_NA12878_H06HDADXX130110.1.ATCACGAT.200k_reads.bam \
    -O C:\Bioinformatics_Project\Data\test_sorted.bam \
    -SO coordinate
Output: test_sorted.bam
Description: Sorts BAM by genomic coordinates

3b. Mark Duplicates
> java -jar C:\bioinformatics_tools\gatk\gatk-package-4.6.2.0-local.jar MarkDuplicates \
    -I C:\Bioinformatics_Project\Data\test_sorted.bam \
    -O C:\Bioinformatics_Project\Data\test_dedup.bam \
    -M C:\Bioinformatics_Project\Data\test_metrics.txt
Outputs: test_dedup.bam, test_metrics.txt
Description: Marks duplicate reads from PCR

3c. Index BAM
> java -jar C:\bioinformatics_tools\gatk\gatk-package-4.6.2.0-local.jar BuildBamIndex \
    -I C:\Bioinformatics_Project\Data\test_dedup.bam
Output: test_dedup.bai

3d. Index Reference Genome for Alignment (BWA)
Before aligning reads, the reference genome must be indexed with BWA. This step only needs to be done once per reference genome.
> bwa index /mnt/c/Bioinformatics_Project/Data/Homo_sapiens_assembly38.fasta
Description:
Creates all the necessary BWA index files (.amb, .ann, .bwt, .pac, .sa) in the same directory as the FASTA.
Required for efficient read alignment using bwa mem.
Approximate Runtime:
Human genome (~3.2 Gb): 50–60 minutes on a standard desktop machine.
Runtime will vary depending on CPU speed, number of cores, and available memory.
Dependencies:
BWA must be installed (on WSL2/Linux):
> sudo apt install bwa
Example Output:
After successful indexing, you should see the following files in the same directory as your FASTA:
Homo_sapiens_assembly38.fasta.amb
Homo_sapiens_assembly38.fasta.ann
Homo_sapiens_assembly38.fasta.bwt
Homo_sapiens_assembly38.fasta.pac
Homo_sapiens_assembly38.fasta.sa

4. Run GATK HaplotypeCaller

4a. Access Windows Files from WSL2
Windows drives are mounted under /mnt/.
Example:C:\Users\Tyler Reich\data → /mnt/c/Users/Tyler\ Reich/data
Navigate to your working directory:
> cd /mnt/c/Users/Tyler\ Reich/data

4b. Minimal Variant Calling Workflow (Samtools + BCFtools)
> samtools mpileup -f reference.fa sample.bam | bcftools call -mv -Ov -o variants.vcf
Description:
	reference.fa → reference genome
	sample.bam → aligned reads
	variants.vcf → output variants file
	-mv → call SNPs and indels
	-Ov → output as VCF

4c. Prepare BAM for GATK
Step 1: Check BAM header for contigs
> samtools view -H test_dedup.bam | grep "@SQ"
Step 2: Extract reference contigs
> samtools faidx Homo_sapiens_assembly38.fasta
> cut -f1,2 Homo_sapiens_assembly38.fasta.fai > ref_contigs.txt
Step 3: Fix BAM header
> samtools reheader <(cat ref_contigs.txt) test_dedup.bam > test_dedup_fixed.bam
Step 4: Add read group information
> samtools addreplacerg \
  -r "@RG\tID:NA12878\tSM:NA12878\tLB:lib1\tPL:illumina\tPU:unit1" \
  -o test_dedup_rg.bam \
  test_dedup_fixed.bam
Step 5: Index BAM
> samtools index test_dedup_rg.bam

4d. Run GATK HaplotypeCaller (GVCF)
> gatk --java-options "-Xmx4g" HaplotypeCaller \
  -R Homo_sapiens_assembly38.fasta \
  -I test_dedup_rg.bam \
  -O test_raw_variants.g.vcf.gz \
  -ERC GVCF \
  --sample-name NA12878
Notes:
-ERC GVCF → emits reference confidence blocks
--sample-name → required if BAM contains multiple samples or missing @RG
Check output
> zcat test_raw_variants.g.vcf.gz | head -20

4e. Genotype GVCF (VCF)
> gatk --java-options "-Xmx16g" GenotypeGVCFs \
  -R Homo_sapiens_assembly38.fasta \
  -V test_raw_variants.g.vcf.gz \
  -O test_variants.vcf.gz
Produces fully processed, single-sample VCF for downstream analysis

5. Generate Genotypes from GVCF

After running HaplotypeCaller in GVCF mode, use GenotypeGVCFs to produce a standard VCF with called genotypes:
> gatk --java-options "-Xmx16g" GenotypeGVCFs \
  -R /mnt/c/Bioinformatics_Project/Data/Homo_sapiens_assembly38.fasta \
  -V /mnt/c/Bioinformatics_Project/Data/test_raw_variants.g.vcf.gz \
  -O /mnt/c/Bioinformatics_Project/Data/test_variants.vcf.gz
Description:
- Converts the per-sample GVCF into a fully genotyped VCF.
- Output (test_variants.vcf.gz) contains SNPs, indels, and genotype likelihoods.
- Includes all standard INFO and FORMAT fields necessary for downstream analysis.
Check output:
> zcat /mnt/c/Bioinformatics_Project/Data/test_variants.vcf.gz | head -20
This will display the header and first few variants, confirming that genotypes have been called.

6. Inspect VCF Summary

Use bcftools or similar tools to explore basic statistics for quality checks:
> bcftools stats /mnt/c/Bioinformatics_Project/Data/test_variants.vcf.gz > variant_stats.txt
Description:
- Produces counts of SNPs, indels, transitions/transversions, depth metrics, etc.
- Useful for confirming the pipeline is working and generating high-quality genotype calls.

#####################################################

Preparing Variant Data for R/qtl Analysis

Once you have generated high-quality genotype likelihoods and VCF files using this GATK workflow, 
the next step is to prepare the data for quantitative trait loci (QTL) mapping. A detailed R/qtl workflow for converting 
VCFs to CSV, adding phenotypes, performing quality control, calculating genotype probabilities, and conducting genome 
scans is provided in the R Markdown file: "rQTL_Demo.Rmd".

This R Markdown file contains step-by-step instructions and annotated code for performing QTL mapping using the 
listeria F2 intercross dataset as a proof-of-concept.


