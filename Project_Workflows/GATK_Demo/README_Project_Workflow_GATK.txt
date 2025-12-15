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
> sudo apt update
> sudo apt install openjdk-11-jdk samtools bcftools bwa -y

Verify installations:
> java --version
> samtools --version
> bcftools --version
> bwa 

1. Accessing Project Files

1a. Access repository:
> cd /mnt/c/Bioinformatics_Project

1b. Download data files
Ensure all input files (FASTQ, BAM, reference genome) are in the data/ directory. Maintain the same file structure.

2. Setting Up the GATK Environment
Go to your home directory
> cd ~
Download GATK
> wget https://github.com/broadinstitute/gatk/releases/download/4.6.0.0/gatk-4.6.0.0.zip
Unzip GATK
> unzip gatk-4.6.0.0.zip
Verify
> gatk --version

3. Running the Workflow (Sample Data)
Sample Data:https://console.cloud.google.com/storage/browser/genomics-public-data/test-data/dna/wgs/hiseq2500/NA12878?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))
Sample Reference Genome: https://console.cloud.google.com/storage/browser/genomics-public-data/references/Homo_sapiens_assembly19_1000genomes_decoy?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))

3a. Sort BAM
> java -jar /mnt/c/bioinformatics_tools/gatk/gatk-package-4.6.2.0-local.jar SortSam \
    -I /mnt/c/Bioinformatics_Project/Data/test-data_dna_wgs_hiseq2500_NA12878_H06HDADXX130110.1.ATCACGAT.200k_reads.bam \
    -O /mnt/c/Bioinformatics_Project/Data/test_sorted.bam \
    -SO coordinate
Output: test_sorted.bam
Description: This step sorts the input BAM file by genomic coordinates (chromosome and position). Coordinate sorting is a required preprocessing step for most downstream analyses, including duplicate marking, indexing, and variant calling with tools such as GATK HaplotypeCaller and samtools.The input BAM file contains aligned sequencing reads, but these reads may be ordered by read name or in an arbitrary order depending on how the aligner wrote the file. Sorting ensures that reads are ordered first by reference contig (e.g., chr1, chr2, etc.) and then by increasing genomic position within each contig.
Why this step is important: 
- Required for indexing: BAM files must be coordinate-sorted before a .bai index can be created.
- Required for duplicate marking: Tools like MarkDuplicates assume coordinate-sorted input.
- Improves performance: Many downstream tools rely on sorted reads to efficiently access genomic regions.
- Ensures compatibility: GATK Best Practices pipelines require coordinate-sorted BAM files.
Parameter details:
-I — Input BAM file containing aligned reads.
-O — Output BAM file that will be sorted by coordinates.
-SO coordinate — Specifies coordinate-based sorting (as opposed to query-name sorting).
After this step, test_sorted.bam is ready for downstream processing steps such as duplicate marking, read group assignment, indexing, and variant calling.

3b. Mark Duplicates
> java -jar C:\bioinformatics_tools\gatk\gatk-package-4.6.2.0-local.jar MarkDuplicates \
    -I C:\Bioinformatics_Project\Data\test_sorted.bam \
    -O C:\Bioinformatics_Project\Data\test_dedup.bam \
    -M C:\Bioinformatics_Project\Data\test_metrics.txt
Outputs: test_dedup.bam, test_metrics.txt
Description: This step identifies and marks duplicate reads in the coordinate-sorted BAM file. Duplicate reads most commonly arise during PCR amplification in library preparation, where multiple sequencing reads originate from the same original DNA fragment rather than representing independent observations. MarkDuplicates does not remove these reads by default. Instead, it adds a duplicate flag in the BAM file, allowing downstream tools (such as GATK HaplotypeCaller) to ignore or appropriately downweight these reads during variant calling.
Why this step is important:
- Prevents bias in variant calling: Duplicate reads can artificially inflate read depth and lead to false-positive variant calls.
- Required by GATK Best Practices: GATK tools expect duplicate marking prior to variant discovery.
- Preserves data: Marking (rather than removing) duplicates allows flexibility for downstream analyses and quality control.
How duplicates are identified:
- For paired-end reads, duplicates are detected based on identical start positions of read pairs.
- For single-end reads, duplicates are identified based on identical alignment start positions.
- When multiple duplicates are detected, the read (or read pair) with the highest sum of base quality scores is retained as the representative, and others are flagged as duplicates.
Output details:
- test_dedup.bam — The BAM file with duplicate reads marked using the SAM flag 0x400.
- test_metrics.txt — A summary report containing duplication statistics, including:
	Total reads examined
	Number and percentage of duplicate reads
	Optical vs. PCR duplicates
	Library-level duplication metrics
After this step, the BAM file is suitable for read group assignment, indexing, and variant calling with GATK.

3c. Index BAM
> java -jar /mnt/c/bioinformatics_tools/gatk/gatk-package-4.6.2.0-local.jar BuildBamIndex \
    -I /mnt/c/Bioinformatics_Project/Data/test_dedup.bam
Output: test_dedup.bai
Description: This step creates a BAM index file (.bai) for the duplicate-marked BAM file. The index enables rapid, random access to specific genomic regions within the BAM file, allowing downstream tools to efficiently retrieve reads aligned to particular chromosomes or coordinates without scanning the entire file.
Why this step is important:
- Required for variant calling: GATK tools, including HaplotypeCaller, require an indexed BAM file to operate.
- Improves performance: Indexing dramatically speeds up analyses by allowing tools to access only the genomic regions of interest.
- Enables compatibility: Many downstream tools (e.g., GATK, samtools, IGV) rely on the presence of a .bai file for proper BAM file handling.
How it works:
- The index records the file offsets corresponding to genomic coordinate ranges.
- Each contig (chromosome) in the BAM file is indexed separately.
- The BAM file must be coordinate-sorted prior to indexing, which is satisfied by the previous sorting step.
After this step, test_dedup.bam is fully prepared for read group assignment and variant calling with GATK.

3d. Index Reference Genome for Alignment (BWA)
Before aligning reads, the reference genome must be indexed with BWA. This step only needs to be done once per reference genome.
> bwa index /mnt/c/Bioinformatics_Project/Data/Homo_sapiens_assembly38.fasta
Description: This step builds the BWA index for the reference genome, which is required before aligning sequencing reads with bwa mem. BWA indexing preprocesses the reference FASTA and generates several auxiliary index files that allow the aligner to rapidly search for matching sequences using the Burrows–Wheeler Transform (BWT).
Why this step is important:
- Mandatory for alignment: bwa mem cannot run without these index files.
- Performance optimization: Indexing enables fast and memory-efficient alignment of millions of short reads.
- One-time operation: The index only needs to be generated once per reference genome and can be reused for all future alignments against the same FASTA.- Creates all the necessary BWA index files (.amb, .ann, .bwt, .pac, .sa) in the same directory as the FASTA.
What this step produces: The following index files are created in the same directory as the reference FASTA:
- .amb — ambiguous bases in the reference
- .ann — reference sequence annotations
- .bwt — Burrows–Wheeler transformed reference
- .pac — packed reference sequence
- .sa — suffix array for fast lookups
Approximate runtime:
- Human genome (~3.2 Gb): ~50–60 minutes on a standard desktop system
- Runtime varies with CPU speed, number of cores, disk performance, and available memory.
Dependencies: BWA must be installed (on WSL2/Linux):
> sudo apt install bwa
Example Output:
After successful indexing of the sample reference genome, you should see the following files in the same directory as your FASTA:
Homo_sapiens_assembly38.fasta.amb
Homo_sapiens_assembly38.fasta.ann
Homo_sapiens_assembly38.fasta.bwt
Homo_sapiens_assembly38.fasta.pac
Homo_sapiens_assembly38.fasta.sa

4. Running GATK HaplotypeCaller

4a. Access Windows Files from WSL2
Navigate to your working directory:
> cd /mnt/c/Bioinformatics_Project/Data

4b. Prepare BAM for GATK Haplotype Caller
Index BAM
> samtools index test_dedup_rg.bam
Description: This step creates an index file (.bai) for the coordinate-sorted BAM file that includes read group information. The BAM index allows GATK tools, including HaplotypeCaller, to efficiently access reads at specific genomic regions without scanning the entire file.
Why this step is important:
- Required by GATK: HaplotypeCaller will not run on an unindexed BAM file.
- Efficient random access: Indexing enables rapid retrieval of reads overlapping genomic intervals.
- Quality control: Ensures the BAM file is properly sorted and formatted for downstream variant calling.
Input: test_dedup_rg.bam — duplicate-marked, coordinate-sorted BAM with read group information
Output: test_dedup_rg.bam.bai — BAM index file
Notes:
- The BAM file must be sorted by coordinate before indexing.
- This step must be rerun if the BAM file is modified or regenerated.
- After indexing, the BAM file is fully prepared for variant calling with GATK HaplotypeCaller.
 - BAM files must be indexed after any modification. Since read group information is added in this step, the resulting BAM must be re-indexed before use with GATK.

4c. Run GATK HaplotypeCaller (GVCF)
> gatk --java-options "-Xmx16g" HaplotypeCaller \
  -R Homo_sapiens_assembly38.fasta \
  -I test_dedup_rg.bam \
  -O test_raw_variants.g.vcf.gz \
  -ERC GVCF \
  --sample-name NA12878
Description: This step performs variant calling using GATK’s HaplotypeCaller, which identifies SNPs and indels by locally re-assembling reads into haplotypes rather than relying solely on pileups. Running HaplotypeCaller in GVCF mode (-ERC GVCF) produces a genomic VCF (gVCF) that contains variant calls and reference confidence information for non-variant sites. This format is required for scalable joint genotyping across multiple samples. The input BAM (test_dedup_rg.bam) must be coordinate-sorted, indexed, and contain valid read group information. The reference genome must match the one used for alignment.
Notes:
Key options explained:
-R – Reference genome FASTA
-I – Input BAM file with read groups
-O – Output gVCF file (compressed)
-ERC GVCF – Emit reference confidence blocks for non-variant sites
--sample-name – Explicitly sets the sample name (useful if the BAM header is missing or contains multiple samples)
--java-options "-Xmx16g" – Allocates up to 16 GB of memory for the Java process
Check output
> zcat test_raw_variants.g.vcf.gz | head -20
This confirms that the gVCF file was generated correctly and contains the expected VCF header and records.

4d. Genotype GVCF (VCF)
> gatk --java-options "-Xmx16g" GenotypeGVCFs \
  -R Homo_sapiens_assembly38.fasta \
  -V test_raw_variants.g.vcf.gz \
  -O test_variants.vcf.gz
Description: This step converts the intermediate gVCF produced by HaplotypeCaller into a standard, fully genotyped VCF file. GenotypeGVCFs evaluates the reference confidence information stored in the gVCF and assigns final genotypes (e.g., 0/0, 0/1, 1/1) at each variant site. For a single-sample workflow, this step finalizes variant calls for downstream analyses such as variant filtering, annotation, and summary statistics. In multi-sample workflows, GenotypeGVCFs can be used to jointly genotype multiple gVCFs, ensuring consistent genotyping across samples.
Key options explained:
-R – Reference genome FASTA (must match alignment and HaplotypeCaller steps)
-V – Input gVCF file
-O – Output compressed VCF file containing finalized variant calls
--java-options "-Xmx16g" – Allocates memory for genotyping (can be reduced for small datasets)
Notes:
- The output VCF contains only variant sites (unlike the gVCF, which includes reference blocks).
- The resulting VCF is suitable for downstream tools such as bcftools, vcftools, and variant annotation pipelines.
Check output:
> bcftools view test_variants.vcf.gz | head
Displays the first few header lines and variant records from the final VCF file to quickly verify that variant calling completed successfully and the file is readable.

5. Inspect VCF Summary

Use bcftools or similar tools to explore basic statistics for quality checks:
> bcftools stats /mnt/c/Bioinformatics_Project/Data/test_variants.vcf.gz > variant_stats.txt
Description: Description:
This step generates a summary of the variants in the final VCF file using bcftools stats. The output includes counts of SNPs and indels, transitions vs. transversions, allele frequencies, depth of coverage, and other basic quality metrics.
These statistics are useful for:
- Performing quick quality control on variant calls
- Verifying that the variant calling pipeline worked correctly
- Identifying potential issues before downstream analyses, such as annotation or filtering
Optional: To quickly view the first few lines of the summary:
> head variant_stats.txt

#####################################################

Preparing Variant Data for R/qtl Analysis

Once you have generated high-quality genotype likelihoods and VCF files using this GATK workflow, 
the next step is to prepare the data for quantitative trait loci (QTL) mapping. A detailed R/qtl workflow for converting 
VCFs to CSV, adding phenotypes, performing quality control, calculating genotype probabilities, and conducting genome 
scans is provided in the R Markdown file: "rQTL_Demo.Rmd".

This R Markdown file contains step-by-step instructions and annotated code for performing QTL mapping using the 
listeria F2 intercross dataset as a proof-of-concept.



