This workflow relies on several essential bioinformatics tools to process genomic data and generate high-quality genotype calls. Each tool is critical for specific steps in the pipeline:

1. GATK (Genome Analysis Toolkit)
Purpose: GATK is used for variant discovery and genotyping. In this workflow, it is used for HaplotypeCaller, GenotypeGVCFs, and other variant processing steps.
Download: https://gatk.broadinstitute.org
Notes: Java 8 or higher is required. We recommend using the latest stable release.

2. Picard
Purpose: Picard provides tools for manipulating BAM files, including marking duplicates, adding read groups, and fixing headers. These steps are required before running GATK.
Download: https://broadinstitute.github.io/picard/
Notes: Picard runs on Java 8+ and is often included with GATK releases.

3. Samtools
Purpose: Samtools is used for inspecting, indexing, and converting BAM files. It is also used for generating pileups when performing minimal variant calling.
Download: http://www.htslib.org/
Notes: Ensure the installed version supports BAM indexing (.bai) and SAM/BAM conversion commands.

Optional Tools
BCFtools: Often used alongside Samtools for variant calling and filtering. http://www.htslib.org/

Installing directly from Linux/WSL2:

# Update package lists
sudo apt update
sudo apt upgrade -y

# Install Samtools and BCFtools
sudo apt install samtools bcftools -y

# Install Java (required for GATK and Picard)
sudo apt install openjdk-11-jdk -y

# Verify Java installation
java -version

# Install GATK (latest stable release)
# Download GATK tar.gz from the official website
wget https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip
unzip gatk-4.3.0.0.zip
cd gatk-4.3.0.0

# Add GATK to PATH (optional, for convenience)
echo 'export PATH=$PATH:/path/to/gatk-4.3.0.0' >> ~/.bashrc
source ~/.bashrc

# Verify GATK installation
gatk --version

# Picard can be downloaded separately or is included in GATK
wget https://github.com/broadinstitute/picard/releases/download/2.27.4/picard.jar
