This workflow relies on several essential bioinformatics tools to process genomic data and generate high-quality genotype calls. Each tool is critical for specific steps in the pipeline:

All commands should be run in a terminal:
Linux/macOS: native terminal
Windows: WSL2 (Windows Subsystem for Linux 2) or Linux-compatible terminal

1. WSL2 (Windows Only)
Purpose: This workflow uses Linux-native tools (samtools, bcftools, bwa) on Windows via WSL2.
Enable: Powershell as Admininstrator
1a. Open PowerShell as Administrator
1b. Enable required Windows features:
  > dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart
  > dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart
1c. Restart your computer.
1d. Enable CPU virtualization in BIOS: AMD -> SVM Mode / AMD-V

2. Ubuntu
Purpose: Ubuntu provides a Linux operating system environment required to run the bioinformatics tools used in this workflow (e.g., samtools, bcftools, bwa, GATK). When used through WSL2 on Windows, Ubuntu allows these Linux-native tools to run reliably and consistently, ensuring compatibility with standard genomics software and reproducibility of results across systems.
Download: Microsoft Store
2a. Open Microsoft Store
2b. Search for Ubuntu 20.04 LTS and click Install
2c. Launch Ubuntu after installation
2d Create Linux User:
  Username: any name (e.g., tyler)
  Password: secure password

The following packages can be installed directly from Linux/WSL2 using Ubuntu (below), but if desired to install separately, follow inscructions 3-5.

3. GATK (Genome Analysis Toolkit)
Purpose: GATK is used for variant discovery and genotyping. In this workflow, it is used for HaplotypeCaller, GenotypeGVCFs, and other variant processing steps.
Download: https://gatk.broadinstitute.org
Notes: Java 8 or higher is required. We recommend using the latest stable release.

4. Picard
Purpose: Picard provides tools for manipulating BAM files, including marking duplicates, adding read groups, and fixing headers. These steps are required before running GATK.
Download: https://broadinstitute.github.io/picard/
Notes: Picard runs on Java 8+ and is often included with GATK releases.

5. Samtools
Purpose: Samtools is used for inspecting, indexing, and converting BAM files. It is also used for generating pileups when performing minimal variant calling.
Download: http://www.htslib.org/
Notes: Ensure the installed version supports BAM indexing (.bai) and SAM/BAM conversion commands.

Optional Tools
BCFtools: Often used alongside Samtools for variant calling and filtering. http://www.htslib.org/

---
Installing directly from Linux/WSL2:

```
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
```
