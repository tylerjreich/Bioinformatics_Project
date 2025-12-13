The workflow requires a reference genome for alignment and variant calling:

Reference genome:
https://console.cloud.google.com/storage/browser/_details/genomics-public-data/references/Homo_sapiens_assembly19_1000genomes_decoy/Homo_sapiens_assembly19_1000genomes_decoy.fasta;tab=live_object

Instructions:
1. Download the FASTA file and associated index/dictionary files (or generate them locally using samtools faidx and Picard CreateSequenceDictionary).
2. Place all files in your project reference_genomes/ directory, e.g.:
C:\Bioinformatics_Project\reference_genomes\Homo_sapiens_assembly19.fasta
C:\Bioinformatics_Project\reference_genomes\Homo_sapiens_assembly19.fasta.fai
C:\Bioinformatics_Project\reference_genomes\Homo_sapiens_assembly19.dict

Notes:
You only need to download these once. They can be reused across multiple runs.
Make sure file paths match the locations specified in your workflow scripts.
