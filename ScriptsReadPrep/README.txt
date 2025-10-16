# ScriptsReadPrep

This folder contains the bash scripts used to prepare the fasta files needed for the downstream identification of the diet items using Kraken2 and the RefSeq database. Raw sequences were preprocessed and quality checked, and reads identified as host-derived, human-derived or microbiome-derived were removed using the 4 script contained in this folder.

# Folder content

## Bash scripts

1. 01_Script_fastp_Remove_Orenil

This script concatenates all raw fasta files per sample from the Illumina run, and uses fastp to preprocess and quality control these raw reads. It further identifies and removes reads mapping to the Oreochromis niloticus reference genome using bwa and SamTools.  

2. 02_Script_Kraken2_Remove_Human

This script identifies and removes reads identified as human, using Kraken2.

3. 03_Script_Kraken2_Remove_MB

This script identifies and removes reads identified as part of the microbiome, using Kraken2.

4. 04_Script_Remove_HostSP

This script identifies and removes reads identified as from the specific host cichlid species, using bwa and SamTools.

# Downstream pipeline

The resulting fasta files are further used to identify diet items in the reads, following the scripts in the folder ScriptsKrakenClassify
