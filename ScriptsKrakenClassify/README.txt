# ScriptsKrakenClassify

This folder contains the bash scripts used to identify the reads in the fasta files, after quality control and non-target sequences removal. The leftover sequences resulting from the scripts in folder ScriptsReadPrep were identified using Kraken2 and the reference database specified in the script. Three different reference databases were used to identify the sequences in paralel, coming from 3 seperate bash scripts, and for each of them the resulting kraken reports were combined using a separate bash script to be used in downstream analyses.

# Folder content

## Bash scripts

1. 01_Script_Kraken2_Classify_RefSeqv205

This script identifies reads using Kraken2 and the RefSeq v205 database.  

2. 01_Script_Kraken2_Classify_LTFocus

This script identifies reads using Kraken2 and a custom-made reference database containing potential targeted diet species and self-sequenced diet items. See paper for more information.

3. 01_Script_Kraken2_Classify_Cichlid

This script identifies reads using Kraken2 and a custom-made database containing the genomes of the Lake Tanganyika cichlids generated in the CichlidX project. See paper for more information.

4. 02_Script_Kraken2_Combine_Reports

This script combines the separate Kraken report files per species into one count matrix which can be read by R.

# Downstream pipeline

The resulting count matrices are used for downstream analyses in R.
