By C.E.T. Huyghe et al.

# 1. Scripts concerning Metagenomic read processing and taxonomic assignment

## 1.1 ScriptsReadPrep
Clean sequences and remove (contamination from) host, human, viral, bacterial sources:

### 1.1.1 Script_fastp_Remove_Orenil
Clean sequence files and remove reads mapping to high quality Oreochromis niloticus genome

### 1.1.2 Script_Kraken2_Remove_Human
Remove reads mapping to human

### 1.1.3 Script_Kraken2_Remove_MB
Remove reads mapping to the standard Kraken2 microbiome database

### 1.1.4 Script_Remove_HostSP
Remove the reads mapping to the genome of the host species (not for reference samples)

## 1.2 ScriptsRefAssemble
Assemble the reference samples:

### 1.2.1 Script_ABySS_Hamburg
Script_ABySS_HamburgDe novo assemble

### 1.2.2 Script_FCS_Adaptor
Remove contaminants like adaptors 

### 1.2.3 Script_FCS_GX
Remove contaminants by other taxa 

## 1.3 ScriptsMitoGenCichlids
Generate Cichlid Genome Kraken2 Database:

### 1.3.1 Script_Download_NCBI_CichlidGenomes
Download genomes batch from NCBI

### 1.3.2 Script_Accession_to_TaxID_CichlidGenomes
Make file with Accession and TaxID of genomes

### 1.3.3 Script_Add_CichlidGenomes_KrakenDB
Add TaxID to genomes and add to Kraken2 database

### 1.3.4 Script_Build_CichlidGenomes_KrakenDB
Build Kraken2 Database

## 1.4 ScriptsFocusDB

Generate LT Focus Genome Kraken2 Database:
### 1.4.1. Script_Add_FocusGenomes_KrakenDB
Add TaxID to genomes and add to Kraken2 database
### 1.4.2. Script_Build_FocusGenomes_KrakenDB
Build Kraken2 Database

## 1.5 ScriptsCichlidDB
Generate mitochondrial genomes of the host cichlid from metagenomes

### 1.5.1 Script_Generate_Cichlid_Mitogenomes
Generate mitogenomes by mapping to reference, filter and consensus

### 1.5.2 Script_Metrics_Cichlid_Mitogenomes
Calculate coverage, mean depth and number of mapped reads

## 1.6 ScriptsKrakenClassify
Classify reads using Kraken2:

### 1.6.1 Script_Kraken2_Classify_RefSeqv205
Use RefSeq v205 reference database

### 1.6.2 Script_Kraken2_Classify_LTFocus
Use the custom Kraken2 database with a focus on LT taxa and sampled references

### 1.6.3 Script_Kraken2_Classify_Cichlid
Use the custom Kraken2 database with all available LT cichlid genomes from Ronco et al.,2021

### 1.6.4 Script_Kraken2_Combine_Reports
Combine reports of the Kraken2 output, do this for each of the above separately


## 1.7 ScriptsCalcNrReads
Other calculations

### 1.7.1 Script_Calc_Nr_Seq
Calculate number of sequences at each step along the workflow


# 2. RScripts concerning Metagenomic analysis of count matrices

## 2.1 RefSeqDataSet
Contains RScripts used for the analyses of the metagenomic dataset, based on the count matrix constructed by read classifications using Kraken2 and the RefSeq v205 reference database

## 2.2 FocusDataSet
Contains RScripts used for the analyses of the metagenomic dataset, based on the count matrix constructed by read classifications using Kraken2 and the custom LT Focus reference database

## 2.3 CichlidDataSet
Contains RScripts used for the analyses of the metagenomic dataset, based on the count matrix constructed by read classifications using Kraken2 and the custom cichlid reference database

## 2.4 Reused_Data
Datasets used for the analyses of the metagenomic dataset, originating from other studies (e.g. Ronco et al., 2021)
