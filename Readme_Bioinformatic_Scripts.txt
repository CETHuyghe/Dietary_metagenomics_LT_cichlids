When having the output from the DNA NovaSeq sequencing run.
By C.E.T. Huyghe (2024)

Clean sequences and remove (contamination from) host, human, viral, bacterial sources:
1. Clean sequence files and remove reads mapping to high quality Oreochromis niloticus genome
"Script_fastp_Remove_Orenil"
2. Remove reads mapping to human
"Script_Kraken2_Remove_Human"
3. Remove reads mapping to the standard Kraken2 microbiome database
"Script_Kraken2_Remove_MB"
4.Remove the reads mapping to the genome of the host species (not for reference samples)
"Script_Remove_HostSP"

Assemble the reference samples:
1. De novo assemble
"Script_ABySS_Hamburg"
2. Remove contaminants like adaptors 
"Script_FCS_Adaptor"
3. Remove contaminants by other taxa 
"Script_FCS_GX"

Generate Cichlid Genome Kraken2 Database:
1. Download genomes batch from NCBI
"Script_Download_NCBI_CichlidGenomes"
2. Make file with Accession and TaxID of genomes
"Script_Accession_to_TaxID_CichlidGenomes"
3. Add TaxID to genomes and add to Kraken2 database
"Script_Add_CichlidGenomes_KrakenDB"
4. Build Kraken2 Database
"Script_Build_CichlidGenomes_KrakenDB"

Generate LT Focus Genome Kraken2 Database:
1. Add TaxID to genomes and add to Kraken2 database
"Script_Add_FocusGenomes_KrakenDB"
2. Build Kraken2 Database
"Script_Build_FocusGenomes_KrakenDB"

Generate mitochondrial genomes of the host cichlid from metagenomes:
1. Generate mitogenomes by mapping to reference, filter and consensus
"Script_Generate_Cichlid_Mitogenomes"
2. Calculate coverage, mean depth and number of mapped reads
"Script_Metrics_Cichlid_Mitogenomes"

Classify reads using Kraken2:
1. Use RefSeq v205 reference database
"Script_Kraken2_Classify_RefSeqv205"
2. Use the custom Kraken2 database with a focus on LT taxa and sampled references
"Script_Kraken2_Classify_LTFocus"
3. Use the custom Kraken2 database with all available LT cichlid genomes from Ronco et al.,2021
"Script_Kraken2_Classify_Cichlid"
4. Combine reports of the Kraken2 output, do this for each of the above separately
"Script_Kraken2_Combine_Reports"

Other calculations:
1. Calculate number of sequences at each step along the workflow
"Script_Calc_Nr_Seq"