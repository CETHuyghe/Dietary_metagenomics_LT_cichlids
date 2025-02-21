###########################################################
# 01. Prepare Focus dataset
###########################################################
# RScript to clean and prepare metagenomic data for analysis
# For reads classified using Kraken2 and the custom focus database, from the digestive system

# Set working directory
setwd("")

# 1.1 Prepare metadata
################################################################################

# 1.1.1 Load metadata
all_metadata <- read.delim(file="CichlidLT_Diet_MetaData_CETH.txt",sep="\t",check.names=FALSE, row.names = 1)

# 1.1.2 Load phylogeny

# Load phylogeny
library(ape)
ActualPhylogeny = "/scicore/home/salzburg/huyghe0000/sinergia/DNA_gut/Kraken_combined/b1_with_Oretan.tre"
pruned.tree = ape::read.tree(ActualPhylogeny)

# 1.2 Prepare abundance matrix
################################################################################
library(plyr)
library(dplyr)

# Read table from combine_kreport.py
matrix <- read.delim(file="CichlidLT_Diet_Focus_DataMatrix_CETH.txt",sep="\t")
# Make rownames NCBI TaxID
rownames(matrix) <- matrix$taxid
# Remove % total reads, combined number of reads (including reads within subtree), combined number of reads (only at this level), taxonomic classification level, name of level
matrix_2 <- as.data.frame(subset(matrix, select = -c(lvl_type,name,X.perc,tot_all,tot_lvl,taxid)))
rm(matrix)
# Remove majority taxid's
matrix_2 <- matrix_2[rowSums(matrix_2[])>0,]
# Remove of each sample the total nr of reads at tax level, only use the nr of reads at that specific level, otherwise more reads in total
matrix_2 <- as.data.frame(t(matrix_2))
matrix_2 <- dplyr::filter(matrix_2, grepl('_lvl', rownames(matrix_2))) 
rownames(matrix_2) <- sub("_lvl", "", rownames(matrix_2))
matrix_2 <- as.data.frame(t(matrix_2))
# Remove taxIDs without reads at that specific level
matrix_2 <- matrix_2[rowSums(matrix_2[])>0,]
# Replace colnames 
colnames <- read.delim(file="CichlidLT_Diet_Focus_DataMatrix_SampleNames_CETH.txt",sep="\t" )
matrix_6 <- matrix_2
rm(matrix_2)
colnames(matrix_6) <- colnames$SampleID

# 1.3 Prepare taxonomy table
################################################################################

# 1.3.1 Construct basic taxonomy table
library(taxonomizr)

# First extract TaxID's
taxaId <- rownames(matrix_6)
# Make table with the whole NCBI taxonomy of the specific TaxID's
#prepareDatabase('accessionTaxa.sql')
TaxID <- getTaxonomy(taxaId,'../taxonomy_ncbi/accessionTaxa.sql')
TaxID <- data.frame(TaxID)
# Only keep Eukaryota for Diet
Euk <- rownames(TaxID[TaxID$superkingdom == 'Eukaryota',])
Tax_table <- TaxID[rownames(TaxID) %in% Euk, ]

# 1.3.2 Filter potential host taxa

# Remove Cichlids
Tax_table_strict <- Tax_table[ ! Tax_table$order %in% "Cichliformes", ]
# Also remove Chordata and Actinopterygi that are not identified to order level, might incl cichlids
Tax_table_inter <- Tax_table_strict[ Tax_table_strict$phylum %in% "Chordata", ]
Tax_table_inter_2 <- Tax_table_inter[is.na(Tax_table_inter$order), ] 
Tax_table_inter_3 <- Tax_table_inter_2[is.na(Tax_table_inter_2$family), ] 
OnlChr <- rownames(Tax_table_inter_3[is.na(Tax_table_inter_3$class), ] )
OnlAct <- rownames(Tax_table_inter_3[Tax_table_inter_3$class == "Actinopteri", ] )
Tax_table_strict <- Tax_table_strict[ ! rownames(Tax_table_strict) %in% OnlChr, ]
Tax_table_strict <- Tax_table_strict[ ! rownames(Tax_table_strict) %in% OnlAct, ]
# Remove taxa not identified lower than Eukaryota
Tax_table_strict0 <- Tax_table_strict[is.na(Tax_table_strict$phylum), ] 
Tax_table_strict00 <- Tax_table_strict0[is.na(Tax_table_strict0$class), ] 
Tax_table_strict000 <- Tax_table_strict00[is.na(Tax_table_strict00$order), ] 

Tax_table_strict <- Tax_table_strict[ ! rownames(Tax_table_strict) %in% rownames(Tax_table_strict000), ]

# 1.3.3 Correct taxonomy table

Tax_table_strict[rownames(subset(Tax_table_strict, order=='Ceratodontiformes')),"class"] <- "Dipneusti"
Tax_table_strict[rownames(subset(Tax_table_strict, order=='Crocodylia')),"class"] <- "Reptilia"
Tax_table_strict[rownames(subset(Tax_table_strict, order=='Opistognathidae')),"class"] <- "Reptilia"
Tax_table_strict[rownames(subset(Tax_table_strict, family=='Lymnaeidae')),"order"] <- "Hygrophila"
Tax_table_strict[rownames(subset(Tax_table_strict, class=='Eustigmatophyceae')),"phylum"] <- "Ochrophyta"
Tax_table_strict[rownames(subset(Tax_table_strict, class=='Dinophyceae')),"order"] <- "Myzozoa"
Tax_table_strict[rownames(subset(Tax_table_strict, genus=='Lates')),"family"] <- "Latidae"
Tax_table_strict[rownames(subset(Tax_table_strict, family=='Latidae')),"order"] <- "Carangaria incertae sedis Latidae"

# 1.3.4 Compare datasets

# Remove spaces in rownames that might cause errors
row.names(matrix_6) <- gsub(' ', '', row.names(matrix_6))
row.names(Tax_table_strict) <- gsub(' ', '', row.names(Tax_table_strict))
# Keep only these rows in abundance table for which we have a taxonomy
matrix_abund_0 <- matrix_6[ rownames(matrix_6) %in% rownames(Tax_table_strict), ]
# Remove the reference samples
refgut <- rownames(all_metadata[all_metadata$TissueID == 'Refgut',])
matrix_abund_0 <- matrix_abund_0[ , ! colnames(matrix_abund_0) %in% refgut]
# Remove hybrid sample
Hybrid <- rownames(all_metadata[all_metadata$SpecimenID == 'TYD7',])
matrix_abund_0 <- matrix_abund_0[ , ! colnames(matrix_abund_0) %in% Hybrid ]
# Remove problematic samples from mitogen tree
rm1 <- rownames(all_metadata[all_metadata$TissueTubeID == 'TXD2',])
rm2 <- rownames(all_metadata[all_metadata$TissueTubeID == 'VIE9',])
rm3 <- rownames(all_metadata[all_metadata$TissueTubeID == 'TZB5',])
matrix_abund_0 <- matrix_abund_0[ , ! colnames(matrix_abund_0) %in% rm1 ]
matrix_abund_0 <- matrix_abund_0[ , ! colnames(matrix_abund_0) %in% rm2 ]
matrix_abund_0 <- matrix_abund_0[ , ! colnames(matrix_abund_0) %in% rm3 ]

# 1.3.5 Filter count matrix by number of reads

# Check total nr of reads per sample
tot_nr_0 <- as.data.frame( colSums(matrix_abund_0))
rownames(tot_nr_0) <- colnames(matrix_abund_0)
tot_nr_0$Sample <- rownames(tot_nr_0)
tot_nr_0$read_count <- tot_nr_0$`colSums(matrix_abund_0)`
# Check total nr of reads per OTU
tot_nr_OTU <- as.data.frame(rowSums(matrix_abund_0))
tot_nr_OTU$Sample <- rownames(tot_nr_OTU)
tot_nr_OTU$read_count <- tot_nr_OTU$`rowSums(matrix_abund_0)`
# Filter samples and taxa
matrix_abund_1 <- matrix_abund_0
matrix_abund_1 <- matrix_abund_1[rowSums(matrix_abund_1) > 10,]
matrix_abund_1 <- matrix_abund_1[,colSums(matrix_abund_1) >= 100]
matrix_abund_2 <- matrix_abund_1[,colSums(matrix_abund_1 > 20) >= 2]

# 1.3.6 Remove taxa taxonomy table

# To taxonomy after filtering
Tax_table_strict_2 <- Tax_table_strict[ rownames(Tax_table_strict) %in% rownames(matrix_abund_2), ]
# Remove our jellyfish and Sponge, too contaminated
Tax_table_strict_2 <- Tax_table_strict_2[ ! rownames(Tax_table_strict_2) %in% c("6050","1914956"), ]
# Remove parasites
# Annelid identified as leech kept as annelid, no intestinal leech parasites
Tax_table_strict_2 <- Tax_table_strict_2[ ! Tax_table_strict_2$phylum %in% c("Acanthocephala","Platyhelminthes","Nematoda"), ]
# Remove protists
Tax_table_strict_2 <- Tax_table_strict_2[ ! Tax_table_strict_2$phylum %in% c("Ciliophora","Myzozoa","Fornicata","Apicomplexa"), ]

# 1.4 Make final Phyloseq object
################################################################################

# 1.4.1 Keep same samples and taxa over all datasets

# Keep only these rows in abundance table for which we have a taxonomy
matrix_abund_2 <- matrix_abund_2[ rownames(matrix_abund_2) %in% rownames(Tax_table_strict_2), ]
# Remove samples with zero counts
matrix_abund_2 <- matrix_abund_2[,colSums(matrix_abund_2) >= 1]
# Create metadata
all_metadata_2 <- all_metadata[ rownames(all_metadata) %in% colnames(matrix_abund_2) ,]

# 1.4.2 Make Phyloseq object
library(phyloseq)
library(vegan)

# Specify datasets needed
otu_mat <- matrix_abund_2
tax_mat <- Tax_table_strict_2
samples_df <- all_metadata_2
# Transform otu table to matrix with numeric values
otu_mat <- as.matrix(otu_mat)
class(otu_mat) <- "numeric"
# Transform taxonomy table to matrix
tax_mat <- as.matrix(tax_mat)
# Transform matrices to phyloseq objects
OTU = phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = phyloseq::tax_table(tax_mat)
samples = phyloseq::sample_data(samples_df)
physeq <- phyloseq::phyloseq(OTU, TAX, samples)
physeq


