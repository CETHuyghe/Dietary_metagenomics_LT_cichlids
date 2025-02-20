###########################################################
# 01. Prepare dataset
###########################################################
# RScript to clean and prepare metagenomic data for analysis
# For reads classified using Kraken2 and the RefSeq v205 database, from the digestive system

# 1.1 Prepare metadata
################################################################################

# Set working directory
setwd("/scicore/home/salzburg/huyghe0000/sinergia/DNA_gut/Kraken_combined/Kraken_NoHostSP_R/")
# Load field data and metadata of the cichlid samples and species
all_metadata <- read.delim(file="../Metadata_DNA_15112023.txt",sep="\t",check.names=FALSE, row.names = 1)

# 1.1.1 Filter samples from metadata
library(plyr)
library(dplyr)

# Remove reference and hybrid samples (diet reference/hybrid/miscluster tree)
all_metadata_1 <- all_metadata[ ! all_metadata$TissueID == 'Refgut',]
all_metadata_1 <- all_metadata_1[ ! all_metadata_1$SpecimenID == 'TYD7',]
all_metadata_1 <- all_metadata_1[ ! all_metadata_1$TissueTubeID %in% c('TXD2','VIE9','TZB5'),]

# 1.1.2 Calculate the relative gut length and Zihler's index

# Calculate species means of standard body length and weigth
SL_vars <- ddply(all_metadata_1, .(SpeciesID), summarize, SL_mean=mean(SL), SL_std_dev=sd(SL))
Weigth_vars <- ddply(all_metadata_1, .(SpeciesID), summarize,Weight_mean=mean(Weight), Weight_std_dev=sd(Weight))
# Load genome-wide phylogenetic tree from Ronco et al. (2021)
ActualPhylogeny = "/scicore/home/salzburg/huyghe0000/sinergia/DNA_gut/Kraken_combined/b1_with_Oretan.tre"
pruned.tree = ape::read.tree(ActualPhylogeny)
# Define the cichlid species part of the analysis
spp <- all_metadata_1[order(all_metadata_1$Tribe),]
species <- unique(spp$SpeciesID)
# Keep only the species of this study in the tree
pruned.tree = drop.tip(pruned.tree, pruned.tree$tip.label[-which(pruned.tree$tip.label %in% species)])
# Calculate Zihler's Index (Zihler) from the field data in this study, using the gut length and weight
LNG_CETH <- all_metadata_1[c("SpeciesID","SL","TL","Weight","GL","Tribe")]
LNG_CETH$Zihler <- (LNG_CETH$GL)/((sqrt(sqrt(LNG_CETH$Weight))))
# Calculate the relative gut length (RLG), with total gut length over standard body length
LNG_CETH$RGL <- LNG_CETH$GL/LNG_CETH$SL
# Calculate species means and standard deviation of Zihler's Index and relative gut length
ZHL_vars_CETH <- ddply(LNG_CETH, .(SpeciesID), summarize, ZHL_mean=mean(Zihler), ZHL_std_dev=sd(Zihler))
RGL_vars_CETH <- ddply(LNG_CETH, .(SpeciesID), summarize, RGL_mean=mean(RGL), RGL_std_dev=sd(RGL))

# 1.1.4 Add morphometric data from Ronco et al. (2021)

# Load Principal Component PC1 and PC2 scores of Body (BDY), Lower Pharyngeal Jaw (LPJ) and Upper Oral Jaw (UOJ) shape from Ronco et al. (2021)
LPj <- read.delim(file="../PCA_pc_scores_LPJ.txt",sep="\t",check.names=FALSE, row.names = 1)
UOJ <- read.delim(file="../PCA_pc_scores_UOJ.txt",sep="\t",check.names=FALSE, row.names = 1)
BDY <- read.delim(file="PCA_pc12_scores_body copy.txt",sep="\t",check.names=FALSE, row.names = 1)
# Calculate their species means and standard deviation
LPJ_vars <- ddply(LPj, .(sp), summarize, LPJ_PC1_mean=mean(LPJ_PC1), LPJ_PC1_std_dev=sd(LPJ_PC1),LPJ_PC2_mean=mean(LPJ_PC2), LPJ_PC2_std_dev=sd(LPJ_PC2))
UOJ_vars <- ddply(UOJ, .(sp), summarize, UOJ_PC1_mean=mean(UOJ_PC1), UOJ_PC1_std_dev=sd(UOJ_PC1),UOJ_PC2_mean=mean(UOJ_PC2), UOJ_PC2_std_dev=sd(UOJ_PC2))
BDY_vars <- ddply(BDY, .(sp), summarize, BDY_PC1_mean=mean(PC1), BDY_PC1_std_dev=sd(PC1),BDY_PC2_mean=mean(PC2), BDY_PC2_std_dev=sd(PC2))

# 1.2 Prepare abundance matrix
################################################################################

# Read table generated using combine_kreport.py
matrix0 <- read.delim(file="../kraken_combined_RefSeq_NoHost_30_report.txt",sep="\t") 
# Rename the rownames by NCBI TaxID
rownames(matrix0) <- matrix0$taxid
# Remove unnecessary columns
# Includes % total reads, combined number of reads (including reads within subtree), combined number of reads (only at this level), taxonomic classification level, name of level
matrix <- as.data.frame(subset(matrix0, select = -c(lvl_type,name,X.perc,tot_all,tot_lvl,taxid)))
rm(matrix0)
# Remove TaxID's with no reads
matrix <- matrix[rowSums(matrix[])>0,]
# Remove of each sample the total nr of reads at tax level, only use the nr of reads at that specific level
matrix <- as.data.frame(t(matrix))
matrix <- dplyr::filter(matrix, grepl('_lvl', rownames(matrix))) 
rownames(matrix) <- sub("_lvl", "", rownames(matrix))
matrix <- as.data.frame(t(matrix))
# Remove taxIDs without reads at that specific level
matrix <- matrix[rowSums(matrix[])>0,]
# Replace colnames by sample names
# Use file with smaple names in same order
colnames <- read.delim(file="../samples_NoHost_30_report.txt",sep="\t" )
colnames(matrix) <- colnames$SampleID

# 1.3 Prepare taxonomy table
################################################################################

# 1.3.1 Construct basic taxonomy table
library(taxonomizr)

# Extract TaxID's from the abundance matrix
taxaId <- rownames(matrix)
# Make table with the whole NCBI taxonomy of the specific TaxID's
TaxID <- getTaxonomy(taxaId,'../taxonomy_ncbi/accessionTaxa.sql')
TaxID <- data.frame(TaxID)
# Only keep Eukaryota for food source identification
Euk <- rownames(TaxID[TaxID$superkingdom == 'Eukaryota',])
Tax_table <- TaxID[rownames(TaxID) %in% Euk, ]

# 1.3.2 Filter potential host taxa

# Remove taxa which might contain host sequences, including cichliformes
Tax_table_strict <- Tax_table[ ! Tax_table$order %in% "Cichliformes", ]
# Identify taxa which might contain host sequences, including unidentified fish/chordata taxa
Tax_table_inter <- Tax_table_strict[ Tax_table_strict$phylum %in% "Chordata", ]
Tax_table_inter_2 <- Tax_table_inter[is.na(Tax_table_inter$order), ] 
Tax_table_inter_3 <- Tax_table_inter_2[is.na(Tax_table_inter_2$family), ] 
OnlChr <- rownames(Tax_table_inter_3[is.na(Tax_table_inter_3$class), ] )
OnlAct <- rownames(Tax_table_inter_3[Tax_table_inter_3$class == "Actinopteri", ] )
# Remove these TaxID's from the dataset
Tax_table_strict <- Tax_table_strict[ ! rownames(Tax_table_strict) %in% OnlChr, ]
Tax_table_strict <- Tax_table_strict[ ! rownames(Tax_table_strict) %in% OnlAct, ]
# Remove taxa not identified lower than Eukaryota
Tax_table_strict0 <- Tax_table_strict[is.na(Tax_table_strict$phylum), ] 
Tax_table_strict00 <- Tax_table_strict0[is.na(Tax_table_strict0$class), ] 
Tax_table_strict000 <- Tax_table_strict00[is.na(Tax_table_strict00$order), ] 
Tax_table_strict <- Tax_table_strict[ ! rownames(Tax_table_strict) %in% rownames(Tax_table_strict000), ]

# 1.3.3 Compare datasets

# Remove spaces in rownames of the tax and abundance table that might cause errors
row.names(matrix_6) <- gsub(' ', '', row.names(matrix_6))
row.names(Tax_table_strict) <- gsub(' ', '', row.names(Tax_table_strict))
# Keep only these rows in abundance table for which we have a taxonomy
matrix_abund_0 <- matrix_6[ rownames(matrix_6) %in% rownames(Tax_table_strict), ]
# Keep only samples from the metadata in the abundance table
matrix_abund_0 <- matrix_abund_0[ ,colnames(matrix_abund_0) %in% rownames(all_metadata_1) ]

# 1.3.4 Filter count matrix by number of reads

# Filter samples and OTU's
matrix_abund_1 <- matrix_abund_0
matrix_abund_1 <- matrix_abund_1[rowSums(matrix_abund_1) > 10,]
matrix_abund_1 <- matrix_abund_1[,colSums(matrix_abund_1) >= 100]
matrix_abund_2 <- matrix_abund_1[,colSums(matrix_abund_1 > 20) >= 2]
# Keep only remaining OTUs and samples in taxonomy
Tax_table_strict_2 <- Tax_table_strict[ rownames(Tax_table_strict) %in% rownames(matrix_abund_2), ]

# 1.3.5 Correct taxonomy table

# Add order to missing data of Actinopterygii
Tax_table_strict_2[rownames(subset(Tax_table_strict_2, family=='Pomacentridae')),"order"] <- "Ovalentaria incertae sedis Pomacentridae"
Tax_table_strict_2[rownames(subset(Tax_table_strict_2, family=='Ambassidae')),"order"] <- "Ovalentaria incertae sedis Ambassidae"
Tax_table_strict_2[rownames(subset(Tax_table_strict_2, family=='Opistognathidae')),"order"] <- "Ovalentaria incertae sedis Opistognathidae"
Tax_table_strict_2[rownames(subset(Tax_table_strict_2, genus=='Lates')),"family"] <- "Latidae"
Tax_table_strict_2[rownames(subset(Tax_table_strict_2, family=='Latidae')),"order"] <- "Carangaria incertae sedis Latidae"
Tax_table_strict_2[rownames(subset(Tax_table_strict_2, family=='Toxotidae')),"order"] <- "Carangaria incertae sedis Toxotidae"
Tax_table_strict_2[rownames(subset(Tax_table_strict_2, family=='Menidae')),"order"] <- "Carangaria incertae sedis Menidae"
Tax_table_strict_2[rownames(subset(Tax_table_strict_2, family=='Lactariidae')),"order"] <- "Carangaria incertae sedis Lactariidae"
Tax_table_strict_2[rownames(subset(Tax_table_strict_2, family=='Sciaenidae')),"order"] <- "Eupercaria incertae sedis Sciaenidae"
Tax_table_strict_2[rownames(subset(Tax_table_strict_2, family=='Moronidae')),"order"] <- "Eupercaria incertae sedis Moronidae"
Tax_table_strict_2[rownames(subset(Tax_table_strict_2, family=='Emmelichthyidae')),"order"] <- "Eupercaria incertae sedis Emmelichthyidae"
Tax_table_strict_2[rownames(subset(Tax_table_strict_2, family=='Pomacanthidae')),"order"] <- "Eupercaria incertae sedis Pomacanthidae"
Tax_table_strict_2[rownames(subset(Tax_table_strict_2, family=='Scatophagidae')),"order"] <- "Eupercaria incertae sedis Scatophagidae"
Tax_table_strict_2[rownames(subset(Tax_table_strict_2, family=='Sillaginidae')),"order"] <- "Eupercaria incertae sedis Sillaginidae"
# Add phylum to Dinophyceae
Tax_table_strict_2[rownames(subset(Tax_table_strict_2, class == 'Dinophyceae')),"phylum"] <- "Myzozoa"
# Reassign Octopoda to Mollusca
Tax_table_strict_2[rownames(subset(Tax_table_strict_2, class=='Cephalopoda')),"order"] <- NA
Tax_table_strict_2[rownames(subset(Tax_table_strict_2, class=='Cephalopoda')),"family"] <- NA
Tax_table_strict_2[rownames(subset(Tax_table_strict_2, class=='Cephalopoda')),"genus"] <- NA
Tax_table_strict_2[rownames(subset(Tax_table_strict_2, class=='Cephalopoda')),"class"] <- NA
# Reassign Anthoathacata to Limnomedusae
Tax_table_strict_2[rownames(subset(Tax_table_strict_2, order == 'Anthoathecata')),"family"] <- "Olindiidae"
Tax_table_strict_2[rownames(subset(Tax_table_strict_2, order == 'Anthoathecata')),"genus"] <- NA
Tax_table_strict_2[rownames(subset(Tax_table_strict_2, order == 'Anthoathecata')),"species"] <- NA
Tax_table_strict_2[rownames(subset(Tax_table_strict_2, order == 'Anthoathecata')),"order"] <- "Limnomedusae"

# 1.3.6 Remove taxa taxonomy table

# Remove Amphibia, Aves, Mammalia, Reptilia and Chondrichthyes
Tax_table_strict_2 <- Tax_table_strict_2[ ! Tax_table_strict_2$class %in% c("Aves","Amphibia","Mammalia","Chondrichthyes","Lepidosauria"), ]
Tax_table_strict_2 <- Tax_table_strict_2[ ! Tax_table_strict_2$order %in% c("Crocodylia","Testudines","Petromyzontiformes","Amphioxiformes","Stolidobranchia","Phlebobranchia"), ]
# Remove Actinopterygii orders not in LT, potentially due to host
Tax_table_strict_2 <- Tax_table_strict_2[ ! Tax_table_strict_2$order %in% c("Atheriniformes","Beloniformes","Blenniiformes","Ovalentaria incertae sedis Pomacentridae","valentaria incertae sedis Ambassidae","Carangaria incertae sedis Toxotidae","Carangaria incertae sedis Menidae","Carangiformes","Pleuronectiformes","Istiophoriformes","Eupercaria incertae sedis Sciaenidae","Eupercaria incertae sedis Moronidae","Spariformes","Labriformes","Centrarchiformes","Acropomatiformes","Lutjaniformes","Priacanthiformes","Gobiiformes","Kurtiformes","Batrachoidiformes","Syngnathiformes","Scombriformes","Holocentriformes","Gadiformes","Myctophiformes","Salmoniformes","Esociformes","Gymnotiformes","Elopiformes","Anguilliformes","Semionotiformes","Acipenseriformes","Ovalentaria incertae sedis Ambassidae","Perciformes","Coelacanthiformes"), ]
# Remove eDNA as ref
Tax_table_strict_2 <- Tax_table_strict_2[ ! Tax_table_strict_2$species %in% "Mollusca environmental sample", ]
# Remove Animalia taxa considered parasitic, analyzed separately
Tax_table_strict_2 <- Tax_table_strict_2[ ! Tax_table_strict_2$phylum %in% c("Acanthocephala","Platyhelminthes","Nematoda"), ]
# Remove Protista, analyzed separately
Tax_table_strict_2 <- Tax_table_strict_2[ ! Tax_table_strict_2$phylum %in% c("Ciliophora","Myzozoa","Fornicata","Apicomplexa"), ]
# Remove Fungi, analyzed separately
Tax_table_strict_2 <- Tax_table_strict_2[ ! Tax_table_strict_2$phylum %in% c("Microsporidia","Mucoromycota","Oomycota","Ascomycota"), ]
Tax_table_strict_2 <- Tax_table_strict_2[ ! Tax_table_strict_2$order %in% "Malasseziales", ]

# 1.4 Make final Phyloseq object
################################################################################

# 1.4.1 Keep same samples and taxa over all datasets

# Keep only these rows in abundance table for which we have a taxonomy
matrix_abund_2 <- matrix_abund_2[ rownames(matrix_abund_2) %in% rownames(Tax_table_strict_2), ]
# Remove samples with zero counts due to the removal of taxa
matrix_abund_2 <- matrix_abund_2[,colSums(matrix_abund_2) >= 1]
# Keep only samples from the abundance table in the metadata
all_metadata_2 <- all_metadata_1[ rownames(all_metadata_1) %in% colnames(matrix_abund_2) ,]

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
# Make phyloseq object
physeq <- phyloseq::phyloseq(OTU, TAX, samples)
physeq





