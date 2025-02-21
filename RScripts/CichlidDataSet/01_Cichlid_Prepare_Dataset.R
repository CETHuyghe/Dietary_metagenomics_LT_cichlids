###########################################################
# 01. Prepare Cichlid dataset
###########################################################
# RScript to clean and prepare metagenomic data for analysis
# For reads classified using Kraken2 and the custom cichlid database, from the digestive system

# Set working directory
setwd("")

# 1.1 Prepare metadata
################################################################################

# 1.1.1 Load metadata

# Load metadata
all_metadata <- read.delim(file="CichlidLT_Diet_MetaData_CETH",sep="\t",check.names=FALSE, row.names = 1)

# 1.2 Prepare abundance matrix
################################################################################

# 1.2.1 Setup matrix

# Read cichlid genome table from combine_kreport.py
matrix <- read.delim(file="CichlidLT_Diet_Cichlidae_DataMatrix_CETH.txt",sep="\t") # Includes 
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
# Replace colnames by sample names same order
colnames <- read.delim(file="CichlidLT_Diet_Cichlidae_DataMatrix_SampleNames_CETH.txt",sep="\t" )
matrix_2 <- matrix_2
colnames(matrix_2) <- colnames$SampleID

# 1.3 Prepare taxonomy table
################################################################################

# 1.3.1 Construct basic taxonomy table
library(taxonomizr)

# First extract TaxID's
taxaId <- rownames(matrix_6)
# Make table with the whole NCBI taxonomy of the specific TaxID's
TaxID <- getTaxonomy(taxaId,'../taxonomy_ncbi/accessionTaxa.sql')
Tax_table <- data.frame(TaxID)

# 1.3.2 Filter potential host taxa

# Remove Cichlids not identified lower than family
Tax_table_strict0 <- Tax_table[is.na(Tax_table$genus), ] 
Tax_table_strict <- Tax_table[ ! rownames(Tax_table) %in% rownames(Tax_table_strict0), ]

# 1.3.3 Compare datasets

# Remove spaces in rownames that might cause errors
row.names(matrix_6) <- gsub(' ', '', row.names(matrix_6))
row.names(Tax_table_strict) <- gsub(' ', '', row.names(Tax_table_strict))
# Keep only these rows in abundance table for which we have a taxonomy
matrix_abund_0 <- matrix_6[ rownames(matrix_6) %in% rownames(Tax_table_strict), ]

# 1.3.4 Remove samples

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

# 1.3.6 Filter out taxa

# Remove 4 outlier cichlid taxa: "278786","158781", "2591303", "158760"
Tax_table_strict[  rownames(Tax_table_strict) %in% c("278786","158781", "2591303", "158760"), ]
Tax_table_strict_1 <- Tax_table_strict[ ! rownames(Tax_table_strict) %in% c("278786","158781", "2591303", "158760"), ]
# Remove non-LT cichlids
Tax_table_strict_1 <- Tax_table_strict_1[ ! Tax_table_strict_1$species %in% c("Benitochromis conjunctus","Heterotilapia buettikoferi","Coptodon bakossiorum" , "Apistogramma diplotaenia", "Steatocranus sp. 'ultraslender' AB-2019", "Pelmatolapia mariae","Tilapia sparmanii", "Tilapia brevimanus","Pharyngochromis acuticeps","Astatoreochromis straeleni","Orthochromis mazimeroensis","Orthochromis indermauri","Thoracochromis brauschi","Stomatepia pindu","Pungu maclareni", "Sarotherodon caroli","Sarotherodon steinbachi","Sarotherodon lohbergeri","Oreochromis malagarasi","Orthochromis malagaraziensis","Orthochromis uvinzae","Astatotilapia flaviijosephi","Astatotilapia stappersii","Serranochromis macrocephalus","Sargochromis carlottae","Pelvicachromis taeniatus","Hemichromis elongatus","Etia nguti","Heterochromis multidens","Australoheros scitulus","Amphilophus zaliosus","Etroplus canarensis","Paratilapia polleni","Ptychochromis oligocanthus"), ]

# 1.3.7 Add tribes

# Add tribes to cichlid taxa
# Load tribe info
all_tribes <- read.delim(file="All_tribes_Cichlids.txt",sep="\t",check.names=FALSE, row.names = 1)
# Add tribe using tribe info
Tax_table_strict_3 <- Tax_table_strict_1
Tax_table_strict_3$tribe <-with(all_tribes, Tribe[match(Tax_table_strict_3$species,rownames(all_tribes))])
# Add manually for more complicated cases, or not as described species
Tax_table_strict_3[rownames(subset(Tax_table_strict_3, genus=='Neolamprologus')),"tribe"] <- "Lamprologini"
Tax_table_strict_3[rownames(subset(Tax_table_strict_3, genus=='Telmatochromis')),"tribe"] <- "Lamprologini"
Tax_table_strict_3[rownames(subset(Tax_table_strict_3, genus=='Chalinochromis')),"tribe"] <- "Lamprologini"
Tax_table_strict_3[rownames(subset(Tax_table_strict_3, genus=='Lepidiolamprologus')),"tribe"] <- "Lamprologini"
Tax_table_strict_3[rownames(subset(Tax_table_strict_3, genus=='Lamprologus')),"tribe"] <- "Lamprologini"
Tax_table_strict_3[rownames(subset(Tax_table_strict_3, genus=='Altolamprologus')),"tribe"] <- "Lamprologini"
Tax_table_strict_3[rownames(subset(Tax_table_strict_3, genus=='Julidochromis')),"tribe"] <- "Lamprologini"
Tax_table_strict_3[rownames(subset(Tax_table_strict_3, genus=='Bathybates')),"tribe"] <- "Bathybatini"
Tax_table_strict_3[rownames(subset(Tax_table_strict_3, genus=='Hemibates')),"tribe"] <- "Bathybatini"
Tax_table_strict_3[rownames(subset(Tax_table_strict_3, genus=='Benthochromis')),"tribe"] <- "Benthochromini"
Tax_table_strict_3[rownames(subset(Tax_table_strict_3, genus=='Xenotilapia')),"tribe"] <- "Ectodini"
Tax_table_strict_3[rownames(subset(Tax_table_strict_3, genus=='Ophthalmotilapia')),"tribe"] <- "Ectodini"
Tax_table_strict_3[rownames(subset(Tax_table_strict_3, genus=='Callochromis')),"tribe"] <- "Ectodini"
Tax_table_strict_3[rownames(subset(Tax_table_strict_3, genus=='Cyathopharynx')),"tribe"] <- "Ectodini"
Tax_table_strict_3[rownames(subset(Tax_table_strict_3, genus=='Ectodus')),"tribe"] <- "Ectodini"
Tax_table_strict_3[rownames(subset(Tax_table_strict_3, genus=='Petrochromis')),"tribe"] <- "Tropheini"
Tax_table_strict_3[rownames(subset(Tax_table_strict_3, genus=='Tropheus')),"tribe"] <- "Tropheini"
Tax_table_strict_3[rownames(subset(Tax_table_strict_3, genus=='Pseudosimochromis')),"tribe"] <- "Tropheini"
Tax_table_strict_3[rownames(subset(Tax_table_strict_3, genus=='Cyprichromis')),"tribe"] <- "Cyprichromini"
Tax_table_strict_3[rownames(subset(Tax_table_strict_3, genus=='Paracyprichromis')),"tribe"] <- "Cyprichromini"
Tax_table_strict_3[rownames(subset(Tax_table_strict_3, genus=='Cyphotilapia')),"tribe"] <- "Cyphotilapiini"
Tax_table_strict_3[rownames(subset(Tax_table_strict_3, genus=='Limnochromis')),"tribe"] <- "Limnochromini"
Tax_table_strict_3[rownames(subset(Tax_table_strict_3, genus=='Oreochromis')),"tribe"] <- "Oreochromini"
Tax_table_strict_3[rownames(subset(Tax_table_strict_3, genus=='Haplochromis')),"tribe"] <- "Haplochromini"
Tax_table_strict_3[rownames(subset(Tax_table_strict_3, genus=='Sarotherodon')),"tribe"] <- "Oreochromis"
Tax_table_strict_3[rownames(subset(Tax_table_strict_3, genus=='Trematochromis')),"tribe"] <- "Haplochromini"
Tax_table_strict_3[rownames(subset(Tax_table_strict_3, genus=='Orthochromis')),"tribe"] <- "Haplochromini"
Tax_table_strict_3[rownames(subset(Tax_table_strict_3, genus=='Perissodus')),"tribe"] <- "Perissodini"
Tax_table_strict_3[rownames(subset(Tax_table_strict_3, genus=='Spathodus')),"tribe"] <- "Eretmodini"
Tax_table_strict_3[rownames(subset(Tax_table_strict_3, genus=='Greenwoodochromis')),"tribe"] <- "Limnochromini"
Tax_table_strict_3[rownames(subset(Tax_table_strict_3, genus=='Astatotilapia')),"tribe"] <- "Haplochromini"
Tax_table_strict_3[rownames(subset(Tax_table_strict_3, genus=='Plecodus')),"tribe"] <- "Perissodini"
Tax_table_strict_3[rownames(subset(Tax_table_strict_3, genus=='Haplotaxodon')),"tribe"] <- "Perissodini"
Tax_table_strict_3[rownames(subset(Tax_table_strict_3, genus=='Eretmodus')),"tribe"] <- "Eretmodini"
Tax_table_strict_3[rownames(subset(Tax_table_strict_3, genus=='Coptodon')),"tribe"] <- "Coptodini"
# Remove non-LT genera
Tax_table_strict_3 <- Tax_table_strict_3[ ! Tax_table_strict_3$genus %in% c("Tilapia", "Orthochromis", "Sarotherodon"), ]
# Relocate tribe column after family
library(plyr)
library(dplyr)
Tax_table_strict_2 <- Tax_table_strict_3 %>% relocate(tribe, .after = family)

# 1.4 Make final Phyloseq object
################################################################################

# 1.4.1 Keep same samples and taxa over all datasets

# Take to abundance
matrix_abund_1 <- matrix_abund_0[ rownames(matrix_abund_0) %in% rownames(Tax_table_strict_2), ]
# Filter samples and taxa
matrix_abund_1 <- matrix_abund_1[rowSums(matrix_abund_1) > 10,]
matrix_abund_2 <- matrix_abund_1[,colSums(matrix_abund_1) >= 20]
# To taxonomy after filtering
Tax_table_strict_2 <- Tax_table_strict_2[ rownames(Tax_table_strict_2) %in% rownames(matrix_abund_2), ]
# Filter metadata
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
# Make phyloseq object
physeq <- phyloseq::phyloseq(OTU, TAX, samples)
physeq





