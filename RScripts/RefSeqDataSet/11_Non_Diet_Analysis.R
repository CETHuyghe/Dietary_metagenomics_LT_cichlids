###########################################################
# 11. Non diet taxa analysis
###########################################################
# RScript to analyze the taxa not identified as diet

# 11.1 Protista
################################################################################

# 11.1.1 Select and filter data

# Filter abundance table before taxa filtering
matrix_abund_prot <- matrix_abund_1[rowSums(matrix_abund_1) >= 800,]
matrix_abund_prot <- matrix_abund_prot[rowSums(matrix_abund_prot > 20) >= 2,]
# Remove reference, hybrid and problematic samples
matrix_abund_prot <- matrix_abund_prot[ ,colnames(matrix_abund_prot) %in% rownames(all_metadata_1) ]
# Remove from taxonomy table
Tax_table_strict_prot <- Tax_table_strict[ rownames(Tax_table_strict) %in% rownames(matrix_abund_prot), ]
# Only take Protista taxa
Tax_table_strict_prot <- Tax_table_strict_prot[ Tax_table_strict_prot$phylum %in% c("Ciliophora","Myzozoa","Fornicata","Apicomplexa"), ]
# Keep only these rows in abundance table for which we have a taxonomy
matrix_abund_prot <- matrix_abund_prot[ rownames(matrix_abund_prot) %in% rownames(Tax_table_strict_prot), ]
# Remove samples with lower than 10 counts
matrix_abund_prot <- matrix_abund_prot[,colSums(matrix_abund_prot) >= 10]
# Keep only samples with Protista in metadata
all_metadata_prot <- all_metadata[ rownames(all_metadata) %in% colnames(matrix_abund_prot) ,]

# 11.1.2 Construct phyloseq object

# Specify datasets needed
otu_mat_prot <- matrix_abund_prot
tax_mat_prot <- Tax_table_strict_prot
samples_df_prot <- all_metadata_prot
# Transform otu table to matrix with numeric values
otu_mat_prot <- as.matrix(otu_mat_prot)
class(otu_mat_prot) <- "numeric"
# Transform taxonomy table to matrix
tax_mat_prot <- as.matrix(tax_mat_prot)
# Transform matrices to phyloseq objects
OTU_prot = phyloseq::otu_table(otu_mat_prot, taxa_are_rows = TRUE)
TAX_prot = phyloseq::tax_table(tax_mat_prot)
samples_prot = phyloseq::sample_data(samples_df_prot)
# Make phyloseq object
physeq_prot <- phyloseq::phyloseq(OTU_prot, TAX_prot, samples_prot)
physeq_prot

# 11.1.3 Define presence/absence at phylum level

# Subset to phylum level
physeq_prot_phylum <- tax_glom(physeq_prot, taxrank=rank_names(physeq_prot)[2], NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
# Need a min of 10 reads per OTU per sample
df_prot_phylum <- as.data.frame(otu_table(physeq_prot_phylum))
df_prot_phylum[df_prot_phylum < 10 ] <- 0
# Remove samples with lower than 10 counts
df_prot_phylum <- df_prot_phylum[,colSums(df_prot_phylum) >= 10]
# Only presence-abscence, to change values to 1
df_prot_phylum[df_prot_phylum > 1 ] <- 1
# Add again to Phyloseq
otu_table(physeq_prot_phylum) <- phyloseq::otu_table(df_prot_phylum, taxa_are_rows = TRUE)

# 11.1.4 Prepare and define bubble plot

# Add variable combining SpeciesID and TissueID for plot
variable1 = as.character(get_variable(physeq_prot_phylum, "SpeciesID"))
variable2 = as.character(get_variable(physeq_prot_phylum, "TissueID"))
sample_data(physeq_prot_phylum)$SpeciesID_TissueID <- mapply(paste0, variable1, variable2,collapse = "_")
# Merge per species, FG vs HG
physeq_prot_phylum_spp_tiss <- merge_samples(physeq_prot_phylum, "SpeciesID_TissueID", fun=sum) 
# Extract otu table
prot_phylum_spp_tiss = as(otu_table(physeq_prot_phylum_spp_tiss), "matrix")
# Rename TaxID's
colnames(prot_phylum_spp_tiss) <- with(Tax_table_strict_prot, phylum[match(colnames(prot_phylum_spp_tiss),rownames(Tax_table_strict_prot))])
# Melt matrix for graph
prot_phylum_spp_tiss_melt <- melt(t(prot_phylum_spp_tiss), id.vars = "SpeciesID_TissueID", variable.name = "Phylum")
prot_phylum_spp_tiss_melt$value <- as.numeric(prot_phylum_spp_tiss_melt$value)
# Put NA so only dots when there are values
prot_phylum_spp_tiss_melt[prot_phylum_spp_tiss_melt == 0] <- NA
# Add metadata
sample_data_phylum_prot <- as.data.frame(sample_data(physeq_prot_phylum))
prot_phylum_spp_tiss_melt$SpeciesID = with(sample_data_phylum_prot, SpeciesID[match(prot_phylum_spp_tiss_melt$Var2,sample_data_phylum_prot$SpeciesID_TissueID)])
prot_phylum_spp_tiss_melt$Tribe = with(sample_data_phylum_prot, Tribe[match(prot_phylum_spp_tiss_melt$Var2,sample_data_phylum_prot$SpeciesID_TissueID)])
prot_phylum_spp_tiss_melt$TissueID = with(sample_data_phylum_prot, TissueID[match(prot_phylum_spp_tiss_melt$Var2,sample_data_phylum_prot$SpeciesID_TissueID)])
# Rename FG HG
prot_phylum_spp_tiss_melt <- prot_phylum_spp_tiss_melt %>%
  mutate(TissueID = recode(TissueID, ContentFG = 'AG', ContentHG = 'PG' ))
# Combine protsite and tissue id for graph
prot_phylum_spp_tiss_melt$prot_tiss <- paste(prot_phylum_spp_tiss_melt$Var1,prot_phylum_spp_tiss_melt$TissueID)
# Order Cichlids in rigth way
prot_cichlid_phylum_host_ordering <- c("Cphgib","Cphfro","Ophven","PcybrN","Cypdwj","Varmoo","Telvit","Neobre")
# Rename column
prot_phylum_spp_tiss_melt$FOO <- prot_phylum_spp_tiss_melt$value
# Define colours by species' tribe
prot_cichlid_phylum_host_cols <- c("Cphgib"="#FDDF13","Cphfro"="#FDDF13","Ophven"="#9AB9D9","PcybrN"="#F04D29","Cypdwj"="#F04D29","Varmoo"="#C588BB","Telvit"="#C588BB","Neobre"="#C588BB")
# Make bubble plot
prot_plot <- ggplot(prot_phylum_spp_tiss_melt, aes(y=prot_tiss,x=factor(SpeciesID, level = prot_cichlid_phylum_host_ordering),colour=Tribe, size=FOO)) +
  geom_point() +
  scale_colour_manual(values=tribe_colours) + 
  scale_size(range = c(4,10)) +
  scale_y_discrete(limits=rev) +
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.label.y=element_blank(),
        axis.ticks.y=element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,color = "black"))

# 11.2 Fungi
################################################################################

# 11.2.1 Select and filter data

# Filter abundance table before taxa filtering
matrix_abund_fung <- matrix_abund_1[rowSums(matrix_abund_1) >= 800,]
matrix_abund_fung <- matrix_abund_fung[rowSums(matrix_abund_fung > 20) >= 2,]
# Remove reference, hybrid and problematic samples
matrix_abund_fung <- matrix_abund_fung[ ,colnames(matrix_abund_fung) %in% rownames(all_metadata_1) ]
# Remove from taxonomy table
Tax_table_strict_fung <- Tax_table_strict[ rownames(Tax_table_strict) %in% rownames(matrix_abund_fung), ]
# Only take Fungi taxa
Tax_table_strict_fung <- Tax_table_strict_fung[ Tax_table_strict_fung$phylum %in% c("Microsporidia","Mucoromycota","Oomycota","Ascomycota","Basidiomycota"), ]
# Keep only these rows in abundance table for which we have a taxonomy
matrix_abund_fung <- matrix_abund_fung[ rownames(matrix_abund_fung) %in% rownames(Tax_table_strict_fung), ]
# Remove samples with lower than 10 counts
matrix_abund_fung <- matrix_abund_fung[,colSums(matrix_abund_fung) >= 10]
# Only keep metadata of the samples with Fungi
all_metadata_fung <- all_metadata[ rownames(all_metadata) %in% colnames(matrix_abund_fung) ,]

# 11.2.2 Construct phyloseq object

# Specify datasets needed
otu_mat_fung <- matrix_abund_fung
tax_mat_fung <- Tax_table_strict_fung
samples_df_fung <- all_metadata_fung
# Transform otu table to matrix with numeric values
otu_mat_fung <- as.matrix(otu_mat_fung)
class(otu_mat_fung) <- "numeric"
# Transform taxonomy table to matrix
tax_mat_fung <- as.matrix(tax_mat_fung)
# Transform matrices to phyloseq objects
OTU_fung = phyloseq::otu_table(otu_mat_fung, taxa_are_rows = TRUE)
TAX_fung = phyloseq::tax_table(tax_mat_fung)
samples_fung = phyloseq::sample_data(samples_df_fung)
# Make phyloseq object
physeq_fung <- phyloseq::phyloseq(OTU_fung, TAX_fung, samples_fung)
physeq_fung

# 11.2.3 Define presence/absence at phylum level

# Subset to phylum level
physeq_fung_phylum <- tax_glom(physeq_fung, taxrank=rank_names(physeq_fung)[2], NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
# Need a min of 10 reads per OTU per sample
df_fung_phylum <- as.data.frame(otu_table(physeq_fung_phylum))
df_fung_phylum[df_fung_phylum < 10 ] <- 0
# Remove samples with lower than 10 counts
df_fung_phylum <- df_fung_phylum[,colSums(df_fung_phylum) >= 10]
# Only presence-abscence, to change values to 1
df_fung_phylum[df_fung_phylum > 1 ] <- 1
# Add again to Phyloseq
otu_table(physeq_fung_phylum) <- phyloseq::otu_table(df_fung_phylum, taxa_are_rows = TRUE)

# 11.2.4 Prepare and define bubble plot

# Add variable combining SpeciesID and TissueID for plot
variable1 = as.character(get_variable(physeq_fung_phylum, "SpeciesID"))
variable2 = as.character(get_variable(physeq_fung_phylum, "TissueID"))
sample_data(physeq_fung_phylum)$SpeciesID_TissueID <- mapply(paste0, variable1, variable2,collapse = "_")
# Merge per species, FG vs HG
physeq_fung_phylum_spp_tiss <- merge_samples(physeq_fung_phylum, "SpeciesID_TissueID", fun=sum) 
# Extract otu table
fung_phylum_spp_tiss = as(otu_table(physeq_fung_phylum_spp_tiss), "matrix")
# Rename TaxID's
colnames(fung_phylum_spp_tiss) <- with(Tax_table_strict_fung, phylum[match(colnames(fung_phylum_spp_tiss),rownames(Tax_table_strict_fung))])
# Melt matrix for graph
fung_phylum_spp_tiss_melt <- melt(t(fung_phylum_spp_tiss), id.vars = "SpeciesID_TissueID", variable.name = "Phylum")
fung_phylum_spp_tiss_melt$value <- as.numeric(fung_phylum_spp_tiss_melt$value)
# Put NA so only dots when there are values
fung_phylum_spp_tiss_melt[fung_phylum_spp_tiss_melt == 0] <- NA
# Add metadata
sample_data_phylum_fung <- as.data.frame(sample_data(physeq_fung_phylum))
fung_phylum_spp_tiss_melt$SpeciesID = with(sample_data_phylum_fung, SpeciesID[match(fung_phylum_spp_tiss_melt$Var2,sample_data_phylum_fung$SpeciesID_TissueID)])
fung_phylum_spp_tiss_melt$Tribe = with(sample_data_phylum_fung, Tribe[match(fung_phylum_spp_tiss_melt$Var2,sample_data_phylum_fung$SpeciesID_TissueID)])
fung_phylum_spp_tiss_melt$TissueID = with(sample_data_phylum_fung, TissueID[match(fung_phylum_spp_tiss_melt$Var2,sample_data_phylum_fung$SpeciesID_TissueID)])
# Rename FG HG
fung_phylum_spp_tiss_melt <- fung_phylum_spp_tiss_melt %>%
  mutate(TissueID = recode(TissueID, ContentFG = 'AG', ContentHG = 'PG' ))
# Combine fungsite and tissue id for graph
fung_phylum_spp_tiss_melt$para_tiss <- paste(fung_phylum_spp_tiss_melt$Var1,fung_phylum_spp_tiss_melt$TissueID)
# Order Cichlids in rigth way
fung_cichlid_phylum_host_ordering <- c("Cunlon","Asplep","Plestr","Spaery","Erecya","Trodub","Petpol","Petort","Ctehor","Astbur","Varmoo","Neotre","Neowal","Neobri","Telvit","Neonig","Altcom")
# Rename column
fung_phylum_spp_tiss_melt$FOO <- fung_phylum_spp_tiss_melt$value
# Make bubble plot
fung_plot <- ggplot(fung_phylum_spp_tiss_melt, aes(y=para_tiss,x=factor(SpeciesID, level = fung_cichlid_phylum_host_ordering),colour=Tribe, size=FOO)) +
  geom_point() +
  scale_colour_manual(values=tribe_colours) + 
  scale_size(range = c(4,10)) +
  scale_y_discrete(limits=rev) +
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.label.y=element_blank(),
        axis.ticks.y=element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,color = "black"))

# 11.3 Invertebrate parasites
################################################################################

# 11.3.1 Select and filter data

# Filter abundance table before taxa filtering
matrix_abund_para <- matrix_abund_1[rowSums(matrix_abund_1) >= 800,]
matrix_abund_para <- matrix_abund_para[rowSums(matrix_abund_para > 20) >= 2,]
# Remove reference, hybrid and problematic samples
matrix_abund_para <- matrix_abund_para[ ,colnames(matrix_abund_para) %in% rownames(all_metadata_1) ]
# Same for taxonomy table
Tax_table_strict_para <- Tax_table_strict[ rownames(Tax_table_strict) %in% rownames(matrix_abund_para), ]
# Only take Animalia parasite taxa
Tax_table_strict_para <- Tax_table_strict_para[ Tax_table_strict_para$phylum %in% c("Acanthocephala","Platyhelminthes","Nematoda","Myzozoa"), ]
# Keep only these rows in abundance table for which we have a taxonomy
matrix_abund_para <- matrix_abund_para[ rownames(matrix_abund_para) %in% rownames(Tax_table_strict_para), ]
# Remove samples with lower than 10 counts
matrix_abund_para <- matrix_abund_para[,colSums(matrix_abund_para) >= 10]
# Keep only these samples in metadata
all_metadata_para <- all_metadata[ rownames(all_metadata) %in% colnames(matrix_abund_para) ,]

# 11.3.2 Construct phyloseq object

# Specify datasets needed
otu_mat_para <- matrix_abund_para
tax_mat_para <- Tax_table_strict_para
samples_df_para <- all_metadata_para
# Transform otu table to matrix with numeric values
otu_mat_para <- as.matrix(otu_mat_para)
class(otu_mat_para) <- "numeric"
# Transform taxonomy table to matrix
tax_mat_para <- as.matrix(tax_mat_para)
# Transform matrices to phyloseq objects
OTU_para = phyloseq::otu_table(otu_mat_para, taxa_are_rows = TRUE)
TAX_para = phyloseq::tax_table(tax_mat_para)
samples_para = phyloseq::sample_data(samples_df_para)
# Make phyloseq object
physeq_para <- phyloseq::phyloseq(OTU_para, TAX_para, samples_para)
physeq_para

# 11.3.3 Define presence/absence at phylum level

# Subset to phylum level
physeq_para_phylum <- tax_glom(physeq_para, taxrank=rank_names(physeq_para)[2], NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
# Need a min of 10 reads per OTU per sample
df_para_phylum <- as.data.frame(otu_table(physeq_para_phylum))
df_para_phylum[df_para_phylum < 10 ] <- 0
# Remove samples with lower than 10 counts
df_para_phylum <- df_para_phylum[,colSums(df_para_phylum) >= 10]
# Only presence-abscence, change values to 1
df_para_phylum[df_para_phylum > 1 ] <- 1
# Add again to Phyloseq
otu_table(physeq_para_phylum) <- phyloseq::otu_table(df_para_phylum, taxa_are_rows = TRUE)

# 11.3.4 Prepare and define bubble plot at phylum level

# Add variable combining SpeciesID and TissueID for plot
variable1 = as.character(get_variable(physeq_para_phylum, "SpeciesID"))
variable2 = as.character(get_variable(physeq_para_phylum, "TissueID"))
sample_data(physeq_para_phylum)$SpeciesID_TissueID <- mapply(paste0, variable1, variable2,collapse = "_")
# Merge per species, FG vs HG
physeq_para_phylum_spp_tiss <- merge_samples(physeq_para_phylum, "SpeciesID_TissueID", fun=sum) 
# Extract otu table
para_phylum_spp_tiss = as(otu_table(physeq_para_phylum_spp_tiss), "matrix")
# Rename TaxID's
colnames(para_phylum_spp_tiss) <- with(Tax_table_strict_para, phylum[match(colnames(para_phylum_spp_tiss),rownames(Tax_table_strict_para))])
# Melt matrix for graph
library(reshape2)
para_phylum_spp_tiss_melt <- melt(t(para_phylum_spp_tiss), id.vars = "SpeciesID_TissueID", variable.name = "Phylum")
para_phylum_spp_tiss_melt$value <- as.numeric(para_phylum_spp_tiss_melt$value)
# Put NA so only dots when there are values
para_phylum_spp_tiss_melt[para_phylum_spp_tiss_melt == 0] <- NA
# Add metadata
sample_data_phylum_para <- as.data.frame(sample_data(physeq_para_phylum))
para_phylum_spp_tiss_melt$SpeciesID = with(sample_data_phylum_para, SpeciesID[match(para_phylum_spp_tiss_melt$Var2,sample_data_phylum_para$SpeciesID_TissueID)])
para_phylum_spp_tiss_melt$Tribe = with(sample_data_phylum_para, Tribe[match(para_phylum_spp_tiss_melt$Var2,sample_data_phylum_para$SpeciesID_TissueID)])
para_phylum_spp_tiss_melt$TissueID = with(sample_data_phylum_para, TissueID[match(para_phylum_spp_tiss_melt$Var2,sample_data_phylum_para$SpeciesID_TissueID)])
# Rename FG HG
para_phylum_spp_tiss_melt <- para_phylum_spp_tiss_melt %>%
  mutate(TissueID = recode(TissueID, ContentFG = 'AG', ContentHG = 'PG' ))
# Combine parasite and tissue id for graph
para_phylum_spp_tiss_melt$para_tiss <- paste(para_phylum_spp_tiss_melt$Var1,para_phylum_spp_tiss_melt$TissueID)
# Order Cichlids in rigth way
para_cichlid_phylum_host_ordering <- c("Boumic","Batgra","Batvit","Gralem","Calple","Cyafur","Xenbou","PcybrN","Cypdwj","Tromoo","Petpol","Simdia","Limdar","Simple","Gnapfe","Astbur","Lamlem")
# Rename column
para_phylum_spp_tiss_melt$FOO <- para_phylum_spp_tiss_melt$value
# Make bubble graph
para_phy_plot <- ggplot(para_phylum_spp_tiss_melt, aes(y=para_tiss,x=factor(SpeciesID, level = para_cichlid_phylum_host_ordering),colour=Tribe, size=FOO)) +
  geom_point() +
  scale_colour_manual(values=tribe_colours) + 
  scale_size(range = c(4,10)) +
  scale_y_discrete(limits=rev) +
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.label.y=element_blank(),
        axis.ticks.y=element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,color = "black"))

# 11.3.5 Define presence/absence at phylum level

# Merge at class level
physeq_para_class <- tax_glom(physeq_para, taxrank=rank_names(physeq_para)[3], NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
df_para_class <- as.data.frame(otu_table(physeq_para_class))
# Need a min of 10 reads per OTU per sample
df_para_class[df_para_class < 10 ] <- 0
# Remove samples with lower than 10 counts
df_para_class <- df_para_class[,colSums(df_para_class) >= 10]
# Only presence-abscence, to change values to 1
df_para_class[df_para_class > 1 ] <- 1
# Add again to Phyloseq
otu_table(physeq_para_class) <- phyloseq::otu_table(df_para_class, taxa_are_rows = TRUE)

# 11.3.5 Prepare and define bubble plot at class level

# Merge SpeciesID and TissueID
variable1 = as.character(get_variable(physeq_para_class, "SpeciesID"))
variable2 = as.character(get_variable(physeq_para_class, "TissueID"))
sample_data(physeq_para_class)$SpeciesID_TissueID <- mapply(paste0, variable1, variable2,collapse = "_")
# Merge per species and FG HG
physeq_para_class_spp_tiss <- merge_samples(physeq_para_class, "SpeciesID_TissueID", fun=sum) 
# combine species
para_class_spp_tiss = as(otu_table(physeq_para_class_spp_tiss), "matrix")
# Rename TaxID's
colnames(para_class_spp_tiss) <- with(Tax_table_strict_para, class[match(colnames(para_class_spp_tiss),rownames(Tax_table_strict_para))])
# Melt data
para_class_spp_tiss_melt <- melt(t(para_class_spp_tiss), id.vars = "SpeciesID_TissueID", variable.name = "class")
# Define as numeric
para_class_spp_tiss_melt$value <- as.numeric(para_class_spp_tiss_melt$value)
# Remove species that have no parasite classes
para_class_spp_tiss_melt[para_class_spp_tiss_melt == 0] <- NA
# Add metadata
sample_data_class_para <- as.data.frame(sample_data(physeq_para_class))
para_class_spp_tiss_melt$SpeciesID = with(sample_data_class_para, SpeciesID[match(para_class_spp_tiss_melt$Var2,sample_data_class_para$SpeciesID_TissueID)])
para_class_spp_tiss_melt$Tribe = with(sample_data_class_para, Tribe[match(para_class_spp_tiss_melt$Var2,sample_data_class_para$SpeciesID_TissueID)])
para_class_spp_tiss_melt$TissueID = with(sample_data_class_para, TissueID[match(para_class_spp_tiss_melt$Var2,sample_data_class_para$SpeciesID_TissueID)])
# Rename FG HG
para_class_spp_tiss_melt <- para_class_spp_tiss_melt %>%
  mutate(TissueID = recode(TissueID, ContentFG = 'AG', ContentHG = 'PG' ))
# Combine metadata
para_class_spp_tiss_melt$para_tiss <- paste(para_class_spp_tiss_melt$Var1,para_class_spp_tiss_melt$TissueID)
# Order the taxa in the desired way
para_class_ordering <- c("Eoacanthocephala AG","Eoacanthocephala PG","Chromadorea AG","Chromadorea PG","Cestoda AG","Cestoda PG","Trematoda AG","Trematoda PG")
para_class_spp_tiss_melt$para_tiss <- factor(para_class_spp_tiss_melt$para_tiss, levels=para_class_ordering)
para_class_spp_tiss_melt <- para_class_spp_tiss_melt[order(para_class_spp_tiss_melt$para_tiss),]
# Rename column
para_class_spp_tiss_melt$FOO <- para_class_spp_tiss_melt$value
# Make bubble plot
para_class_plot <- ggplot(para_class_spp_tiss_melt, aes(y=para_tiss,x=factor(SpeciesID, level = para_cichlid_phylum_host_ordering),colour=Tribe, size=FOO)) +
  geom_point() +
  scale_colour_manual(values=tribe_colours) + 
  theme_minimal() + 
  scale_size(range = c(4,10)) +
  scale_y_discrete(limits=rev) +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,color = "black"))

# 11.4 Combined non-diet taxa
################################################################################

# 11.4.1 Make combined bubble plot

pdf(file = "Plot_NonDiet.pdf",
    width = 10,
    height = 10) 
grid.arrange(prot_plot, fung_plot, para_phy_plot, para_class_plot, ncol=1 )
dev.off()








