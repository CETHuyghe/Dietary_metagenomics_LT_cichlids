###########################################################
# 08. Genus level analysis
###########################################################
# RScript to analyze the metagenomic diet data at genus level 

# 8.1 Prepare data to genus level analyses
################################################################################

# 8.1.1 Subset phyloseq object

# Make physeq with identification to genus level
physeq_genus <- tax_glom(physeq, taxrank=rank_names(physeq)[6], NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
df_physeq_genus <- as.data.frame(otu_table(physeq_genus))

# 8.1.2 Filter data

# Remove hits with fewer than 20 reads
df_physeq_genus[(df_physeq_genus) < 20 ] <- 0
# Remove genera with fewer than 100 reads
df_physeq_genus <- df_physeq_genus[rowSums(df_physeq_genus) >= 100,]
# Remove hits with fewer than 0.5 % abundance in sample
df_physeq_genus[apply(df_physeq_genus,2,function(x){x/sum(x)}) < 0.005 ] <- 0
# Add filtered table to phyloseq
otu_table(physeq_genus) <- phyloseq::otu_table(df_physeq_genus, taxa_are_rows = TRUE)
# Remove samples with fewer than 100 reads
physeq_genus = prune_samples(sample_sums(physeq_genus)>=100, physeq_genus)

# 8.2 Investigate data
################################################################################

# 8.2.1 Check data

physeq_genus
# Check number of cichlid species
length(unique(sample_data(physeq_genus)$SpeciesID))
# Check min, max, mean, median and sd in number of reads per sample
summary(colSums(otu_table(physeq_genus)))
sd(colSums(otu_table(physeq_genus)))

# 8.3 Diet genus composition
################################################################################

# 8.3.1 FOO

# Merge samples by Specimen
physeq_genus_Specimen0 <- physeq_genus
physeq_genus_Specimen <- merge_samples(physeq_genus_Specimen0,"SpecimenID",fun=sum)
# Take relative value
physeq_genus_Specimen1 <- phyloseq::transform_sample_counts(physeq_genus_Specimen, function(x) (x / sum(x)*100))
# Merge by species
sample_data(physeq_genus_Specimen1)$SpeciesID <-with(all_metadata_2, SpeciesID[match(rownames(sample_data(physeq_genus_Specimen1)),all_metadata_2$SpecimenID)])
physeq_genus_Specimen_spp <- merge_samples(physeq_genus_Specimen1, "SpeciesID", fun=mean) 
# Take relative value
physeq_genus_Specimen_spp_rel <- phyloseq::transform_sample_counts(physeq_genus_Specimen_spp, function(x) (x / sum(x)*100))
genus_Specimen_spp_rel = as(otu_table(physeq_genus_Specimen_spp_rel), "matrix")
# Rename TaxID's
colnames(genus_Specimen_spp_rel) <- with(Tax_table_strict_2, genus[match(colnames(genus_Specimen_spp_rel),rownames(Tax_table_strict_2))])
# Check overall abundance per genus
colSums(genus_Specimen_spp_rel)/nrow(genus_Specimen_spp_rel)
# Calculate FOO
genus_Specimen_spp_rel2 <- genus_Specimen_spp_rel
genus_Specimen_spp_rel2[genus_Specimen_spp_rel2 >= 1] <- 1
genus_Specimen_spp_rel2[genus_Specimen_spp_rel2 < 1] <- 0
# FOO
colSums(genus_Specimen_spp_rel2)
# FOO%
colSums(genus_Specimen_spp_rel2)/nrow(genus_Specimen_spp_rel2)


