###########################################################
# 09. Species level analysis
###########################################################
# RScript to analyze the metagenomic diet data at species level 

# 9.1 Prepare data to species level analyses
################################################################################

# 9.1.1 Subset phyloseq object

# Make physeq with identification to species level
physeq_species <- tax_glom(physeq, taxrank=rank_names(physeq)[7], NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
df_physeq_species <- as.data.frame(otu_table(physeq_species))

# 9.1.2 Filter data

# Need at least 20 reads for a hit
df_physeq_species[(df_physeq_species) < 20 ] <- 0
# Need a min of 100 total reads per diet species
df_physeq_species <- df_physeq_species[rowSums(df_physeq_species) >= 100,]
# Need 0.5 % abundance reads per diet species per sample for hit
df_physeq_species[apply(df_physeq_species,2,function(x){x/sum(x)}) < 0.005 ] <- 0
# Add to physeq_species
otu_table(physeq_species) <- phyloseq::otu_table(df_physeq_species, taxa_are_rows = TRUE)
# Remove Samples that have less than 100 reads left
physeq_species = prune_samples(sample_sums(physeq_species)>=100, physeq_species)

# 9.2 Investigate data
################################################################################

# 9.2.1 Check data

physeq_species
# Check number of cichlid species
length(unique(sample_data(physeq_species)$SpeciesID))
# Check min, max, mean, median and sd in number of reads per sample
summary(colSums(otu_table(physeq_species)))
sd(colSums(otu_table(physeq_species)))

# 9.3 Diet genus composition
################################################################################

# 9.3.1 FOO

# Merge samples by Specimen
physeq_species_Specimen0 <- physeq_species
physeq_species_Specimen <- merge_samples(physeq_species_Specimen0,"SpecimenID",fun=sum)
# Take relative values
physeq_species_Specimen1 <- phyloseq::transform_sample_counts(physeq_species_Specimen, function(x) (x / sum(x)*100))
# Merge by cichlid species
sample_data(physeq_species_Specimen1)$SpeciesID <-with(all_metadata_2, SpeciesID[match(rownames(sample_data(physeq_species_Specimen1)),all_metadata_2$SpecimenID)])
physeq_species_Specimen_spp <- merge_samples(physeq_species_Specimen1, "SpeciesID", fun=mean) 
# Take relative values
physeq_species_Specimen_spp_rel <- phyloseq::transform_sample_counts(physeq_species_Specimen_spp, function(x) (x / sum(x)*100))
species_Specimen_spp_rel = as(otu_table(physeq_species_Specimen_spp_rel), "matrix")
# Rename TaxID's
colnames(species_Specimen_spp_rel) <- with(Tax_table_strict_2, species[match(colnames(species_Specimen_spp_rel),rownames(Tax_table_strict_2))])
# Check abundance per species
colSums(species_Specimen_spp_rel)/nrow(species_Specimen_spp_rel)

