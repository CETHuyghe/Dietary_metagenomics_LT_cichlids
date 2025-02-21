###########################################################
# 07. Family level analysis
###########################################################
# RScript to analyze the metagenomic diet data at family level 

# 7.1 Prepare data to family level analyses
################################################################################

# 7.1.1 Subset phyloseq object

# Make physeq with identification to family level
physeq_family <- tax_glom(physeq, taxrank=rank_names(physeq)[5], NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
df_physeq_family <- as.data.frame(otu_table(physeq_family))

# 7.1.2 Filter data

# Hits with fewer than 20 reads removed
df_physeq_family[(df_physeq_family) < 20 ] <- 0
# Remove families 
df_physeq_family <- df_physeq_family[rowSums(df_physeq_family) >= 200,]
# Remove hits with fewer than 0.5 % abundance in sample
df_physeq_family[apply(df_physeq_family,2,function(x){x/sum(x)}) < 0.005 ] <- 0
# Add filtered table to physeq_family
otu_table(physeq_family) <- phyloseq::otu_table(df_physeq_family, taxa_are_rows = TRUE)
# Remove Samples that have fewer than 100 reads left
physeq_family = prune_samples(sample_sums(physeq_family)>=100, physeq_family)
physeq_family

# 7.2 Investigate data and calculate diversity
################################################################################

# 7.2.1 Check data

# Check number of cichlid species
length(unique(sample_data(physeq_family)$SpeciesID))
# Check min, max, mean, median and sd reads per sample
summary(colSums(otu_table(physeq_family)))
sd(colSums(otu_table(physeq_family)))

# 7.3 Diet family composition
################################################################################

# 7.3.1 FOO

# Merge samples by Specimen
physeq_family_Specimen <- merge_samples(physeq_family,"SpecimenID",fun=sum)
# Take rel
physeq_family_Specimen_spp_rel <- phyloseq::transform_sample_counts(physeq_family_Specimen_spp, function(x) (x / sum(x)*100))
family_Specimen_spp_rel = as(otu_table(physeq_family_Specimen_spp_rel), "matrix")
# Rename TaxID's
colnames(family_Specimen_spp_rel) <- with(Tax_table_strict_2, family[match(colnames(family_Specimen_spp_rel),rownames(Tax_table_strict_2))])
# Check overall abundance per family
colSums(family_Specimen_spp_rel)/nrow(family_Specimen_spp_rel)
# Check FOO
family_Specimen_spp_rel2 <- family_Specimen_spp_rel
family_Specimen_spp_rel2[family_Specimen_spp_rel2 >= 1] <- 1
family_Specimen_spp_rel2[family_Specimen_spp_rel2 < 1] <- 0
# FOO
colSums(family_Specimen_spp_rel2)
# FOO%
colSums(family_Specimen_spp_rel2)/nrow(family_Specimen_spp_rel2)



