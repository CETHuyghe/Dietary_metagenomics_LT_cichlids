###########################################################
# 02. Tribe level analysis Focus dataset
###########################################################
# RScript to analyze the metagenomic cichlid diet data at tribe level 

# 2.1 Prepare data to tribe level
################################################################################

# 2.1.1 Subset phyloseq object

# Make physeq with identification to tribe level
physeq_tribe_cichl <- tax_glom(physeq, taxrank=rank_names(physeq)[6], NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
df_physeq_tribe_cichl <- as.data.frame(otu_table(physeq_tribe_cichl))

# 2.1.2 Filter data

# Remove hits with fewer than 20 reads
df_physeq_tribe_cichl[(df_physeq_tribe_cichl) < 20 ] <- 0
# Remove taxa with fewer than 20 reads
df_physeq_tribe_cichl <- df_physeq_tribe_cichl[rowSums(df_physeq_tribe_cichl) >= 20,]
# Remove hits with fewer than 0.5 % abundance in sample
df_physeq_tribe_cichl[apply(df_physeq_tribe_cichl,2,function(x){x/sum(x)}) < 0.005 ] <- 0
# Add to phyloseq object
otu_table(physeq_tribe_cichl) <- phyloseq::otu_table(df_physeq_tribe_cichl, taxa_are_rows = TRUE)
# Remove Samples that have less than 100 reads left
physeq_tribe_cichl = prune_samples(sample_sums(physeq_tribe_cichl)>=20, physeq_tribe_cichl)

# 2.1.3 Check data

# Check number of species with cichlid DNA
summary(unique(sample_data(physeq_tribe_cichl)$SpeciesID))
# Check nr reads per sample
summary(colSums(otu_table(physeq_tribe_cichl)))
sd(colSums(otu_table(physeq_tribe_cichl)))

# 2.2 Occurrence foreign cichlid DNA
################################################################################

# 2.2.1 Combine specimen diet abundance

# Merge samples by Specimen
physeq_tribe_cichl_Specimen0 <- physeq_tribe_cichl
physeq_tribe_cichl_Specimen <- merge_samples(physeq_tribe_cichl_Specimen0,"SpecimenID",fun=sum)
# Make dataset for total nr reads
physeq_tribe_cichl_Specimen0_nr <- as.data.frame(t(otu_table(physeq_tribe_cichl_Specimen0)))
physeq_tribe_cichl_Specimen0_nr$rowsum <- rowSums(physeq_tribe_cichl_Specimen0_nr)
# Add to metadata
all_metadata_2$rowsum_tribe_cichl <-with(physeq_tribe_cichl_Specimen0_nr, rowsum[match(rownames(all_metadata_2),rownames(physeq_tribe_cichl_Specimen0_nr))])

# 2.2.2 Calculate FOO%

#Take otu table per specimen
tribe_cichl_Specimen = as(otu_table(physeq_tribe_cichl_Specimen), "matrix")
# Make copy where we keep original TaxIDs
tribe_cichl_Specimen0 <- as.data.frame(t(tribe_cichl_Specimen))
# Take FOO per specimen to put back into PhyloSeq object
tribe_cichl_Specimen0[apply(tribe_cichl_Specimen0,2,function(x){x/sum(x)}) >= 0.005] <- 1
tribe_cichl_Specimen0[tribe_cichl_Specimen0 > 1] <- 0
rowSums(tribe_cichl_Specimen0)
# Make Phyloseq object with FOO
physeq_tribe_cichl_FOO_Specimen <- physeq_tribe_cichl_Specimen
otu_table(physeq_tribe_cichl_FOO_Specimen) <- phyloseq::otu_table(tribe_cichl_Specimen0, taxa_are_rows = TRUE)
sample_data(physeq_tribe_cichl_FOO_Specimen)$SpeciesID <- with(all_metadata_2, SpeciesID[match(rownames(sample_data(physeq_tribe_cichl_FOO_Specimen)),all_metadata_2$SpecimenID)])
# Merge specimens per species
physeq_tribe_cichl_FOO_Specimen_spp <- merge_samples(physeq_tribe_cichl_FOO_Specimen, "SpeciesID", fun=sum) 
# Extract otu_table
tribe_cichl_Specimen_spp = as(otu_table(physeq_tribe_cichl_FOO_Specimen_spp), "matrix")
# rename TaxIDs
colnames(tribe_cichl_Specimen_spp) <- with(Tax_table_strict_2, tribe[match(colnames(tribe_cichl_Specimen_spp),rownames(Tax_table_strict_2))])
# Get nr specimens per species to calculate %
physeq_tribe_cichl_FOO_Specimen_sppFreq = as(sample_data(physeq_tribe_cichl_FOO_Specimen), "data.frame")
tribe_cichl_specimen_freq <- table(physeq_tribe_cichl_FOO_Specimen_sppFreq$SpeciesID)
# Calculate FOO%
tribe_cichl_Specimen_spp_FOOPerc <- t(apply(tribe_cichl_Specimen_spp, 2, "/", tribe_cichl_specimen_freq))

# 2.2.2 Plot FOO%
library(reshape2)
library(ggplot2)
library(ape)

# Prepare data
bar_tribe_cichl_sp_FOO <- melt(t(tribe_cichl_Specimen_spp), id.vars = "SpeciesID", variable.name = "tribe_cichl")
bar_tribe_cichl_sp_FOO$value <- as.numeric(bar_tribe_cichl_sp_FOO$value)
bar_tribe_cichl_sp_FOO$Var2 <- factor(bar_tribe_cichl_sp_FOO$Var2)
bar_tribe_cichl_sp_FOO$tribe_cichl <- factor(bar_tribe_cichl_sp_FOO$Var1)
bar_tribe_cichl_sp_FOO$Tribe = with(all_metadata_2, Tribe[match(bar_tribe_cichl_sp_FOO$Var2,all_metadata_2$SpeciesID)])
bar_tribe_cichl_sp_FOO$Tribe <- factor(bar_tribe_cichl_sp_FOO$Tribe)
# Remove tribe of host cichlid
bar_tribe_cichl_sp_FOO$value <- replace(bar_tribe_cichl_sp_FOO$value, which(bar_tribe_cichl_sp_FOO$Var1 == bar_tribe_cichl_sp_FOO$Tribe), NA)
bar_tribe_cichl_sp_FOO$value[bar_tribe_cichl_sp_FOO$value == 0] <- NA
bar_tribe_cichl_sp_FOO$value[bar_tribe_cichl_sp_FOO$Var1 == "Haplochromini" & bar_tribe_cichl_sp_FOO$Tribe == "Tropheini"] <- NA
bar_tribe_cichl_sp_FOO$value[bar_tribe_cichl_sp_FOO$Var1 == "Tropheini" & bar_tribe_cichl_sp_FOO$Tribe == "Haplochromini"] <- NA
# Define tribe plot colors
tribe_colours <- c("Bathybatini"="#242626", "Benthochromini"="#AE262A", "Boulengerochromini"="#59595C", "Cyphotilapiini"="#FDDF13","Cyprichromini"="#F04D29","Ectodini"="#9AB9D9","Eretmodini"="#682E7A","Haplochromini"="darkgreen","Lamprologini"="#C588BB","Limnochromini"="#535CA9", "Perissodini"="orange","Trematocarini"="#959170", "Tropheini"="#86C773","Oreochromini"="grey")
# Order cichlid species by tree position
# Load tree
ActualPhylogeny = "b1_with_Oretan.tre"
pruned.tree = ape::read.tree(ActualPhylogeny)
# Filter tree by species in this dataset
spp <- all_metadata_2[order(all_metadata_2$Tribe),]
species <- unique(spp$SpeciesID)
pruned.tree = drop.tip(pruned.tree, pruned.tree$tip.label[-which(pruned.tree$tip.label %in% species)])
# Define order of consumed cichlid tribes
Tribe_order <- c("Oreochromini","Boulengerochromini","Trematocarini","Bathybatini","Cyphotilapiini","Limnochromini","Ectodini","Cyprichromini","Perissodini","Benthochromini","Eretmodini","Tropheini","Haplochromini","Lamprologini")
# Put NA so only dots when there are values
bar_tribe_cichl_sp_FOO[bar_tribe_cichl_sp_FOO == 0] <- NA
# Make bubble plot
pdf(file = "Plot_Cichlid_FOO_24092024.pdf",
    width = 21,
    height = 5) 
ggplot(bar_tribe_cichl_sp_FOO, aes(y=factor(Var1,level=Tribe_order),x=factor(Var2, level =pruned.tree$tip.label),colour=Var1, size=value)) +
  geom_point() +
  scale_colour_manual(values=tribe_colours) + 
  scale_size(range = c(4,10)) +
  scale_x_discrete(limits=rev) +
  scale_y_discrete(limits=rev) +
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.label.y=element_blank(),
        axis.ticks.y=element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,color = "black",size=12),
        axis.text.y = element_text( vjust = 0.5, color = "black",size=12),
        legend.text = element_text(size=12),color = "black") 
dev.off()

