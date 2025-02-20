###########################################################
# 02. Phylum level initial analysis
###########################################################
# RScript to analyze the metagenomic diet data at phylum level 

# 2.1 Prepare data to phylum level
################################################################################

# 2.1.1 Subset phyloseq object

# Subset phyloseq object at phylum level
physeq_phylum <- tax_glom(physeq, taxrank=rank_names(physeq)[2], NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
# Filter otu table at phylum level
df_physeq_phylum <- as.data.frame(otu_table(physeq_phylum))

# 2.1.2 Filter data

# Remove hits with fewer than 20 counts
df_physeq_phylum[(df_physeq_phylum) < 20 ] <- 0
# Remove phyla with fewer than 500 total counts
df_physeq_phylum <- df_physeq_phylum[rowSums(df_physeq_phylum) >= 500,]
# Remove hits with fewer than 0.5 % count abundance in a sample
df_physeq_phylum[apply(df_physeq_phylum,2,function(x){x/sum(x)}) < 0.005 ] <- 0
# Add again to the phyloseq object
otu_table(physeq_phylum) <- phyloseq::otu_table(df_physeq_phylum, taxa_are_rows = TRUE)
# Remove Samples that have fewer than 100 counts left
physeq_phylum = prune_samples(sample_sums(physeq_phylum)>=100, physeq_phylum)
physeq_phylum

# 2.2 Calculate number of reads and diversity
################################################################################

# 2.2.1 Number of total reads

# Investigate min, max, mean, median and standard deviation in number of reads per sample
summary(colSums(otu_table(physeq_phylum)))
sd(colSums(otu_table(physeq_phylum)))

# 2.2.2 Diversity Plot

# Plot the Shannon and Chao1 diversity indices
div_phl_plot <- plot_richness(physeq_phylum, x='SpeciesID', measures=c("Shannon","Chao1"), color="Tribe") +
  geom_point(size=2) +
  geom_boxplot(aes(fill=Tribe),alpha=0.7,size=0.5 ) +
  theme_minimal() +
  theme(axis.title=element_text(size=10),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=10),
        axis.text.y = element_text( size=8),
        legend.position="none")+
  scale_colour_manual(values=tribe_colours) +
  scale_fill_manual(values=tribe_colours) + 
  scale_x_discrete(limits=rev(pruned.tree$tip.label))

# 2.3 Diet phylum abundance 
################################################################################

# 2.3.1 Abundance barplots with total number of reads per sample

# Define colours
phyla_cols = c("Annelida"="lightpink3","Arthropoda"="peru", "Bacillariophyta"="darkseagreen2", "Chlorophyta"="palegreen3","Algae"= "#2f7d44" , "Chordata"="#224e66","Cnidaria"="ghostwhite","Mollusca"="brown","Nematoda"="indianred4","Platyhelminthes"="lightpink3","Porifera"="khaki1","Streptophyta"="mediumseagreen","Myzozoa"="gray38","Rotifera"="lightsalmon4")
# Define ordering of species and tribe by order in the phylogenetic tree
Species_ordering <- rev(pruned.tree$tip.label)
Tribe_ordering <-unique(with(all_metadata_2, Tribe[match(Species_ordering,all_metadata_2$SpeciesID)]))
# Stacked barplot of absolute read counts
pdf(file = "Plot_AbsBar_Phylum_19092024.pdf",
    width = 21,
    height = 10) 
phyloseq::plot_bar(physeq_phylum, fill="phylum") + 
  geom_bar(aes(fill = phylum), stat="identity", position="stack") +
  labs(x = "", y = "Total number of counts") +
  theme(axis.title=element_text(size=12, color = "black"),
        strip.text.x = element_text(),
        strip.background = element_rect(color="grey", fill="white", size=1.5, linetype="solid")) +
  facet_wrap(factor(Tribe, level=Tribe_ordering)~factor(SpeciesID, level=Species_ordering), scales=("free"),nrow=5) +
  scale_x_discrete(label=NULL)  + 
  theme(panel.background = element_blank()) +
  scale_fill_manual(values=phyla_cols)+
  xlab("Sample")
dev.off()

# 2.3.2 Main Figure 1: diet abundance per species

# 2.3.2.1 Merge samples per species

# Merge samples by Specimen
# Take the sum of samples belonging to the same individual
physeq_phylum_Specimen <- merge_samples(physeq_phylum,"SpecimenID",fun=sum)
# Then merge by species
# Take relative abundance for equal contribution each individual
physeq_phylum_Specimen1 <- phyloseq::transform_sample_counts(physeq_phylum_Specimen, function(x) (x / sum(x)*100))
# Add SpeciesID to sample data
sample_data(physeq_phylum_Specimen1)$SpeciesID <-with(all_metadata_2, SpeciesID[match(rownames(sample_data(physeq_phylum_Specimen1)),all_metadata_2$SpecimenID)])
# Merge by SpeciesID
physeq_phylum_Specimen_spp <- merge_samples(physeq_phylum_Specimen1, "SpeciesID", fun=mean) 
# Add metadata
sample_data(physeq_phylum_Specimen_spp)$d15N_mean <-with(all_metadata_2, d15N_mean[match(rownames(sample_data(physeq_phylum_Specimen_spp)),all_metadata_2$SpeciesID)])
sample_data(physeq_phylum_Specimen_spp)$d15N_std_dev <-with(all_metadata_2, d15N_std_dev[match(rownames(sample_data(physeq_phylum_Specimen_spp)),all_metadata_2$SpeciesID)])
sample_data(physeq_phylum_Specimen_spp)$Sampling <-with(all_metadata_2, Sampling[match(rownames(sample_data(physeq_phylum_Specimen_spp)),all_metadata_2$SpeciesID)])
sample_data(physeq_phylum_Specimen_spp)$Tribe <-with(all_metadata_2, Tribe[match(rownames(sample_data(physeq_phylum_Specimen_spp)),all_metadata_2$SpeciesID)])
sample_data(physeq_phylum_Specimen_spp)$BreedingType <-with(all_metadata_2, BreedingType[match(rownames(sample_data(physeq_phylum_Specimen_spp)),all_metadata_2$SpeciesID)])
sample_data(physeq_phylum_Specimen_spp)$BreedingMode <-with(all_metadata_2, BreedingMode[match(rownames(sample_data(physeq_phylum_Specimen_spp)),all_metadata_2$SpeciesID)])
sample_data(physeq_phylum_Specimen_spp)$d13C_mean <-with(all_metadata_2, d13C_mean[match(rownames(sample_data(physeq_phylum_Specimen_spp)),all_metadata_2$SpeciesID)])
sample_data(physeq_phylum_Specimen_spp)$d13C_std_dev <-with(all_metadata_2, d13C_std_dev[match(rownames(sample_data(physeq_phylum_Specimen_spp)),all_metadata_2$SpeciesID)])
sample_data(physeq_phylum_Specimen_spp)$FoodCat <-with(all_metadata_2, FoodCat[match(rownames(sample_data(physeq_phylum_Specimen_spp)),all_metadata_2$SpeciesID)])
sample_data(physeq_phylum_Specimen_spp)$Habitat <-with(all_metadata_2, Habitat[match(rownames(sample_data(physeq_phylum_Specimen_spp)),all_metadata_2$SpeciesID)])
# Take relative values
physeq_phylum_Specimen_spp_rel <- phyloseq::transform_sample_counts(physeq_phylum_Specimen_spp, function(x) (x / sum(x)*100))

# 2.3.2.2 Prepare data for stacked barplot

# Extract otu_table from Phyloseq object as count table
phylum_Specimen_spp_rel = as(otu_table(physeq_phylum_Specimen_spp_rel), "matrix")
# Rename TaxID's by their phylum
colnames(phylum_Specimen_spp_rel) <- with(Tax_table_strict_2, phylum[match(colnames(phylum_Specimen_spp_rel),rownames(Tax_table_strict_2))])
# Melt count table for stacking 
bar_phylum_sp_otu <- melt(phylum_Specimen_spp_rel, id.vars = "SpeciesID", variable.name = "Phylum")
# Define numeric and factor
bar_phylum_sp_otu$value <- as.numeric(bar_phylum_sp_otu$value)
bar_phylum_sp_otu$Var2 <- factor(bar_phylum_sp_otu$Var2)
bar_phylum_sp_otu$Var1 <- factor(bar_phylum_sp_otu$Var1)
# Make table with the metadata of the species
Metadata_spp <- all_metadata_2[c("SpeciesID", "Tribe","BreedingType","d15N_mean", "d13C_mean", "FoodCat","Habitat")]
Metadata_spp <- Metadata_spp[!duplicated(Metadata_spp), ]
rownames(Metadata_spp) <- Metadata_spp$SpeciesID
Metadata_spp$d15N_mean <- as.numeric(Metadata_spp$d15N_mean)
Metadata_spp$d13C_mean <- as.numeric(Metadata_spp$d13C_mean)

# 2.3.2.3 Add additional data for plot

# Add table with the phyla identified as a food source by the literature
Phyla_Lit <- read.delim(file="../Lit_Diet_Phylum_08022024.txt",sep="\t",check.names=FALSE, row.names = 1)
# Melt the literature table
Phyla_Lit_melt <- melt(t(Phyla_Lit), id.vars = "SpeciesID", variable.name = "Phylum")
# Define numeric
Phyla_Lit_melt$value <- as.numeric(Phyla_Lit_melt$value)

# 2.3.2.4 Define phylogenetic tree figure

# Make basis figure with phylogenetic tree
p <- ggtree(pruned.tree, ladderize = FALSE)
# Assign cichlid tribes to tips tree
p$data$Tribe <- with(Metadata_spp, Tribe[match(p$data$label,rownames(Metadata_spp))])
# Check tribes in tree
p + geom_tree(size=1.5,aes(colour=Tribe)) + scale_colour_manual(values=tribe_colours) + geom_text(aes(label=node), hjust=-.3) + geom_treescale()
# Define tribes of other branches
p$data$Tribe[p$data$parent >= 85 & p$data$parent < 94 ] <- "Tropheini"
p$data$Tribe[p$data$parent == 84 & p$data$node == 85 ] <- "Tropheini"
p$data$Tribe[p$data$parent == 83 & p$data$node == 95 ] <- "Eretmodini"
p$data$Tribe[p$data$parent == 96 & p$data$node == 98 ] <- "Cyprichromini"
p$data$Tribe[p$data$parent >= 101 & p$data$parent < 108 ] <- "Ectodini"
p$data$Tribe[p$data$parent == 100 & p$data$node == 109 ] <- "Limnochromini"
p$data$Tribe[p$data$parent == 100 & p$data$node == 101 ] <- "Ectodini"
p$data$Tribe[p$data$parent == 99 & p$data$node == 110 ] <- "Cyphotilapiini"
p$data$Tribe[p$data$parent >= 62 & p$data$parent < 80 ] <- "Lamprologini"
p$data$Tribe[p$data$parent == 61 & p$data$node == 62 ] <- "Lamprologini"
p$data$Tribe[p$data$parent == 112 & p$data$node == 114 ] <- "Trematocarini"
p$data$Tribe[p$data$parent == 114 & p$data$node == 115 ] <- "Trematocarini"
p$data$Tribe[p$data$parent == 112 & p$data$node == 113 ] <- "Bathybatini"

# 2.3.2.5 Construct main figure
library(ggimage)

pdf(file = "Plot_abund_phyla_23092024.pdf",   # The directory you want to save the file in
    width = 20, # The width of the plot in inches
    height = 10) # The height of the plot in inches
p + geom_tree(size=1.3,aes(colour=Tribe)) + 
  scale_colour_manual(values=tribe_colours) + 
  geom_treescale()+
  geom_tiplab(geom="text",offset=0,size=3,as_ylab=FALSE) +
  geom_fruit(data=bar_phylum_sp_otu, geom=geom_bar,mapping=aes(y=Var1, x=value,fill = Var2), stat="identity", width = 0.9,pwidth=2,offset=0.7)+
  scale_fill_manual(values=phyla_cols) + 
  new_scale_colour() +
  geom_fruit(data=Phyla_Lit_melt,geom=geom_point,mapping=aes(y=Var1,x=Var2,colour=Var2, size=value), pwidth = 0.3, offset=-1.24,grid.params=list()) +
  scale_colour_manual(values=phyla_cols) +
  scale_size_continuous(range = c(-1, 5))  + 
  new_scale_fill() +
  geom_fruit(data=Metadata_spp,geom=geom_tile,mapping=aes(y=SpeciesID,fill=d13C_mean),width = 0.6, offset=-2.1) + 
  scale_fill_distiller(direction=1) + new_scale_fill() + 
  new_scale_fill() +
  geom_fruit(data=Metadata_spp,geom=geom_tile,mapping=aes(y=SpeciesID,fill=d15N_mean),width = 0.6, offset=0.045) + 
  scale_fill_distiller(direction=1,palette='YlOrBr') +  
  geom_tiplab(aes(image=paste0("../paintings_fish/",label,'.png')),geom="image",offset=2.25,size=.03)  
dev.off()

# 2.3.3 Investigate the abundance of phyla over all species

# Make dataset with most abundant phylum in each species
main_phylum <- as.data.frame(colnames(phylum_Specimen_spp_rel)[apply(phylum_Specimen_spp_rel,1,which.max)])
# Modify names
rownames(main_phylum) <- rownames(phylum_Specimen_spp_rel)
main_phylum$phylum <- main_phylum$`colnames(phylum_Specimen_spp_rel)[apply(phylum_Specimen_spp_rel, 1, which.max)]`
# Check which phyla are main
unique(main_phylum$phylum)
# Check main phyla in how many species
length(which(main_phylum$phylum == "Arthropoda"))
length(which(main_phylum$phylum == "Bacillariophyta"))
length(which(main_phylum$phylum == "Chordata"))
length(which(main_phylum$phylum == "Streptophyta"))
length(which(main_phylum$phylum == "Chlorophyta"))
length(which(main_phylum$phylum == "Cnidaria"))
# Check the relative abundance of phyla over all species
colSums(phylum_Specimen_spp_rel)/nrow(phylum_Specimen_spp_rel)
# Investigate occurrence of phyla in each species
phylum_Specimen_spp_rel2 <- phylum_Specimen_spp_rel
# Defined to occur in species if present more than 0.5 %
phylum_Specimen_spp_rel2[phylum_Specimen_spp_rel2 >= 0.5] <- 1
phylum_Specimen_spp_rel2[phylum_Specimen_spp_rel2 < 0.5] <- 0
# Check in how many species the phyla occur
colSums(phylum_Specimen_spp_rel2)
# Check the occurence of phyla in species over all species
colSums(phylum_Specimen_spp_rel2)/nrow(phylum_Specimen_spp_rel2)

# 2.4 The Frequency of Occurence (FOO%)
################################################################################

# 2.4.1 Calculate the total occurrence over all specimens

# Take otu table per specimen
phylum_Specimen = as(otu_table(physeq_phylum_Specimen), "matrix")
# Make copy where we keep original TaxIDs
phylum_Specimen0 <- phylum_Specimen
# Rename TaxID's
colnames(phylum_Specimen) <- with(Tax_table_strict_2, phylum[match(colnames(phylum_Specimen),rownames(Tax_table_strict_2))])
# Make dataset with main phylum in each individual
main_phylum_specimen <- as.data.frame(colnames(phylum_Specimen)[apply(phylum_Specimen,1,which.max)])
rownames(main_phylum_specimen) <- rownames(phylum_Specimen)
main_phylum_specimen$phylum <- main_phylum_specimen$`colnames(phylum_Specimen)[apply(phylum_Specimen, 1, which.max)]`
# Look at taxa that are the main food source over all individuals
unique(main_phylum_specimen$phylum)
# Investigate these phyla further
length(which(main_phylum_specimen$phylum == "Arthropoda"))
length(which(main_phylum_specimen$phylum == "Bacillariophyta"))
length(which(main_phylum_specimen$phylum == "Chordata"))
length(which(main_phylum_specimen$phylum == "Streptophyta"))
length(which(main_phylum_specimen$phylum == "Chlorophyta"))
length(which(main_phylum_specimen$phylum == "Cnidaria"))
length(which(main_phylum_specimen$phylum == "Annelida"))
length(which(main_phylum_specimen$phylum == "Mollusca"))
# Look at the Frequency of Occurrence (FOO) per specimen
phylum_Specimen2 <- as.data.frame(t(phylum_Specimen))
# Assign present if more than 0.5 % abundance
phylum_Specimen2[apply(phylum_Specimen2,2,function(x){x/sum(x)}) >= 0.005] <- 1
phylum_Specimen2[phylum_Specimen2 > 1] <- 0

# 2.4.2 Calculate the FOO% of diet phyla over all specimens

# Take the relative FOO over all individuals in this study
rowSums(phylum_Specimen2)/193
# Take FOO per specimen to put back into PhyloSeq object
phylum_Specimen0 <- as.data.frame(t(phylum_Specimen0))
phylum_Specimen0[apply(phylum_Specimen0,2,function(x){x/sum(x)}) >= 0.005] <- 1
phylum_Specimen0[phylum_Specimen0 > 1] <- 0
rowSums(phylum_Specimen0)

# 2.4.3 Calculate the FOO% over all species

# Make Phyloseq object with FOO
physeq_phylum_FOO_Specimen <- physeq_phylum_Specimen
otu_table(physeq_phylum_FOO_Specimen) <- phyloseq::otu_table(phylum_Specimen0, taxa_are_rows = TRUE)
sample_data(physeq_phylum_FOO_Specimen)$SpeciesID <- with(all_metadata_2, SpeciesID[match(rownames(sample_data(physeq_phylum_FOO_Specimen)),all_metadata_2$SpecimenID)])
# Merge specimens with FOO per species by taking sum
physeq_phylum_FOO_Specimen_spp <- merge_samples(physeq_phylum_FOO_Specimen, "SpeciesID", fun=sum) 
# Extract otu_table
phylum_Specimen_spp = as(otu_table(physeq_phylum_FOO_Specimen_spp), "matrix")
# rename TaxIDs
colnames(phylum_Specimen_spp) <- with(Tax_table_strict_2, phylum[match(colnames(phylum_Specimen_spp),rownames(Tax_table_strict_2))])
# Get nr specimens per species to calculate %
physeq_phylum_FOO_Specimen_sppFreq = as(sample_data(physeq_phylum_FOO_Specimen), "data.frame")
phylum_specimen_freq <- table(physeq_phylum_FOO_Specimen_sppFreq$SpeciesID)
# Calculate FOO%
phylum_Specimen_spp_FOOPerc <- t(apply(phylum_Specimen_spp, 2, "/", phylum_specimen_freq))

# 2.4.4 Make plot with FOO% of phyla over all species

# Melt FOO% data
bar_phylum_sp_FOO <- melt(phylum_Specimen_spp_FOOPerc, id.vars = "SpeciesID", variable.name = "Phylum")
# Add metadata
bar_phylum_sp_FOO$value <- as.numeric(bar_phylum_sp_FOO$value)
bar_phylum_sp_FOO$Phylum <- factor(bar_phylum_sp_FOO$Var1)
bar_phylum_sp_FOO$Tribe = with(all_metadata_2, Tribe[match(bar_phylum_sp_FOO$Var2,all_metadata_2$SpeciesID)])
# Make FOO% barplot
pdf(file = "Plot_FOO_phylum_19092024.pdf",
    width = 21,
    height = 10) 
ggplot(data=bar_phylum_sp_FOO,aes(x=Phylum,y=value,fill=Phylum) ) +
  geom_bar(position="dodge",stat = "identity") +
  geom_col(color = "black") +
  ylab("FOO%") +
  facet_grid(.~Var2, scales = "free", switch = "x", space = "free_x") +
  facet_wrap(factor(Tribe, level=Tribe_ordering)~factor(Var2, level=Species_ordering),scales = "free",nrow=5) +
  scale_fill_manual(values=phyla_cols) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) 
dev.off()

# 2.5 The Bray-Curtis distance
################################################################################

# 2.5.1 Calculate BC distances

# Take species physeq
otu_phylum_spp <- as.data.frame(otu_table(physeq_phylum_Specimen_spp_rel))
# Calculate Bray-Curtis distance
dist.BrCu.spp_phylum<-as.matrix(vegdist(otu_phylum_spp,method="bray"))
# Cluster species based on Ward
clust.BrCu.spp_phylum<-hclust(vegdist(otu_phylum_spp,method="bray"), method="ward.D") #agglomerative clustering using complete linkage

# 2.5.2 Construct heatmap
library(pheatmap)

# Define metadata for heatmap
df3 <- as.matrix(Metadata_spp)
df_phylum_spp = data.frame(
  Tribe = factor(df3[,c("Tribe")])
)
heatmap_colour = list(
  Tribe = tribe_colours)
# Plot heatmap
pdf(file = "Plot_BrCu_heatmap_19092024.pdf",
    width = 15,
    height = 13) 
pheatmap(1-dist.BrCu.spp_phylum,cluster_cols= clust.BrCu.spp_phylum,cluster_rows= clust.BrCu.spp_phylum,show_rownames=TRUE,show_colnames=T,annotation_col=df_phylum_spp,annotation_row=df_phylum_spp,annotation_colors=heatmap_colour)
dev.off()
