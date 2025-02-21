###########################################################
# 02. Phylum level analysis Focus dataset
###########################################################
# RScript to analyze the metagenomic diet data at phylum level 

# 2.1 Prepare data to phylum level
################################################################################

# 2.1.1 Subset phyloseq object

# Make physeq with identification to phylum level
physeq_phylum <- tax_glom(physeq, taxrank=rank_names(physeq)[2], NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
df_physeq_phylum <- as.data.frame(otu_table(physeq_phylum))

# 2.1.2 Filter data

df_physeq_phylum[(df_physeq_phylum) < 20 ] <- 0
df_physeq_phylum <- df_physeq_phylum[rowSums(df_physeq_phylum) >= 500,]
df_physeq_phylum[apply(df_physeq_phylum,2,function(x){x/sum(x)}) < 0.005 ] <- 0
# Add to physeq_phylum
otu_table(physeq_phylum) <- phyloseq::otu_table(df_physeq_phylum, taxa_are_rows = TRUE)
# Remove Samples that have less than 100 reads left
physeq_phylum = prune_samples(sample_sums(physeq_phylum)>=100, physeq_phylum)

# 2.2 Diet abundance
################################################################################

# 2.2.1 Combine species diet abundance

# Merge samples by Specimen
physeq_phylum_Specimen0 <- physeq_phylum
physeq_phylum_Specimen <- merge_samples(physeq_phylum_Specimen0,"SpecimenID",fun=sum)
# Take relative values
physeq_phylum_Specimen1 <- phyloseq::transform_sample_counts(physeq_phylum_Specimen, function(x) (x / sum(x)*100))
# Merge by SpeciesID
sample_data(physeq_phylum_Specimen1)$SpeciesID <-with(all_metadata_2, SpeciesID[match(rownames(sample_data(physeq_phylum_Specimen1)),all_metadata_2$SpecimenID)])
physeq_phylum_Specimen_spp <- merge_samples(physeq_phylum_Specimen1, "SpeciesID", fun=mean) 
# Take relative values
physeq_phylum_Specimen_spp_rel <- phyloseq::transform_sample_counts(physeq_phylum_Specimen_spp, function(x) (x / sum(x)*100))
# Filter phylogenetic tree by only species in this study
species <- unique(sample_data(physeq_phylum)$SpeciesID)
pruned.tree = drop.tip(pruned.tree, pruned.tree$tip.label[-which(pruned.tree$tip.label %in% species)])

# 2.2.2 Define diet composition barplots
library(ggplot2)
library(phytools)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(reshape2)

# Extract OTU table
phylum_Specimen_spp_rel = as(otu_table(physeq_phylum_Specimen_spp_rel), "matrix")
# Rename TaxID's
colnames(phylum_Specimen_spp_rel) <- with(Tax_table_strict_2, phylum[match(colnames(phylum_Specimen_spp_rel),rownames(Tax_table_strict_2))])
# Melt dataset for plot
bar_phylum_sp_otu <- melt(phylum_Specimen_spp_rel, id.vars = "SpeciesID", variable.name = "Phylum")
bar_phylum_sp_otu$value <- as.numeric(bar_phylum_sp_otu$value)
bar_phylum_sp_otu$Var2 <- factor(bar_phylum_sp_otu$Var2)
bar_phylum_sp_otu$Var1 <- factor(bar_phylum_sp_otu$Var1)
# Add metadata
Metadata_spp <- all_metadata_2[c("SpeciesID", "Tribe","BreedingType","d15N_mean", "d13C_mean", "FoodCat","Habitat")]
Metadata_spp <- Metadata_spp[!duplicated(Metadata_spp), ]
rownames(Metadata_spp) <- Metadata_spp$SpeciesID
Metadata_spp$d15N_mean <- as.numeric(Metadata_spp$d15N_mean)
Metadata_spp$d13C_mean <- as.numeric(Metadata_spp$d13C_mean)
# Make basis figure with phylogenetic tree
p <- ggtree(pruned.tree,ladderize = FALSE)
# Define colours of the tribes
tribe_colours <- c("Bathybatini"="#242626", "Benthochromini"="#AE262A", "Boulengerochromini"="#59595C", "Cyphotilapiini"="#FDDF13","Cyprichromini"="#F04D29","Ectodini"="#9AB9D9","Eretmodini"="#682E7A","Haplochromini"="darkgreen","Lamprologini"="#C588BB","Limnochromini"="#535CA9", "Perissodini"="orange","Trematocarini"="#959170", "Tropheini"="#86C773","Oreochromini"="grey")
# Assign tribes to tree
p$data$Tribe <- with(Metadata_spp, Tribe[match(p$data$label,rownames(Metadata_spp))])
p + geom_tree(size=1.5,aes(colour=Tribe)) + scale_colour_manual(values=tribe_colours) + geom_text(aes(label=node), hjust=-.3) + geom_treescale()
p$data$Tribe[p$data$parent >= 85 & p$data$parent < 94 ] <- "Tropheini"
p$data$Tribe[p$data$parent == 84 & p$data$node == 85 ] <- "Tropheini"
p$data$Tribe[p$data$parent == 83 & p$data$node == 95 ] <- "Eretmodini"
p$data$Tribe[p$data$parent == 96 & p$data$node == 98 ] <- "Cyprichromini"
p$data$Tribe[p$data$parent >= 101 & p$data$parent < 108 ] <- "Ectodini"
p$data$Tribe[p$data$parent == 100 & p$data$node == 101 ] <- "Ectodini"
p$data$Tribe[p$data$parent == 100 & p$data$node == 109 ] <- "Limnochromini"
p$data$Tribe[p$data$parent == 99 & p$data$node == 110 ] <- "Cyphotilapiini"
p$data$Tribe[p$data$parent >= 62 & p$data$parent < 80 ] <- "Lamprologini"
p$data$Tribe[p$data$parent == 61 & p$data$node == 62 ] <- "Lamprologini"
p$data$Tribe[p$data$parent >= 114 & p$data$node == 115 ] <- "Trematocarini"
p$data$Tribe[p$data$parent == 112 & p$data$node == 114 ] <- "Trematocarini"
p$data$Tribe[p$data$parent == 112 & p$data$node == 113 ] <- "Bathybatini"
# Define colors phyla
phyla_cols = c("Annelida"="lightpink3","Arthropoda"="peru", "Bacillariophyta"="darkseagreen2", "Chlorophyta"="palegreen3","Algae"= "#2f7d44" , "Rhodophyta" ="darkseagreen"  , "Chordata"="#224e66","Cnidaria"="ghostwhite","Mollusca"="brown","Nematoda"="indianred4","Platyhelminthes"="lightpink3","Porifera"="khaki1","Streptophyta"="mediumseagreen","Myzozoa"="gray38","Rotifera"="lightsalmon4")
# Make stacked barplot with phylogeny of cichlids
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
  scale_fill_distiller(direction=1,palette='YlOrBr') 

# 2.3 Normalize and transform data
################################################################################

# 2.3.1 Perform DESeq and vst
library(DESeq2)

# Preparation DESeq object
phylum_otu <- as.data.frame(t(otu_table(physeq_phylum_Specimen)+1))
phylum_sample <- as.data.frame(sample_data(physeq_phylum_Specimen))
phylum_tax <- as.data.frame(tax_table(physeq_phylum_Specimen))
# Factor
phylum_sample$SpecimenID<-factor(phylum_sample$SpecimenID)
phylum_sample$TissueID<-factor(phylum_sample$TissueID)
phylum_sample$Sex<-factor(phylum_sample$Sex)
phylum_sample$SpeciesID<-factor(phylum_sample$SpeciesID)
# coldata needs to be ordered same way as counts for deseq
phylum_sample <- phylum_sample[order(rownames(phylum_sample)),]
phylum_otu <- phylum_otu[,order(colnames(phylum_otu))]
# Add metadata
phylum_sample$SpeciesID = with(all_metadata_2, SpeciesID[match(rownames(phylum_sample),SpecimenID)])
phylum_sample$Tribe = with(all_metadata_2, Tribe[match(rownames(phylum_sample),SpecimenID)])
# Make a DeSeq object
dds <- DESeqDataSetFromMatrix(countData=phylum_otu, 
                              colData=phylum_sample, 
                              design=~SpeciesID) # Specifies how the counts from each OTU depend on our variables in the metadata, this case SpeciesID
# Run DeSeq on object 
dds2 <- DESeq(dds)
# Check results
res <- results(dds2)
summary (res)
# Then we need to transform the data
# vst function will perform variance stabilizing transformation
vsdata <- varianceStabilizingTransformation(dds2)

# 2.4 Perform principal component analysis
################################################################################

# 2.4 PCA

bf_pca_phylum <- data.frame(t(assay(vsdata)))
# Perform PCA
vsdata_res_phylum <- prcomp((bf_pca_phylum), scale=F)
summary(vsdata_res_phylum)
# Add metadata
bf_pca_phylum$Tribe = with(all_metadata_2, Tribe[match(rownames(bf_pca_phylum),all_metadata_2$SpecimenID)])
bf_pca_phylum$Sex = with(all_metadata_2, Sex[match(rownames(bf_pca_phylum),all_metadata_2$SpecimenID)])
bf_pca_phylum$SpeciesID = with(all_metadata_2, SpeciesID[match(rownames(bf_pca_phylum),all_metadata_2$SpecimenID)])
# Make PCA plot with TaxIDs
autoplot(vsdata_res_phylum,  data = bf_pca_phylum,colour='Tribe' ,loadings = TRUE, loadings.colour = 'gray24',loadings.label = T, loadings.label.size = 3,x = 1,y=2) +
  coord_fixed() +
  scale_color_manual(values=tribe_colours) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  ggtitle("PCA Diet Phylum lvl") +
  theme(plot.title = element_text(hjust = 0.5)) 
# Make PCA plot Supplementary Fig
autoplot(vsdata_res_phylum,  data = bf_pca_phylum,colour='Tribe' ,loadings = TRUE, loadings.colour = 'gray24',loadings.label = F,  size=4,scale=F) +
  coord_fixed() +
  scale_color_manual(values=tribe_colours) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  ggtitle("PCA Diet Phylum lvl") +
  theme(axis.text.x = element_text(color = "black", size = 11, face = "plain"),
        axis.text.y = element_text(color = "black", size = 11),  
        axis.title.x = element_text(color = "black", size = 14, vjust = -1, hjust=0.45),
        axis.title.y = element_text(color = "black", size = 14, vjust = 2),
        legend.title = element_text(size=15), legend.text = element_text(size=12),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + 
  theme(aspect.ratio = 1.1)


