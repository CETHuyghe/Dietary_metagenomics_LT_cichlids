###########################################################
# 10. All ranks analysis
###########################################################
# RScript to analyze the metagenomic diet data over multiple taxonomic ranks 

# 10.1 Gut length plot
################################################################################

# 10.1.1 Plot the Zihler's Index and relative gut length
library(ggplot2)
library(ggimage)

# Define tribe colours
tribe_colours <- c("Bathybatini"="#242626", "Benthochromini"="#AE262A", "Boulengerochromini"="#59595C", "Cyphotilapiini"="#FDDF13","Cyprichromini"="#F04D29","Ectodini"="#9AB9D9","Eretmodini"="#682E7A","Haplochromini"="darkgreen","Lamprologini"="#C588BB","Limnochromini"="#535CA9", "Perissodini"="orange","Trematocarini"="#959170", "Tropheini"="#86C773","Oreochromini"="grey")
# Define Zihler's Index boxplot
zhl_plot <- ggplot(LNG_CETH, aes(SpeciesID,Zihler, fill=Tribe)) +
  geom_boxplot() +
  theme_minimal() +
  theme(axis.title=element_text(size=12),
        axis.title.y = element_text(vjust = 2.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=10, color = "black"),
        axis.text.y = element_text( size=10, color = "black"),
        legend.position="none")+
  scale_fill_manual(values=tribe_colours) + 
  scale_x_discrete(limits=rev(pruned.tree$tip.label)) +
  xlab("Chichlid species") +
  ylab("Zihler's Index for gut length")
# Define relative gut length boxplot
lng_plot <- ggplot(LNG_CETH, aes(SpeciesID,RGL, fill=Tribe)) +
  geom_boxplot() +
  theme_minimal() +
  theme(axis.title=element_text(size=12),
        axis.title.y = element_text(vjust = 2.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=10, color = "black"),
        axis.text.y = element_text( size=10, color = "black"),
        legend.position="none")+
  scale_fill_manual(values=tribe_colours) + 
  scale_x_discrete(limits=rev(pruned.tree$tip.label)) +
  xlab("Chichlid species") +
  ylab("Relative gut length")

# Combine boxplots
library(ggimage)
pdf(file = "Plot_length_RGL_ZHL.pdf",   # The directory you want to save the file in
    width = 22, # The width of the plot in inches
    height = 10) # The height of the plot in inches
grid.arrange(zhl_plot,lng_plot)
dev.off()


# 10.2 Number of reads per rank
################################################################################

# 10.2.1 Calculate total nr reads per sample over all ranks

# Make dataset for total nr reads
physeq_phylum_df <- as.data.frame(t(otu_table(physeq_phylum)))
physeq_phylum_df$rowsum <- rowSums(physeq_phylum_df)
all_metadata_2$rowsum_phylum <-with(physeq_phylum_df, rowsum[match(rownames(all_metadata_2),rownames(physeq_phylum_df))])
physeq_class_df <- as.data.frame(t(otu_table(physeq_class)))
physeq_class_df$rowsum <- rowSums(physeq_class_df)
all_metadata_2$rowsum_class <-with(physeq_class_df, rowsum[match(rownames(all_metadata_2),rownames(physeq_class_df))])
physeq_order_df <- as.data.frame(t(otu_table(physeq_order)))
physeq_order_df$rowsum <- rowSums(physeq_order_df)
all_metadata_2$rowsum_order <-with(physeq_order_df, rowsum[match(rownames(all_metadata_2),rownames(physeq_order_df))])
physeq_family_df <- as.data.frame(t(otu_table(physeq_family)))
physeq_family_df$rowsum <- rowSums(physeq_family_df)
all_metadata_2$rowsum_family <-with(physeq_family_df, rowsum[match(rownames(all_metadata_2),rownames(physeq_family_df))])
physeq_genus_df <- as.data.frame(t(otu_table(physeq_genus)))
physeq_genus_df$rowsum <- rowSums(physeq_genus_df)
all_metadata_2$rowsum_genus <-with(physeq_genus_df, rowsum[match(rownames(all_metadata_2),rownames(physeq_genus_df))])
physeq_species_df <- as.data.frame(t(otu_table(physeq_species)))
physeq_species_df$rowsum <- rowSums(physeq_species_df)
all_metadata_2$rowsum_species <-with(physeq_species_df, rowsum[match(rownames(all_metadata_2),rownames(physeq_species_df))])
# Make seperate dataset with all total number of reads per sample
all_rowsums <- all_metadata_2[ c("rowsum_phylum","rowsum_class","rowsum_order","rowsum_family","rowsum_genus","rowsum_species") ]
all_rowsums$SampleID <- rownames(all_rowsums)

# 10.2.2 Barplot nr reads

# Melt for abundance plot
Reads_melt <- melt(all_rowsums, id.vars = "SampleID", variable.name = "Rank")
# Add metadata
Reads_melt$Tribe <-with(all_metadata_2, Tribe[match(Reads_melt$SampleID,rownames(all_metadata_2))])
Reads_melt$SpeciesID <-with(all_metadata_2, SpeciesID[match(Reads_melt$SampleID,rownames(all_metadata_2))])
Reads_melt$TissueTubeID <-with(all_metadata_2, TissueTubeID[match(Reads_melt$SampleID,rownames(all_metadata_2))])
Reads_melt$spot <- rownames(Reads_melt)
# Define factors
Reads_melt$SampleID <- as.factor(Reads_melt$SampleID)
Reads_melt$Rank <- as.factor(Reads_melt$Rank)
Reads_melt$Tribe <- as.factor(Reads_melt$Tribe)
Reads_melt$SpeciesID <- as.factor(Reads_melt$SpeciesID)
Reads_melt$spot <- as.factor(Reads_melt$spot)
# Combine SampleID and taxonomic Rank
varread1 = as.character(get_variable(Reads_melt, "SampleID"))
varread2 = as.character(get_variable(Reads_melt, "Rank"))
Reads_melt$SampleID_Rank <- mapply(paste0, varread1, varread2,collapse = "_")
# Order by species ID, then rank
reads_ordering <- c("rowsum_phylum","rowsum_class","rowsum_order","rowsum_family","rowsum_genus","rowsum_species")
Reads_melt$SampleID_Rank <- factor(Reads_melt$SampleID_Rank, levels=reads_ordering)
Reads_melt <- Reads_melt[order(Reads_melt$SampleID_Rank),]
# Rename column
Reads_melt$`Taxonomic rank` <- Reads_melt$Rank
# Rename tax ranks
Reads_melt <- Reads_melt %>%
  mutate(`Taxonomic rank` = recode(`Taxonomic rank`, rowsum_phylum = 'Phylum', rowsum_class = 'Class', rowsum_order = 'Order', rowsum_family = 'Family', rowsum_genus = 'Genus', rowsum_species = 'Species' ))
# Define colours for plot
reads_col <- c("lightskyblue","darkgreen","darkgoldenrod1","chocolate","firebrick","gray27")
# Make barplot
pdf(file = "Plot_NrReads_20092024.pdf",
    width = 21,
    height = 10) 
ggplot(data=Reads_melt, aes(x=SampleID, y=value, group=`Taxonomic rank`)) +
  xlab("Sample") + ylab("Number of total reads") +
  scale_fill_manual(values=reads_col) + 
  geom_bar(position="dodge",aes(fill=`Taxonomic rank`),stat = "identity") +
  facet_grid(.~SampleID, scales = "free", switch = "x", space = "free_x") +
  facet_wrap(Tribe~SpeciesID, scales=("free"),nrow=6) +
  theme_classic() +
  theme(axis.text.x = element_blank()) 
dev.off()

# 10.3 Abundance plot phylum, class and order
################################################################################

# 10.3.1 Plot diet composition per species

pdf(file = "Plot_abundances_PhClOr.pdf",
    width = 21,
    height = 10) 
p + geom_tree(size=1.3,aes(colour=Tribe)) + 
  scale_colour_manual(values=tribe_colours) + 
  geom_treescale()+
  geom_tiplab(geom="text",offset=0,size=3,as_ylab=FALSE) +
  geom_fruit(data=bar_phylum_sp_otu, geom=geom_bar,
             mapping=aes(y=Var1, x=value,fill = Var2), 
             stat="identity", width = 0.9,pwidth=2,offset=0.5)+
  scale_fill_manual(values=phyla_cols) +
  new_scale_fill() +
  geom_fruit(data=bar_class_sp_otu, 
             geom=geom_bar,mapping=aes(y=Var1, x=value,fill = Var2), 
             stat="identity", width = 0.9,pwidth=2,offset=0.12)+
  scale_fill_manual(values=class_cols) +
  new_scale_fill() +
  geom_fruit(data=bar_order_sp_otu, geom=geom_bar,
             mapping=aes(y=Var1, x=value,fill = order), 
             stat="identity", width = 0.9,pwidth=2,offset=0.12)+
  scale_fill_manual(values=order_cols) 
dev.off()

# 10.4 Diversity plot
################################################################################

# 10.4.1 Combine diversity plots phylum, class and order

pdf(file = "Plot_diversity.pdf",
    width = 22,
    height = 12) 
grid.arrange(div_phl_plot, div_phl_class, div_phl_order, ncol =1)
dev.off()