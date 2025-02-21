###########################################################
# 06. Order level analysis
###########################################################
# RScript to analyze the metagenomic diet data at order level 

# 6.1 Prepare data to order level analyses
################################################################################

# 6.1.1 Subset phyloseq object

# Make physeq with identification to order level
physeq_order <- tax_glom(physeq, taxrank=rank_names(physeq)[4], NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
df_physeq_order <- as.data.frame(otu_table(physeq_order))

# 6.1.2 Filter data

# Remove hits with fewer than 20 counts
df_physeq_order[(df_physeq_order) < 20 ] <- 0
# Remove orders with fewer than 500 reads
df_physeq_order <- df_physeq_order[rowSums(df_physeq_order) >= 500,]
# Remove hits with fewer than 0.5 % relative abundance in a sample
df_physeq_order[apply(df_physeq_order,2,function(x){x/sum(x)}) < 0.005 ] <- 0
# Add to filtered data back to phyloseq object
otu_table(physeq_order) <- phyloseq::otu_table(df_physeq_order, taxa_are_rows = TRUE)
# Remove Samples that have less than 100 reads left
physeq_order = prune_samples(sample_sums(physeq_order)>=100, physeq_order)

# 6.2 Investigate data and calculate diversity
################################################################################

# 6.2.1 Check data

physeq_order
# Check nr cichlid species
length(unique(sample_data(physeq_order)$SpeciesID))
# Check min, max, mean, median and sd per sample
summary(colSums(otu_table(physeq_order)))
sd(colSums(otu_table(physeq_order)))

# 6.2.1 Calculate diversity

# Plot diversity
div_phl_order <- plot_richness(physeq_order, x='SpeciesID', measures=c("Shannon","Chao1"), color="Tribe") +
  geom_point(size=2) +
  geom_boxplot(aes(fill=Tribe),alpha=0.7,size=0.5 ) +
  theme_minimal() +
  theme(axis.title=element_text(size=10),
        legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=10),axis.text.y = element_text( size=8))+
  scale_colour_manual(values=tribe_colours) +
  scale_fill_manual(values=tribe_colours) + 
  scale_x_discrete(limits=rev(pruned.tree$tip.label))


# 6.3 Diet order composition
################################################################################

# 6.3.1 Merge data per species

# Merge samples by Specimen
physeq_order_Specimen <- merge_samples(physeq_order,"SpecimenID",fun=sum)
# Take relative of each sample for equal contribution when merging
physeq_order_Specimen1 <- phyloseq::transform_sample_counts(physeq_order_Specimen, function(x) (x / sum(x)*100))
# Then merge by species
sample_data(physeq_order_Specimen1)$SpeciesID <-with(all_metadata_2, SpeciesID[match(rownames(sample_data(physeq_order_Specimen1)),all_metadata_2$SpecimenID)])
physeq_order_Specimen_spp <- merge_samples(physeq_order_Specimen1, "SpeciesID", fun=mean) 
# Add metadata
sample_data(physeq_order_Specimen_spp)$d15N_mean <-with(all_metadata_2, d15N_mean[match(rownames(sample_data(physeq_order_Specimen_spp)),all_metadata_2$SpeciesID)])
sample_data(physeq_order_Specimen_spp)$d15N_std_dev <-with(all_metadata_2, d15N_std_dev[match(rownames(sample_data(physeq_order_Specimen_spp)),all_metadata_2$SpeciesID)])
sample_data(physeq_order_Specimen_spp)$Sampling <-with(all_metadata_2, Sampling[match(rownames(sample_data(physeq_order_Specimen_spp)),all_metadata_2$SpeciesID)])
sample_data(physeq_order_Specimen_spp)$Tribe <-with(all_metadata_2, Tribe[match(rownames(sample_data(physeq_order_Specimen_spp)),all_metadata_2$SpeciesID)])
sample_data(physeq_order_Specimen_spp)$BreedingType <-with(all_metadata_2, BreedingType[match(rownames(sample_data(physeq_order_Specimen_spp)),all_metadata_2$SpeciesID)])
sample_data(physeq_order_Specimen_spp)$BreedingMode <-with(all_metadata_2, BreedingMode[match(rownames(sample_data(physeq_order_Specimen_spp)),all_metadata_2$SpeciesID)])
sample_data(physeq_order_Specimen_spp)$d13C_mean <-with(all_metadata_2, d13C_mean[match(rownames(sample_data(physeq_order_Specimen_spp)),all_metadata_2$SpeciesID)])
sample_data(physeq_order_Specimen_spp)$d13C_std_dev <-with(all_metadata_2, d13C_std_dev[match(rownames(sample_data(physeq_order_Specimen_spp)),all_metadata_2$SpeciesID)])
sample_data(physeq_order_Specimen_spp)$FoodCat <-with(all_metadata_2, FoodCat[match(rownames(sample_data(physeq_order_Specimen_spp)),all_metadata_2$SpeciesID)])
sample_data(physeq_order_Specimen_spp)$Habitat <-with(all_metadata_2, Habitat[match(rownames(sample_data(physeq_order_Specimen_spp)),all_metadata_2$SpeciesID)])
# Take relative value
physeq_order_Specimen_spp_rel <- phyloseq::transform_sample_counts(physeq_order_Specimen_spp, function(x) (x / sum(x)*100))

# 6.3.2 Calculate abundance and occurrence orders 

# Extract otu table
order_Specimen_spp_rel = as(otu_table(physeq_order_Specimen_spp_rel), "matrix")
# Rename TaxID's
colnames(order_Specimen_spp_rel) <- with(Tax_table_strict_2, order[match(colnames(order_Specimen_spp_rel),rownames(Tax_table_strict_2))])
# Check total abundance per species
colSums(order_Specimen_spp_rel)/nrow(order_Specimen_spp_rel)
# Calculate FOO
order_Specimen_spp_rel2 <- order_Specimen_spp_rel
order_Specimen_spp_rel2[order_Specimen_spp_rel2 >= 1] <- 1
order_Specimen_spp_rel2[order_Specimen_spp_rel2 < 1] <- 0
# Check FOO of orders
colSums(order_Specimen_spp_rel2)
# Check FOO% of orders
colSums(order_Specimen_spp_rel2)/nrow(order_Specimen_spp_rel2)

# 6.3.3 Plot abundances

# Define colours for orders
order_cols = c("Cyprinodontiformes"="powderblue","Carangaria incertae sedis Latidae"="royalblue4","Synbranchiformes"="slategray3","Cypriniformes"="royalblue1","Siluriformes"="turquoise3", "Clupeiformes"="lightcyan", "Osteoglossiformes"="lightskyblue" ,"Polypteriformes"="mediumblue","Decapoda"="sandybrown","Amphipoda"="peachpuff1", "Calanoida"="coral","Cyclopoida"="coral3","Podocopida"="coral4","Diplostraca"= "tan1" , "Diptera"="wheat", "Lepidoptera"="lightyellow","Hymenoptera"="oldlace","Odonata"="tan","Ixodida"="tan4","Araneae"="tan3","Architaenioglossa"="violetred4","Rhynchobdellida"="mediumpurple3","Limnomedusae"="ghostwhite","Spongillida"="khaki1","Fabales"="springgreen3", "Rosales"="springgreen2", "Malpighiales"="springgreen4" ,"Fagales"="springgreen","Solanales"="green3","Saxifragales"="mediumseagreen","Poales"="greenyellow","Arecales"="chartreuse2","Asparagales"="darkgreen", "Alismatales"="chartreuse3", "Ceratophyllales"="forestgreen" ,"Nymphaeales"="lightgreen","Sphaeropleales"="darkseagreen4","Chlamydomonadales"="darkseagreen1", "Oedogoniales"="darkseagreen3","Chlorellales"="darkseagreen2","Cladophorales"="lightseagreen","Bacillariales"= "olivedrab4" , "Naviculales"="olivedrab2", "Anaulales"="olivedrab1","Thalassiosirales"="olivedrab3","Eustigmatales"="darkolivegreen")
# Plot stacked barplots
plot_bar(physeq_order_Specimen_spp_rel, fill="order") + 
  geom_bar(aes(fill = order), stat="identity", position="stack") +
  labs(x = "", y = "Relative Abundance\n") +
  theme(strip.text.x = element_text(size = 12, color = "black" ),
        strip.background = element_rect(color="grey", fill="white", size=1.5, linetype="solid")) +
  scale_x_discrete()  + 
  theme(panel.background = element_blank()) +
  scale_fill_manual(values=order_cols) +
  scale_x_discrete(limits=rev(Species_ordering)) + 
  coord_flip()
# Plot with tree
# Extract otu_table
order_Specimen_spp_rel = as(otu_table(physeq_order_Specimen_spp_rel), "matrix")
# Rename TaxID's
colnames(order_Specimen_spp_rel) <- with(Tax_table_strict_2, order[match(colnames(order_Specimen_spp_rel),rownames(Tax_table_strict_2))])
# Make table for abundance tables
bar_order_sp_otu <- melt(order_Specimen_spp_rel, id.vars = "SpeciesID", variable.name = "order")
bar_order_sp_otu$value <- as.numeric(bar_order_sp_otu$value)
bar_order_sp_otu$Var2 <- factor(bar_order_sp_otu$Var2)
bar_order_sp_otu$order <- factor(bar_order_sp_otu$Var2)
bar_order_sp_otu$Var1 <- factor(bar_order_sp_otu$Var1)


