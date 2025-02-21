###########################################################
# 04. Order level initial analysis Focus dataset
###########################################################
# RScript to analyze the metagenomic diet data at order level 

# 4.1 Prepare data to order level
################################################################################

# 4.1.1 Subset phyloseq object

# Make physeq with identification to order level
physeq_order <- tax_glom(physeq, taxrank=rank_names(physeq)[4], NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
df_physeq_order <- as.data.frame(otu_table(physeq_order))

# 4.1.2 Filter data

df_physeq_order[(df_physeq_order) < 20 ] <- 0
rowSums(df_physeq_order)
df_physeq_order <- df_physeq_order[rowSums(df_physeq_order) >= 500,]
df_physeq_order[apply(df_physeq_order,2,function(x){x/sum(x)}) < 0.005 ] <- 0
# Add to physeq_order
otu_table(physeq_order) <- phyloseq::otu_table(df_physeq_order, taxa_are_rows = TRUE)
# Remove Samples that have fewer than 100 reads left
physeq_order = prune_samples(sample_sums(physeq_order)>=100, physeq_order)

# 4.1.3 Check data

physeq_order
physeq_order_samples = data.frame(sample_data(physeq_order))
summary(unique(physeq_order_samples$SpeciesID))
# Check nr reads per sample
summary(colSums(otu_table(physeq_order)))
sd(colSums(otu_table(physeq_order)))

# 4.2 Diet abundance
################################################################################

# 4.2.1 Combine species diet abundance

# Merge samples by Specimen
physeq_order_Specimen0 <- physeq_order
physeq_order_Specimen <- merge_samples(physeq_order_Specimen0,"SpecimenID",fun=sum)
# Take relative value
physeq_order_Specimen1 <- phyloseq::transform_sample_counts(physeq_order_Specimen, function(x) (x / sum(x)*100))
# Merge by species
sample_data(physeq_order_Specimen1)$SpeciesID <-with(all_metadata_2, SpeciesID[match(rownames(sample_data(physeq_order_Specimen1)),all_metadata_2$SpecimenID)])
physeq_order_Specimen_spp <- merge_samples(physeq_order_Specimen1, "SpeciesID", fun=mean) 
# Take relative values
physeq_order_Specimen_spp_rel <- phyloseq::transform_sample_counts(physeq_order_Specimen_spp, function(x) (x / sum(x)*100))

# 4.2.2 Define diet composition barplots

# Rename TaxID's
order_Specimen_spp_rel = as(otu_table(physeq_order_Specimen_spp_rel), "matrix")
colnames(order_Specimen_spp_rel) <- with(Tax_table_strict_2, order[match(colnames(order_Specimen_spp_rel),rownames(Tax_table_strict_2))])
unique(colnames(order_Specimen_spp_rel))
# Make table for abundance tables
bar_order_sp_otu <- melt(order_Specimen_spp_rel, id.vars = "SpeciesID", variable.name = "order")
bar_order_sp_otu$value <- as.numeric(bar_order_sp_otu$value)
bar_order_sp_otu$Var2 <- factor(bar_order_sp_otu$Var2)
bar_order_sp_otu$order <- factor(bar_order_sp_otu$Var2)
bar_order_sp_otu$Var1 <- factor(bar_order_sp_otu$Var1)
# Add metadata
Metadata_spp <- all_metadata_2[c("SpeciesID", "Tribe","BreedingType","d15N_mean", "d13C_mean", "FoodCat","Habitat")]
Metadata_spp <- Metadata_spp[!duplicated(Metadata_spp), ]
rownames(Metadata_spp) <- Metadata_spp$SpeciesID
Metadata_spp$d15N_mean <- as.numeric(Metadata_spp$d15N_mean)
Metadata_spp$d13C_mean <- as.numeric(Metadata_spp$d13C_mean)
# Define colours orders
order_cols = c("Cyprinodontiformes"="powderblue","Carangaria incertae sedis Latidae"="royalblue4","Synbranchiformes"="slategray3","Cypriniformes"="royalblue1","Siluriformes"="turquoise3", "Clupeiformes"="lightcyan", "Osteoglossiformes"="lightskyblue" ,"Polypteriformes"="mediumblue","Ceratodontiformes"="deepskyblue4","Characiformes"="turquoise","Decapoda"="sandybrown", "Podocopida"="coral4","Diplostraca"= "tan1" , "Diptera"="wheat","Trichoptera"="#82643C", "Ephemeroptera"="lightyellow","Hemiptera"="oldlace","Odonata"="tan","Ixodida"="tan4","Araneae"="tan3","Unionida"="violetred","Architaenioglossa"="violetred4","Tubificida"="mediumpurple3","Limnomedusae"="ghostwhite","Spongillida"="khaki1","Hygrophila"="thistle1", "Philodinida"="gray20", "Alismatales"="springgreen4" ,"Fagales"="springgreen","Solanales"="green3","Saxifragales"="mediumseagreen","Poales"="greenyellow","Arecales"="chartreuse2","Asparagales"="darkgreen", "Alismatales"="chartreuse3", "Ranunculales"="forestgreen" ,"Nymphaeales"="lightgreen","Zygnematales"="olivedrab4","Charales"="darkseagreen1", "Oedogoniales"="darkseagreen3","Chlorellales"="darkseagreen2","Cladophorales"="lightseagreen", "Naviculales"="olivedrab2", "Anaulales"="olivedrab1","Bonnemaisoniales"="olivedrab3","Eustigmatales"="darkolivegreen")
# Make stacked barplot
p + geom_tree(size=1.3,aes(colour=Tribe)) + 
  scale_colour_manual(values=tribe_colours) + 
  geom_treescale()+
  geom_tiplab(geom="text",offset=0,size=3,as_ylab=FALSE) +
  geom_fruit(data=bar_order_sp_otu, geom=geom_bar,mapping=aes(y=Var1, x=value,fill = order), stat="identity", width = 0.9,pwidth=2,offset=0.12)+
  scale_fill_manual(values=order_cols)

# 4.3 Diet abundance plots at phylum, class and order
################################################################################

# 4.3.1 Plot barplots

pdf(file = "Plot_abundances_focus_PhClOr.pdf",
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

