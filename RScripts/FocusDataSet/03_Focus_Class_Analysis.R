###########################################################
# 03. Class level initial analysis Focus dataset
###########################################################
# RScript to analyze the metagenomic diet data at class level 

# 3.1 Prepare data to class level
################################################################################

# 3.1.1 Subset phyloseq object

# Make physeq with identification to class level
physeq_class <- tax_glom(physeq, taxrank=rank_names(physeq)[3], NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
df_physeq_class <- as.data.frame(otu_table(physeq_class))

# 3.1.2 Filter data

df_physeq_class[(df_physeq_class) < 20 ] <- 0
df_physeq_class <- df_physeq_class[rowSums(df_physeq_class) >= 500,]
df_physeq_class[apply(df_physeq_class,2,function(x){x/sum(x)}) < 0.005 ] <- 0
# Add to physeq_class
otu_table(physeq_class) <- phyloseq::otu_table(df_physeq_class, taxa_are_rows = TRUE)
# Remove Samples that have fewer than 100 reads
physeq_class = prune_samples(sample_sums(physeq_class)>=100, physeq_class)

# 3.1.3 Check data

physeq_class
# Check nr reads per sample
summary(colSums(otu_table(physeq_class)))
sd(colSums(otu_table(physeq_class)))
# Check nr cichlid species
unique(sample_data(physeq_class)$SpeciesID)


# 3.2 Diet abundance
################################################################################

# 3.2.1 Combine species diet abundance

# Merge samples by Specimen
physeq_class_Specimen0 <- physeq_class
physeq_class_Specimen <- merge_samples(physeq_class_Specimen0,"SpecimenID",fun=sum)
# Take relative values
physeq_class_Specimen1 <- phyloseq::transform_sample_counts(physeq_class_Specimen, function(x) (x / sum(x)*100))
# Merge by species
sample_data(physeq_class_Specimen1)$SpeciesID <-with(all_metadata_2, SpeciesID[match(rownames(sample_data(physeq_class_Specimen1)),all_metadata_2$SpecimenID)])
physeq_class_Specimen_spp <- merge_samples(physeq_class_Specimen1, "SpeciesID", fun=mean) 

# 3.2.2 Define diet composition barplots

# Take rel
physeq_class_Specimen_spp_rel <- phyloseq::transform_sample_counts(physeq_class_Specimen_spp, function(x) (x / sum(x)*100))
# Add tree
class_Specimen_spp_rel = as(otu_table(physeq_class_Specimen_spp_rel), "matrix")
# Rename TaxID's
colnames(class_Specimen_spp_rel) <- with(Tax_table_strict_2, class[match(colnames(class_Specimen_spp_rel),rownames(Tax_table_strict_2))])
# Make table for abundance tables
bar_class_sp_otu <- melt(class_Specimen_spp_rel, id.vars = "SpeciesID", variable.name = "class")
bar_class_sp_otu$value <- as.numeric(bar_class_sp_otu$value)
bar_class_sp_otu$Var2 <- factor(bar_class_sp_otu$Var2)
bar_class_sp_otu$Var1 <- factor(bar_class_sp_otu$Var1)
# Add metadata
Metadata_spp <- all_metadata_2[c("SpeciesID", "Tribe","BreedingType","d15N_mean", "d13C_mean", "FoodCat","Habitat")]
Metadata_spp <- Metadata_spp[!duplicated(Metadata_spp), ]
rownames(Metadata_spp) <- Metadata_spp$SpeciesID
Metadata_spp$d15N_mean <- as.numeric(Metadata_spp$d15N_mean)
Metadata_spp$d13C_mean <- as.numeric(Metadata_spp$d13C_mean)
# Define class colours
class_cols = c("Clitellata"="mediumpurple3","Malacostraca"="sandybrown","Hexanauplia"="coral","Ostracoda"="coral4","Insecta"="lightcoral", "Arachnida"="tan3", "Branchiopoda"="peachpuff1" ,"Eustigmatophyceae"="olivedrab1","Bacillariophyceae"="olivedrab","Mediophyceae"="olivedrab2", "Coscinodiscophyceae"="olivedrab3","Chlorophyceae"="seagreen1","Trebouxiophyceae"="seagreen3","Ulvophyceae"= "palegreen" , "Florideophyceae"="lightseagreen","","Dipneusti"="lightskyblue","Actinopteri"="#224e66", "Cladistia"="royalblue1","Hydrozoa"="ghostwhite","Gastropoda"="brown","Bivalvia"="deeppink3","Eurotatoria"="mediumpurple4","Demospongiae"="khaki1","Charophyceae"="seagreen3","Zygnemophyceae"="seagreen1","Magnoliopsida"="mediumseagreen")
# Make stacked barplot per species
p + geom_tree(size=1.3,aes(colour=Tribe)) + 
  scale_colour_manual(values=tribe_colours) + 
  geom_treescale()+
  geom_tiplab(geom="text",offset=0,size=3,as_ylab=FALSE) +
  geom_fruit(data=bar_class_sp_otu, geom=geom_bar,mapping=aes(y=Var1, x=value,fill = Var2), stat="identity", width = 0.9,pwidth=2,offset=0.12)+
  scale_fill_manual(values=class_cols)
