###########################################################
# 05. Class level initial analysis
###########################################################
# RScript to analyze the metagenomic diet data at class level 

# 5.1 Prepare data to class level analyses
################################################################################

# 5.1.1 Subset phyloseq object

# Make physeq with identification to class level
physeq_class <- tax_glom(physeq, taxrank=rank_names(physeq)[3], NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
# Extract otu table to filter
df_physeq_class <- as.data.frame(otu_table(physeq_class))

# 5.1.2 Filter data

# Remove hits with fewer than 20 reads
df_physeq_class[(df_physeq_class) < 20 ] <- 0
# Check number of reads per class
rowSums(df_physeq_class)
# Remove classes with fewer than 500 reads
df_physeq_class <- df_physeq_class[rowSums(df_physeq_class) >= 500,]
# Remove hits with fewer than 0.05 % abundance in a sample
df_physeq_class[apply(df_physeq_class,2,function(x){x/sum(x)}) < 0.005 ] <- 0
# Add to physeq_class
otu_table(physeq_class) <- phyloseq::otu_table(df_physeq_class, taxa_are_rows = TRUE)
# Remove Samples that have fewer than 100 reads left
physeq_class = prune_samples(sample_sums(physeq_class)>=100, physeq_class)

# 5.2 Investigate data and calculate diversity
################################################################################

# 5.2.1 Check data

# Check nr of samples and classes
physeq_class
# Check the min, max, mean, median and sd of reads per sample
summary(colSums(otu_table(physeq_class)))
sd(colSums(otu_table(physeq_class)))
# Check nr cichlid species
length(unique(sample_data(physeq_class)$SpeciesID))

# 5.2.2 Calculate diversity

# Plot Shannon and Chao1 diversity
div_phl_class <- plot_richness(physeq_class, x='SpeciesID', measures=c("Shannon","Chao1"), color="Tribe") +
  geom_point(size=2) +
  geom_boxplot(aes(fill=Tribe),alpha=0.7,size=0.5 ) +
  theme_minimal() +
  theme(axis.title=element_text(size=10),
        legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=10),axis.text.y = element_text( size=8))+
  scale_colour_manual(values=tribe_colours) +
  scale_fill_manual(values=tribe_colours) + 
  scale_x_discrete(limits=rev(pruned.tree$tip.label))

# 5.3 Diet class composition
################################################################################

# 5.3.1 Merge data per species

# Sum samples by Specimen
physeq_class_Specimen0 <- physeq_class
physeq_class_Specimen <- merge_samples(physeq_class_Specimen0,"SpecimenID",fun=sum)
# Then calculate relative value, for equal contribution each specimen
physeq_class_Specimen1 <- phyloseq::transform_sample_counts(physeq_class_Specimen, function(x) (x / sum(x)*100))
# Add and merge by species
sample_data(physeq_class_Specimen1)$SpeciesID <-with(all_metadata_2, SpeciesID[match(rownames(sample_data(physeq_class_Specimen1)),all_metadata_2$SpecimenID)])
physeq_class_Specimen_spp <- merge_samples(physeq_class_Specimen1, "SpeciesID", fun=mean) 
# Add metadata
sample_data(physeq_class_Specimen_spp)$d15N_mean <-with(all_metadata_2, d15N_mean[match(rownames(sample_data(physeq_class_Specimen_spp)),all_metadata_2$SpeciesID)])
sample_data(physeq_class_Specimen_spp)$d15N_std_dev <-with(all_metadata_2, d15N_std_dev[match(rownames(sample_data(physeq_class_Specimen_spp)),all_metadata_2$SpeciesID)])
sample_data(physeq_class_Specimen_spp)$Sampling <-with(all_metadata_2, Sampling[match(rownames(sample_data(physeq_class_Specimen_spp)),all_metadata_2$SpeciesID)])
sample_data(physeq_class_Specimen_spp)$Tribe <-with(all_metadata_2, Tribe[match(rownames(sample_data(physeq_class_Specimen_spp)),all_metadata_2$SpeciesID)])
sample_data(physeq_class_Specimen_spp)$BreedingType <-with(all_metadata_2, BreedingType[match(rownames(sample_data(physeq_class_Specimen_spp)),all_metadata_2$SpeciesID)])
sample_data(physeq_class_Specimen_spp)$BreedingMode <-with(all_metadata_2, BreedingMode[match(rownames(sample_data(physeq_class_Specimen_spp)),all_metadata_2$SpeciesID)])
sample_data(physeq_class_Specimen_spp)$d13C_mean <-with(all_metadata_2, d13C_mean[match(rownames(sample_data(physeq_class_Specimen_spp)),all_metadata_2$SpeciesID)])
sample_data(physeq_class_Specimen_spp)$d13C_std_dev <-with(all_metadata_2, d13C_std_dev[match(rownames(sample_data(physeq_class_Specimen_spp)),all_metadata_2$SpeciesID)])
sample_data(physeq_class_Specimen_spp)$FoodCat <-with(all_metadata_2, FoodCat[match(rownames(sample_data(physeq_class_Specimen_spp)),all_metadata_2$SpeciesID)])
sample_data(physeq_class_Specimen_spp)$Habitat <-with(all_metadata_2, Habitat[match(rownames(sample_data(physeq_class_Specimen_spp)),all_metadata_2$SpeciesID)])
# Take relative values
physeq_class_Specimen_spp_rel <- phyloseq::transform_sample_counts(physeq_class_Specimen_spp, function(x) (x / sum(x)*100))

# 5.3.2 Calculate abundance and occurrence classes 

class_Specimen_spp_rel = as(otu_table(physeq_class_Specimen_spp_rel), "matrix")
# Rename TaxID's
colnames(class_Specimen_spp_rel) <- with(Tax_table_strict_2, class[match(colnames(class_Specimen_spp_rel),rownames(Tax_table_strict_2))])
# Calculate most abundant taxa per species
main_class <- as.data.frame(colnames(class_Specimen_spp_rel)[apply(class_Specimen_spp_rel,1,which.max)])
rownames(main_class) <- rownames(class_Specimen_spp_rel)
main_class$class <- main_class$`colnames(class_Specimen_spp_rel)[apply(class_Specimen_spp_rel, 1, which.max)]`
unique(main_class$class)
# Check these classes
length(which(main_class$class == "Malacostraca"))
length(which(main_class$class == "Actinopteri"))
length(which(main_class$class == "Ostracoda"))
length(which(main_class$class == "Magnoliopsida"))
length(which(main_class$class == "Bacillariophyceae"))
length(which(main_class$class == "Hexanauplia"))
length(which(main_class$class == "Insecta"))
length(which(main_class$class == "Ulvophyceae"))
length(which(main_class$class == "Hydrozoa"))
# Check overall abundance per Class
colSums(class_Specimen_spp_rel)/nrow(class_Specimen_spp_rel)
# Calculate overal occurence in a species
class_Specimen_spp_rel2 <- class_Specimen_spp_rel
class_Specimen_spp_rel2[class_Specimen_spp_rel2 >= 0.5] <- 1
class_Specimen_spp_rel2[class_Specimen_spp_rel2 < 0.5] <- 0
# Check in how many species the classes occur
colSums(class_Specimen_spp_rel2)
# Then relative over all species
colSums(class_Specimen_spp_rel2)/nrow(class_Specimen_spp_rel2)

# 5.3.3 Plot FOO%

#Take otu table per specimen
class_Specimen = as(otu_table(physeq_class_Specimen), "matrix")
# Make copy where we keep original TaxIDs
class_Specimen0 <- as.data.frame(t(class_Specimen))
# Rename TaxID's
colnames(class_Specimen) <- with(Tax_table_strict_2, class[match(colnames(class_Specimen),rownames(Tax_table_strict_2))])
# Look at FOO in species by specimen
class_Specimen2 <- as.data.frame(t(class_Specimen))
class_Specimen2[apply(class_Specimen2,2,function(x){x/sum(x)}) >= 0.005] <- 1
class_Specimen2[class_Specimen2 > 1] <- 0
rowSums(class_Specimen2)/ncol(class_Specimen2)
# Take FOO per species to put back into PhyloSeq object
class_Specimen0[apply(class_Specimen0,2,function(x){x/sum(x)}) >= 0.005] <- 1
class_Specimen0[class_Specimen0 > 1] <- 0
rowSums(class_Specimen0)
# Make Phyloseq object with FOO
physeq_class_FOO_Specimen <- physeq_class_Specimen
otu_table(physeq_class_FOO_Specimen) <- phyloseq::otu_table(class_Specimen0, taxa_are_rows = TRUE)
sample_data(physeq_class_FOO_Specimen)$SpeciesID <- with(all_metadata_2, SpeciesID[match(rownames(sample_data(physeq_class_FOO_Specimen)),all_metadata_2$SpecimenID)])
# Merge specimens per species
physeq_class_FOO_Specimen_spp <- merge_samples(physeq_class_FOO_Specimen, "SpeciesID", fun=sum) 
# Extract otu_table
class_Specimen_spp = as(otu_table(physeq_class_FOO_Specimen_spp), "matrix")
# rename TaxIDs
colnames(class_Specimen_spp) <- with(Tax_table_strict_2, class[match(colnames(class_Specimen_spp),rownames(Tax_table_strict_2))])
# Get nr specimens per species to calculate %
physeq_class_FOO_Specimen_sppFreq = as(sample_data(physeq_class_FOO_Specimen), "data.frame")
class_specimen_freq <- table(physeq_class_FOO_Specimen_sppFreq$SpeciesID)
# Calculate FOO%
class_Specimen_spp_FOOPerc <- t(apply(class_Specimen_spp, 2, "/", class_specimen_freq))
# Melt and prepare data for barplot
bar_class_sp_FOO <- melt(class_Specimen_spp_FOOPerc, id.vars = "SpeciesID", variable.name = "class")
bar_class_sp_FOO$value <- as.numeric(bar_class_sp_FOO$value)
bar_class_sp_FOO$Var2 <- factor(bar_class_sp_FOO$Var2)
bar_class_sp_FOO$class <- factor(bar_class_sp_FOO$Var1)
bar_class_sp_FOO$Tribe = with(all_metadata_2, Tribe[match(bar_class_sp_FOO$Var2,all_metadata_2$SpeciesID)])
bar_class_sp_FOO$Tribe <- factor(bar_class_sp_FOO$Tribe)
# Define colours for the classes
class_cols = c("Clitellata"="mediumpurple3","Malacostraca"="sandybrown","Hexanauplia"="coral","Ostracoda"="coral4","Insecta"="lightcoral", "Arachnida"="tan3", "Branchiopoda"="peachpuff1" ,"Eustigmatophyceae"="olivedrab1","Bacillariophyceae"="olivedrab","Mediophyceae"="olivedrab2", "Coscinodiscophyceae"="olivedrab3","Chlorophyceae"="seagreen1","Trebouxiophyceae"="seagreen3","Ulvophyceae"= "palegreen" , "Actinopteri"="#224e66", "Cladistia"="royalblue1","Hydrozoa"="ghostwhite","Gastropoda"="brown","Demospongiae"="khaki1","Magnoliopsida"="mediumseagreen")
# Make a barplot for FOO% per species
pdf(file = "Plot_FOO_Class_19092024.pdf",
    width = 21,
    height = 10) 
ggplot(data=bar_class_sp_FOO,aes(x=class,y=value,fill=class) ) +
  geom_bar(position="dodge",stat = "identity")  +
  ylab("FOO%") +
  facet_grid(.~Var2, scales = "free", switch = "x", space = "free_x") +
  facet_wrap(factor(Tribe, level=Tribe_ordering)~factor(Var2, level=Species_ordering),scales = "free",nrow=7) +
  scale_fill_manual(values=class_cols) +
  geom_col_pattern(
    aes(pattern=class, pattern_angle=class, pattern_spacing=class, fill=class,pattern_colour=class), 
    colour          = 'black', 
    pattern_density = 0.5
  ) +
  scale_pattern_colour_manual(values=class_cols) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) 
dev.off()

# 5.3.4 Plot abundances

# Plot stacked barplot per species 
plot_bar(physeq_class_Specimen_spp_rel, fill="class") + 
  geom_bar(aes(fill = class), stat="identity", position="stack") +
  labs(x = "", y = "Relative Abundance\n") +
  theme(strip.text.x = element_text(size = 12, color = "black" ),
        strip.background = element_rect(color="grey", fill="white", size=1.5, linetype="solid")) +
  scale_x_discrete()  + 
  theme(panel.background = element_blank()) +
  scale_fill_manual(values=class_cols) +
  scale_x_discrete(limits=rev(Species_ordering)) + 
  coord_flip()
# Make plot with tree and ordering as before by taxa
class_Specimen_spp_rel = as(otu_table(physeq_class_Specimen_spp_rel), "matrix")
# Rename TaxID's
colnames(class_Specimen_spp_rel) <- with(Tax_table_strict_2, class[match(colnames(class_Specimen_spp_rel),rownames(Tax_table_strict_2))])
# Make table for abundance tables
bar_class_sp_otu <- melt(class_Specimen_spp_rel, id.vars = "SpeciesID", variable.name = "class")
bar_class_sp_otu$value <- as.numeric(bar_class_sp_otu$value)
bar_class_sp_otu$Var2 <- factor(bar_class_sp_otu$Var2)
bar_class_sp_otu$Var1 <- factor(bar_class_sp_otu$Var1)
# Plot stacked barplot with tree
p + geom_tree(size=1.3,aes(colour=Tribe)) + 
  scale_colour_manual(values=tribe_colours) + 
  geom_treescale()+
  geom_tiplab(geom="text",offset=0,size=3,as_ylab=FALSE) +
  geom_fruit(data=bar_class_sp_otu, 
             geom=geom_bar,mapping=aes(y=Var1, x=value,fill = Var2), 
             stat="identity", width = 0.9,pwidth=2,offset=0.12)+
  scale_fill_manual(values=class_cols)



