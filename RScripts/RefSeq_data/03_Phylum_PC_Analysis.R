###########################################################
# 03. Phylum level PC analysis
###########################################################
# RScript to further analyze the metagenomic diet data at phylum level using PCA 

# 3.1 Normalise and transform data at phylum level
################################################################################

# 3.1.1 Use DESeq2 to normalize data per specimen
library(DESeq2)

# DeSeq2, Negative Binomial Method
# Make dataframes of Phyloseq object with summed samples per individual
phylum_otu <- as.data.frame(t(otu_table(physeq_phylum_Specimen)+1))
phylum_sample <- as.data.frame(sample_data(physeq_phylum_Specimen))
phylum_tax <- as.data.frame(tax_table(physeq_phylum_Specimen))
# Factor sample/metadata
phylum_sample$SpecimenID<-factor(phylum_sample$SpecimenID)
phylum_sample$TissueID<-factor(phylum_sample$TissueID)
phylum_sample$Sex<-factor(phylum_sample$Sex)
phylum_sample$SpeciesID<-factor(phylum_sample$SpeciesID)
# coldata needs to be ordered same way as counts for deseq
phylum_sample <- phylum_sample[order(rownames(phylum_sample)),]
phylum_otu <- phylum_otu[,order(colnames(phylum_otu))]
# Add metadata for DeSeq2
phylum_sample$SpeciesID = with(all_metadata_2, SpeciesID[match(rownames(phylum_sample),SpecimenID)])
phylum_sample$Tribe = with(all_metadata_2, Tribe[match(rownames(phylum_sample),SpecimenID)])
# Make a DeSeq object
dds <- DESeqDataSetFromMatrix(countData=phylum_otu, 
                              colData=phylum_sample, 
                              design=~SpeciesID) # Specifies how the counts from each OTU depend on our variables in the metadata, this case SpeciesID
dds
# Run DeSeq on object 
dds2 <- DESeq(dds)
# Check results
res <- results(dds2)
summary(res)
# LFC: Log Fold Cutoff
res[order(res$padj),]

# 3.1.2 Variance stabilizing transformation

# Then we need to transform the raw count data
# vst function will perform variance stabilizing transformation
vsdata <- varianceStabilizingTransformation(dds2)

# 3.2 Perform principal component analysis
################################################################################

# 3.2.1 Make PCA graph all specimens

# Make count data from normalised data
bf_pca_phylum <- data.frame(t(assay(vsdata)))
# Calculate PCA
vsdata_res_phylum <- prcomp((bf_pca_phylum), scale=F)
# Add metadata
bf_pca_phylum$Tribe = with(all_metadata_2, Tribe[match(rownames(bf_pca_phylum),all_metadata_2$SpecimenID)])
bf_pca_phylum$Sex = with(all_metadata_2, Sex[match(rownames(bf_pca_phylum),all_metadata_2$SpecimenID)])
bf_pca_phylum$SpeciesID = with(all_metadata_2, SpeciesID[match(rownames(bf_pca_phylum),all_metadata_2$SpecimenID)])
# Make an initial PCA plot with vectors for phyla driving the ordination
autoplot(vsdata_res_phylum,  data = bf_pca_phylum,colour='Tribe' ,loadings = TRUE, loadings.colour = 'gray24',loadings.label = T,  size=4,scale=F) +
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
# Define PCA without phyla names
Plot_PCA_vectors <- autoplot(vsdata_res_phylum, data = bf_pca_phylum,colour='Tribe' ,loadings = TRUE, loadings.colour = 'gray24',loadings.label = F,  size=4,scale=F) +
  coord_fixed() +
  scale_color_manual(values=tribe_colours) +
  xlab("PC1 (38.0%)") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x = element_text(color = "black", size = 11, face = "plain"),
        axis.text.y = element_text(color = "black", size = 11),  
        axis.title.y = element_blank(),
        legend.title = element_text(size=15), 
        legend.text = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0.2,0.2,0.5,0.2), "cm")) + 
  theme(aspect.ratio = 1)

# 3.2.2 Make PCA plot per species

# Extract PC's from PCA analysis
vsdata_res_phylum_x <- as.data.frame(vsdata_res_phylum$x)
# Add metadata
vsdata_res_phylum_x$Tribe = with(all_metadata_2, Tribe[match(rownames(vsdata_res_phylum_x),all_metadata_2$SpecimenID)])
vsdata_res_phylum_x$Sex = with(all_metadata_2, Sex[match(rownames(vsdata_res_phylum_x),all_metadata_2$SpecimenID)])
vsdata_res_phylum_x$SpeciesID = with(all_metadata_2, SpeciesID[match(rownames(vsdata_res_phylum_x),all_metadata_2$SpecimenID)])
# Make polygons of PC1 and PC2 per species
hull_cyl_phylum <- vsdata_res_phylum_x %>%
  group_by(SpeciesID) %>%
  dplyr::slice(chull(PC1, PC2))
hull_cyl_phylum$Tribe = with(all_metadata_2, Tribe[match(hull_cyl_phylum$SpeciesID,all_metadata_2$SpeciesID)])
# Make species list ordered by tribe for loop
spp <- vsdata_res_phylum_x[order(vsdata_res_phylum_x$Tribe),]
species <- unique(spp$SpeciesID)
# Define shape for sex
sex_shapes <- c("F"=25, "M"=24,"O"=20)
vsdata_res_phylum_x$Sex[is.na(vsdata_res_phylum_x$Sex)] <- "O"
# Loop for polygon in PCA per species
n <- list()
for (i in Species_ordering){
  n[[i]] <- ggplot() +
    geom_point(data=vsdata_res_phylum_x, aes(PC1, PC2, color = factor(SpeciesID)), size = 4, color="gray90") +
    coord_fixed() +
    theme_bw() +  
    aes(fill = Tribe) +
    theme(legend.position = "none") + 
    geom_polygon(data = hull_cyl_phylum[hull_cyl_phylum$SpeciesID == i, ], aes(PC1, PC2,colour=Tribe), alpha = 0.5) +
    geom_point(data=vsdata_res_phylum_x[vsdata_res_phylum_x$SpeciesID == i, ],aes(PC1, PC2,color=Tribe,shape=Sex), size=4) +
    scale_fill_manual(values=tribe_colours) +
    scale_color_manual(values=tribe_colours) +
    scale_shape_manual(values=sex_shapes) +
    ggtitle(i) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
}

Plot_pca_perspecies <- grid.arrange(grobs=n,ncol=8)

# 3.2.3 Make one PCA with a polygon for each species

ggplot() +
  geom_polygon(data = hull_cyl_phylum, aes(PC1, PC2,group=factor(SpeciesID)), alpha = 0.5) +
  coord_fixed() +
  theme_bw() +  
  aes(fill = Tribe, color=Tribe) +
  theme(legend.position = "none") + 
  scale_fill_manual(values=tribe_colours) +
  scale_color_manual(values=tribe_colours) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())

# 3.2.4 Make PCA plot with species mean and range

# Make graph with mean and range of each species within the PCA
pca_graph_mmm <- ddply(vsdata_res_phylum_x, .(SpeciesID), summarize, PCA1_mean=mean(PC1), PCA1_min=min(PC1),PCA1_max=max(PC1),PCA2_mean=mean(PC2), PCA2_min=min(PC2), PCA2_max=max(PC2))
pca_graph_mmm$Tribe = with(all_metadata_2, Tribe[match(pca_graph_mmm$SpeciesID,SpeciesID)])
# Plot the PCA with species ranges
Plot_PCA_ranges <- ggplot(pca_graph_mmm, aes(x=PCA1_mean, y=PCA2_mean, color=Tribe)) + 
  theme_minimal() +
  geom_point(size=4)+
  theme_bw() +
  xlab("PC1 (38.0%)") +
  scale_color_manual(values=tribe_colours) +
  geom_errorbar(aes(ymin=PCA2_min, ymax=PCA2_max), size=0.6) +
  geom_errorbarh(aes( xmin=PCA1_min, xmax=PCA1_max), size=0.6) +
  theme(axis.text.x = element_text(color = "black", size = 11, face = "plain"),
        axis.text.y = element_text(color = "black", size = 11),  
        axis.title.y = element_blank(),
        legend.title = element_text(size=15), 
        legend.text = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0.2,0.2,0.5,0.2), "cm")) + 
  theme(aspect.ratio = 1)

# 3.2.5 Plot species in PC3 and PC4

# Make graph with mean and range of each species within the PCA, PC3PC4
pca_graph_mmm_pc34 <- ddply(vsdata_res_phylum_x, .(SpeciesID), summarize, PCA3_mean=mean(PC3), PCA3_min=min(PC3),PCA3_max=max(PC3),PCA4_mean=mean(PC4), PCA4_min=min(PC4), PCA4_max=max(PC4))
pca_graph_mmm_pc34$Tribe = with(all_metadata_2, Tribe[match(pca_graph_mmm_pc34$SpeciesID,SpeciesID)])
# Plot the PCA with species ranges
pdf(file = "Plot_PCA_PC3-4_Phylum_rev_14012025.pdf",
    width = 10,
    height = 10) 
ggplot(pca_graph_mmm_pc34, aes(x=PCA3_mean, y=PCA4_mean, color=Tribe)) + 
  theme_minimal() +
  geom_point(size=6)+
  theme_bw() +
  xlab("PC3 (14%)") +
  ylab("PC4 (10%)") +
  scale_color_manual(values=tribe_colours) +
  geom_errorbar(aes(ymin=PCA4_min, ymax=PCA4_max), size=0.6) +
  geom_errorbarh(aes( xmin=PCA3_min, xmax=PCA3_max), size=0.6) +
  theme(axis.text.x = element_text(color = "black", size = 11, face = "plain"),
        axis.text.y = element_text(color = "black", size = 11),  
        legend.title = element_text(size=15), 
        legend.text = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0.2,0.2,0.5,0.2), "cm")) + 
  theme(aspect.ratio = 1.1)
dev.off()

# 3.2.6 Compare normalized PC1-4 distance to genetic distance

# PC1-4 distance from main figure calculated below
PC_dist <- read.csv("PC1_4_dist_norm.csv")
# Define variables for plot
all(phylo_dist_melt_NoOrAs$Var1 == PC_dist$Var1)
all(phylo_dist_melt_NoOrAs$Var2 == PC_dist$Var2)
PC_dist$Tribe <- phylo_dist_melt_NoOrAs$Tribe
# Plot PC1-4 distance with genetic distance
plot_dist_pca_gen <- ggplot() +
  geom_point(data=PC_dist, aes(gen_dist, value,colour=Tribe), size=4) +
  scale_color_manual(values=tribe_colours) + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  xlab("Genetic distance") +
  #ylab("Distance PCA diet (PC1-4)") +
  theme(axis.text.x = element_text(color = "black", size = 11, face = "plain"),
        axis.text.y = element_text(color = "black", size = 11),  
        #axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size=15), 
        legend.text = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        aspect.ratio = 1,
        plot.margin = unit(c(0.2,0.2,0.5,0.2), "cm"))+ 
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1))

# 3.2.7 Plot main figure 3

# Plot PCA's with genetic distance
Plot_PCA_1 <- grid.arrange(Plot_PCA_ranges,Plot_PCA_vectors,nrow=1)
pdf(file = "Plot_PCA_Phylum_rev_v2_15012025.pdf",
    width = 12,
    height = 12) 
grid.arrange(Plot_PCA_vectors,Plot_PCA_ranges,plot_dist_pca_gen,nrow=2)
dev.off()

# 3.3 Perform PCA only for LT radiation
################################################################################

# 3.3.1 Normalise and transform data without Astbur and Orenil

# Remove Oretan and Astbur since not part of LT adaptive radiation
OreAst <- rownames(phylum_sample[ phylum_sample$SpeciesID %in% c("Oretan","Astbur"), ])
phylum_sample_noOrAs <- phylum_sample[ ! rownames(phylum_sample) %in% OreAst, ]
phylum_otu_noOrAs <- phylum_otu[ ,! colnames(phylum_otu) %in% OreAst ]
# Make a DeSeq object
dds_noOrAs <- DESeqDataSetFromMatrix(countData=phylum_otu_noOrAs, 
                                     colData=phylum_sample_noOrAs, 
                                     design=~SpeciesID)
# Perform DeSeq
dds2_noOrAs <- DESeq(dds_noOrAs)
# Perform vst
vsdata_noOrAs <- varianceStabilizingTransformation(dds2_noOrAs)
# Extract normalised count data
bf_pca_phylum_noOrAs <- data.frame(t(assay(vsdata_noOrAs)))

# 3.3.2 Perform PCA

# Perform PCA analysis
vsdata_res_phylum_noOrAs <- prcomp((bf_pca_phylum_noOrAs), scale=F)
# Add metadata
bf_pca_phylum_noOrAs$Tribe = with(all_metadata_2, Tribe[match(rownames(bf_pca_phylum_noOrAs),all_metadata_2$SpecimenID)])
bf_pca_phylum_noOrAs$Sex = with(all_metadata_2, Sex[match(rownames(bf_pca_phylum_noOrAs),all_metadata_2$SpecimenID)])
bf_pca_phylum_noOrAs$SpeciesID = with(all_metadata_2, SpeciesID[match(rownames(bf_pca_phylum_noOrAs),all_metadata_2$SpecimenID)])
# Make PCA plot
PCA_plot_isopes <- autoplot(vsdata_res_phylum_noOrAs,  data = bf_pca_phylum_noOrAs,colour='Tribe' ,loadings = F,  size=4,scale=F) +
  coord_fixed() +
  scale_color_manual(values=tribe_colours) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.position = "none") +
  ggtitle("PCA Diet Phylum lvl") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(aspect.ratio = 1.5)

# 3.4 Compare PCA to species' phenotype
################################################################################

# 3.4.1 Mean values per species

# Extract PCA
vsdata_res_phylum_x_noOrAs <- as.data.frame(vsdata_res_phylum_noOrAs$x)
vsdata_res_phylum_x_noOrAs$SpeciesID = with(all_metadata_2, SpeciesID[match(rownames(vsdata_res_phylum_x_noOrAs),all_metadata_2$SpecimenID)])
# Make PCA dataframe with species means PC1 and PC2
pca_graph <- ddply(vsdata_res_phylum_x_noOrAs, .(SpeciesID), summarize, PCA1_mean=mean(PC1), PCA1_min=min(PC1),PCA1_max=max(PC1),PCA2_mean=mean(PC2), PCA2_min=min(PC2), PCA2_max=max(PC2))
# Add species means of other metadata
pca_graph$d15N_mean = with(all_metadata_2, d15N_mean[match(pca_graph$SpeciesID,SpeciesID)])
pca_graph$d15N_std_dev = with(all_metadata_2, d15N_std_dev[match(pca_graph$SpeciesID,SpeciesID)])
pca_graph$d13C_mean = with(all_metadata_2, d13C_mean[match(pca_graph$SpeciesID,SpeciesID)])
pca_graph$d13C_std_dev = with(all_metadata_2, d13C_std_dev[match(pca_graph$SpeciesID,SpeciesID)])
pca_graph$Tribe = with(all_metadata_2, Tribe[match(pca_graph$SpeciesID,SpeciesID)])
pca_graph$FoodCat = with(all_metadata_2, FoodCat[match(pca_graph$SpeciesID,SpeciesID)])
pca_graph$LPJ_PC1_mean = with(LPJ_vars, LPJ_PC1_mean[match(pca_graph$SpeciesID,sp)])
pca_graph$LPJ_PC1_std_dev = with(LPJ_vars, LPJ_PC1_std_dev[match(pca_graph$SpeciesID,sp)])
pca_graph$LPJ_PC2_mean = with(LPJ_vars, LPJ_PC2_mean[match(pca_graph$SpeciesID,sp)])
pca_graph$LPJ_PC2_std_dev = with(LPJ_vars, LPJ_PC2_std_dev[match(pca_graph$SpeciesID,sp)])
pca_graph$UOJ_PC1_mean = with(UOJ_vars, UOJ_PC1_mean[match(pca_graph$SpeciesID,sp)])
pca_graph$UOJ_PC1_std_dev = with(UOJ_vars, UOJ_PC1_std_dev[match(pca_graph$SpeciesID,sp)])
pca_graph$UOJ_PC2_mean = with(UOJ_vars, UOJ_PC2_mean[match(pca_graph$SpeciesID,sp)])
pca_graph$UOJ_PC2_std_dev = with(UOJ_vars, UOJ_PC2_std_dev[match(pca_graph$SpeciesID,sp)])
pca_graph$SpeciesID = with(all_metadata_2, SpeciesID[match(pca_graph$SpeciesID,SpeciesID)])
pca_graph$BDY_PC1_mean = with(BDY_vars, BDY_PC1_mean[match(pca_graph$SpeciesID,sp)])
pca_graph$BDY_PC1_std_dev = with(BDY_vars, BDY_PC1_std_dev[match(pca_graph$SpeciesID,sp)])
pca_graph$BDY_PC2_mean = with(BDY_vars, BDY_PC2_mean[match(pca_graph$SpeciesID,sp)])
pca_graph$BDY_PC2_std_dev = with(BDY_vars, BDY_PC2_std_dev[match(pca_graph$SpeciesID,sp)])
pca_graph$SL_mean = with(SL_vars, SL_mean[match(pca_graph$SpeciesID,SpeciesID)])
pca_graph$SL_std_dev = with(SL_vars, SL_std_dev[match(pca_graph$SpeciesID,SpeciesID)])
pca_graph$Weight_mean = with(Weigth_vars, Weight_mean[match(pca_graph$SpeciesID,SpeciesID)])
pca_graph$ZHL_mean_CETH = with(ZHL_vars_CETH, ZHL_mean[match(pca_graph$SpeciesID,SpeciesID)])
pca_graph$ZHL_std_dev_CETH = with(ZHL_vars_CETH, ZHL_std_dev[match(pca_graph$SpeciesID,SpeciesID)])
pca_graph$RGL_mean_CETH = with(RGL_vars_CETH, RGL_mean[match(pca_graph$SpeciesID,SpeciesID)])
pca_graph$RGL_std_dev_CETH = with(RGL_vars_CETH, RGL_std_dev[match(pca_graph$SpeciesID,SpeciesID)])

# 3.4.2 Comparison species' PCA and metadata

# Compare PCA axes with other species variables
# Make list of variables to compare to PC1 and PC2
pc_vars <- c("d15N_mean","d13C_mean","LPJ_PC1_mean","LPJ_PC2_mean","UOJ_PC1_mean","UOJ_PC2_mean","BDY_PC1_mean","BDY_PC2_mean","ZHL_mean_CETH","RGL_mean_CETH","SL_mean","Weight_mean")
# Make loop for PC1 vs metadata
n_pc1 <- list()
for (i in pc_vars){
  n_pc1[[i]] = ggplot(data = pca_graph, aes_string(x="PCA1_mean", y=i, color="Tribe"))+ # aes_string because of string variable i
    geom_point() + 
    theme_bw() +
    xlab("PC1 (37.56%)") +
    scale_color_manual(values=tribe_colours) +
    theme(legend.position = "none") + 
    theme(axis.text.x = element_text(color = "black", size = 11, face = "plain"),
          axis.text.y = element_text(color = "black", size = 11),  
          axis.title.x = element_text(color = "black", size = 15, hjust = .5, vjust = 0),
          axis.title.y = element_text(color = "black", size = 15, vjust = 2),
          legend.title = element_text(size=15), legend.text = element_text(size=12),
          panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    stat_smooth(method = "lm", col = "darkgrey") 
}
iso_pc1_plot <- grid.arrange(grobs=n_pc1,ncol=4)
# Make loop for PC2 vs metadata
n_pc2 <- list()
for (i in pc_vars){
  n_pc2[[i]] = ggplot(data = pca_graph, aes_string(x="PCA2_mean", y=i, color="Tribe"))+
    geom_point() +
    theme_bw() +
    xlab("PC2 (15.84%)") +
    scale_color_manual(values=tribe_colours) +
    theme(legend.position = "none") + 
    theme(axis.text.x = element_text(color = "black", size = 11, face = "plain"),
          axis.text.y = element_text(color = "black", size = 11),  
          axis.title.x = element_text(color = "black", size = 15, hjust = .5, vjust = 0),
          axis.title.y = element_text(color = "black", size = 15, vjust = 2),
          legend.title = element_text(size=15), legend.text = element_text(size=12),
          panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    stat_smooth(method = "lm", col = "darkgrey") 
}
iso_pc2_plot <- grid.arrange(grobs=n_pc2,ncol=4)
iso_pc12_plot <- grid.arrange(iso_pc1_plot,iso_pc2_plot,ncol=1)
pdf(file = "Plot_PCA_variables_19092024.pdf",
    width = 15,
    height = 10)
grid.arrange(PCA_plot_isopes,iso_pc12_plot, ncol=2 )
dev.off()

# 3.4.3 Perform PGLS with PCA
library(caper)

# pgls taking into account phylogeny
# Do this to avoid abnormal error for BDY_PC2 and PC1: 52ERROR: ABNORMAL_TERMINATION_IN_LNSRCH
pca_graph0 <- pca_graph
pca_graph0$BDY_PC2_mean <- pca_graph0$BDY_PC2_mean*10
# Then create comparative dataset of species
compset = comparative.data(pruned.tree, pca_graph0, SpeciesID)
# Look at pgls results for PCA vs metadata
# Each time trait as a result of diet, in this order
print(summary(pgls(d15N_mean ~ PCA1_mean, compset, lambda= 'ML')))
print(summary(pgls(d15N_mean ~ PCA2_mean, compset, lambda= 'ML')))
print(summary(pgls(d13C_mean ~ PCA1_mean, compset, lambda= 'ML')))
print(summary(pgls(d13C_mean ~ PCA2_mean, compset, lambda= 'ML')))
print(summary(pgls(LPJ_PC1_mean ~ PCA1_mean, compset, lambda= 'ML')))
print(summary(pgls(LPJ_PC1_mean ~ PCA2_mean, compset, lambda= 'ML')))
print(summary(pgls(LPJ_PC2_mean ~ PCA1_mean, compset, lambda= 'ML')))
print(summary(pgls(LPJ_PC2_mean ~ PCA2_mean, compset, lambda= 'ML')))
print(summary(pgls(UOJ_PC1_mean ~ PCA1_mean, compset, lambda= 'ML')))
print(summary(pgls(UOJ_PC1_mean ~ PCA2_mean, compset, lambda= 'ML')))
print(summary(pgls(UOJ_PC2_mean ~ PCA1_mean, compset, lambda= 'ML')))
print(summary(pgls(UOJ_PC2_mean ~ PCA2_mean, compset, lambda= 'ML')))
print(summary(pgls(BDY_PC1_mean ~ PCA1_mean, compset, lambda= 'ML')))
print(summary(pgls(BDY_PC1_mean ~ PCA2_mean, compset, lambda= 'ML')))
print(summary(pgls(BDY_PC2_mean ~ PCA1_mean, compset, lambda= 'ML')))
print(summary(pgls(BDY_PC2_mean ~ PCA2_mean, compset, lambda= 'ML')))
print(summary(pgls(ZHL_mean_CETH ~ PCA1_mean, compset, lambda= 'ML')))
print(summary(pgls(ZHL_mean_CETH ~ PCA2_mean, compset, lambda= 'ML')))
print(summary(pgls(RGL_mean_CETH ~ PCA1_mean, compset, lambda= 'ML')))
print(summary(pgls(RGL_mean_CETH ~ PCA2_mean, compset, lambda= 'ML')))
print(summary(pgls(SL_mean ~ PCA1_mean, compset, lambda= 'ML')))
print(summary(pgls(SL_mean ~ PCA2_mean, compset, lambda= 'ML')))
print(summary(pgls(Weight_mean ~ PCA1_mean, compset, lambda= 'ML')))
print(summary(pgls(Weight_mean ~ PCA2_mean, compset, lambda= 'ML')))


