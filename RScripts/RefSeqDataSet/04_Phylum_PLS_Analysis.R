###########################################################
# 04. Phylum level PLS analysis
###########################################################
# RScript to further analyze the metagenomic diet data at phylum level using PLS 

# 4.1 Normalise and transform data at phylum level
################################################################################

# 4.1.1 Keep only species values without Oretan and Astbur

# Do transformations of species together
physeq_phylum_Specimen_spp_rel
# Make dataframes out of Physeq
phylum_otu_spp <- as.data.frame(t(otu_table(physeq_phylum_Specimen_spp_rel))+1)
phylum_sample_spp <- as.data.frame(t(sample_data(physeq_phylum_Specimen_spp_rel)))
# Remove Oretan and Astbur from comparative analysis
phylum_otu_spp <- phylum_otu_spp[ ! colnames(phylum_otu_spp) %in% c("Oretan","Astbur")]
phylum_sample_spp <- phylum_sample_spp[ ! colnames(phylum_sample_spp) %in% c("Oretan","Astbur")]
phylum_sample_spp <- as.data.frame(t(phylum_sample_spp))
# Order the two datasets similarly
phylum_sample_spp <- phylum_sample_spp[order(rownames(phylum_sample_spp)),]
phylum_otu_spp <- phylum_otu_spp[,order(colnames(phylum_otu_spp))]
# Need integers for DESeq
phylum_otu_spp <- round(phylum_otu_spp, digits = 0)
# Add tribe
phylum_sample_spp$Tribe <- as.factor(phylum_sample_spp$Tribe)
# Remove Orenil and Astbur from phylogeny
pruned.tree_noOrAs = drop.tip(pruned.tree, pruned.tree$tip.label[-which(pruned.tree$tip.label %in% rownames(pls_diet))])

# 4.1.2 Perform DESeq

# Make DESeq object
dds_spp <- DESeqDataSetFromMatrix(countData=phylum_otu_spp, 
                                  colData=phylum_sample_spp, 
                                  design=~Tribe)
# Perform DESeq
dds2_spp <- DESeq(dds_spp)
# Do vst
vsdata_spp <- varianceStabilizingTransformation(dds2_spp)
# Make matrix
pls_diet <- as.matrix(t(assay(vsdata_spp)))
# Replace colnames by diet phylum
colnames(pls_diet) <- with(Tax_table_strict_2, phylum[match(colnames(pls_diet),rownames(Tax_table_strict_2))])
pls_diet <- as.matrix(pls_diet)

# 4.2 Two-block PLS and Isotopes
################################################################################
library(geomorph)

# 4.2.1 Without accounting for phylogenetic non-independence

# Make matrix
pls_diet_iso <- as.matrix(t(assay(vsdata_spp)))
# Replace colnames by diet phylum
colnames(pls_diet_iso) <- with(Tax_table_strict_2, phylum[match(colnames(pls_diet_iso),rownames(Tax_table_strict_2))])
pls_diet_iso <- as.matrix(pls_diet_iso)
# Get isotope means
iso_subset <- pca_graph[c("SpeciesID", "d15N_mean", "d13C_mean" )]
rownames(iso_subset) <- iso_subset$SpeciesID
iso_subset <- iso_subset[c( "d15N_mean", "d13C_mean" )]
# Scale isotopes to eachother
pls_iso <- scale(iso_subset)
# Order in the same way
pls_diet_iso <- pls_diet_iso[match(rownames(pls_iso), rownames(pls_diet_iso)),]
all(rownames(pls_iso) == rownames(pls_diet_iso))
# Run 2-block PLS
f_iso = two.b.pls(pls_diet_iso, pls_iso) 
# Summary
summary(f_iso)
# Check loading of the projections
f_iso$left.pls.vectors
f_iso$right.pls.vectors
# Get projected scores
f1_iso= data.frame( "XScores"=f_iso$ XScores[,1], "YScores"=f_iso$ YScores[,1])
f1_iso$Tribe = with(all_metadata_2, Tribe[match(rownames(f1_iso),SpeciesID)])
# R^2 for the PLS fit
summary(lm(f1_iso$YScores ~ f1_iso$XScores))
# pGLS
f1_iso$SpeciesID <- rownames(f1_iso)
compset_iso = comparative.data(pruned.tree, f1_iso, SpeciesID)
# Look at pgls results for PCA vs metadata
# YScores need to be before XScores
print(summary(pgls(YScores ~ XScores, compset_iso, lambda= 'ML')))
# Plot PLS
plot_pls_iso <- ggplot(data = f1_iso, aes(XScores, YScores,colour=Tribe,size=3)) +
  geom_point() +
  theme_bw() +  
  scale_color_manual(values=tribe_colours) +
  ylab("Stable isotope projection") +
  xlab("food source projection") +
  stat_smooth(method = "lm", col = "darkgray",se = FALSE,size=1, linetype=2) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text.x = element_text(color = "black", size = 11, face = "plain"),
        axis.text.y = element_text(color = "black", size = 11),  
        axis.title.x = element_text(color = "black", size = 15, hjust = .5, vjust = 0),
        axis.title.y = element_text(color = "black", size = 15, vjust = 2),
        legend.title = element_text(size=15), legend.text = element_text(size=12),
        legend.position="none",
        plot.margin = unit(c(0.2,0.2,0.2,1), "cm"))

# 4.2.2 With accounting for phylogenetic non-independence

# Order datasets by phylogeny
pls_diet_ord_iso <- pls_diet_iso[ match(pruned.tree_noOrAs$tip.label, rownames(pls_diet_iso)),]
pls_ord_iso <- pls_iso[ match(pruned.tree_noOrAs$tip.label, rownames(pls_iso)),]
# Check
all(pruned.tree_noOrAs$tip.label == rownames(pls_diet_ord_iso)) 
all(pruned.tree_noOrAs$tip.label == rownames(pls_ord_iso))
# Do PLS
fp_iso = phylo.integration(  A = pls_diet_ord_iso,  A2 = pls_ord_iso,  phy = pruned.tree_noOrAs  )
# Summary
summary(fp_iso)
# Check loading of the projections
fp_iso$left.pls.vectors
fp_iso$right.pls.vectors
# Get projected scores
f1p_iso= data.frame( "XScores"=fp_iso$ XScores[,1], "YScores"=fp_iso$ YScores[,1])
f1p_iso$Tribe = with(all_metadata_2, Tribe[match(rownames(f1p_iso),SpeciesID)])
# R^2 for the PLS fit
summary(lm(f1p_iso$YScores ~ f1p_iso$XScores))
# plot
plot_pls_phy_iso <- ggplot(data = f1p_iso, aes(XScores, YScores,colour=Tribe,size=3)) +
  geom_point() +
  theme_bw() +  
  scale_color_manual(values=tribe_colours) +
  ylab("Stable isotope projection") +
  xlab("food source projection") +
  stat_smooth(method = "lm", col = "darkgray",se = FALSE,size=1, linetype=2) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text.x = element_text(color = "black", size = 11, face = "plain"),
        axis.text.y = element_text(color = "black", size = 11),  
        axis.title.x = element_text(color = "black", size = 15, hjust = .5, vjust = 0),
        axis.title.y = element_text(color = "black", size = 15, vjust = 2),
        legend.title = element_text(size=15), legend.text = element_text(size=12),
        legend.position="none",
        plot.margin = unit(c(0.2,0.2,0.2,1), "cm"))

# 4.3 Two-block PLS and Body Morphology
################################################################################

# 4.3.1 PLS without accounting for phylogenetic non-independence

# Load body shape data by Ronco et al
load("sp_means_body.Rdata")
pls_bdy <- sp_means_body
# Subset to species study
pls_bdy <- pls_bdy[ ,, dimnames(pls_bdy)[[3]] %in% rownames(pls_diet) ]
head(sp_means_body_subset)
# Order in the same way
pls_diet <- pls_diet[match(dimnames(pls_bdy)[[3]], rownames(pls_diet)),]
all(dimnames(pls_bdy)[[3]] == rownames(pls_diet))
# Run 2-block PLS
f = two.b.pls(pls_diet, pls_bdy) 
# Plot
P = plot(f)
# Summary
summary(f)
# Check loading of the projections
f$left.pls.vectors
f$right.pls.vectors
# Get projected scores
f1= data.frame( "XScores"=f$ XScores[,1], "YScores"=f$ YScores[,1])
f1$Tribe = with(all_metadata_2, Tribe[match(rownames(f1),SpeciesID)])
# R^2 for the PLS fit
summary(lm(f1$YScores ~ f1$XScores))
# pGLS
f1$SpeciesID <- rownames(f1)
compset_bdy = comparative.data(pruned.tree, f1, SpeciesID)
print(summary(pgls(YScores ~ XScores , compset_bdy, lambda= 'ML')))

# 4.3.2 Plot PLS without accounting for phylogenetic non-independence

# Visualize shape at minimum and maximum PLS scores.
minx <- min(P$plot_args$y)
maxx <- max(P$plot_args$y)
preds <- shape.predictor(P$A2, 
                         x = P$plot.args$y,
                         min = minx, max = maxx)
# Plot shape
plotRefToTarget(mshape(P$A2), preds$min)
plotRefToTarget(mshape(P$A2), preds$max)
# Define shapes
minshape = preds$min
maxshape = preds$max
# Plot PLS
plot_pls_bdy <- ggplot(data = f1, aes(XScores, YScores,colour=Tribe,size=3)) +
  geom_point() +
  theme_bw() +  
  scale_color_manual(values=tribe_colours) +
  ylab("body shape projection") +
  xlab("food source projection") +
  stat_smooth(method = "lm", col = "darkgray",se = FALSE,size=1, linetype=2) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text.x = element_text(color = "black", size = 11, face = "plain"),
        axis.text.y = element_text(color = "black", size = 11),  
        axis.title.x = element_text(color = "black", size = 15, hjust = .5, vjust = 0),
        axis.title.y = element_text(color = "black", size = 15, vjust = 2),
        legend.title = element_text(size=15), legend.text = element_text(size=12),
        legend.position="none",
        plot.margin = unit(c(0.2,0.2,0.2,1), "cm"))
# outlines to draw the shapes
outline1 = c(16,17,18,19,16)
outline2 = c(1,2,3,4,5,6,8,9,11,12,15,1)
outline3 = c(14,10,7)
# Plot the max shape
pdf(file = "Plot_BDY_NoPhy_max_10012025.pdf",
    width = 7,
    height = 7) 
plot(maxshape[,1], maxshape[,2], xlim=c(-0.4,0.4), ylim=c(-0.4,0.4), pch=16, cex=0.3, axes=F, xlab="", ylab="", asp=1)
polygon(maxshape[c(outline1),],lwd=4)
lines(maxshape[c(outline2),],lwd=4)
lines(maxshape[c(outline3),],lwd=4)
dev.off()
# Plot min shape
pdf(file = "Plot_BDY_NoPhy_min_10012025.pdf",
    width = 7,
    height = 7) 
plot(minshape[,1],minshape[,2], xlim=c(-0.4,0.4), ylim=c(-0.4,0.4), pch=16, cex=0.3, axes=F, xlab="", ylab="", asp=1)
polygon(minshape[c(outline1),],lwd=4)
lines(minshape[c(outline2),],lwd=4)
lines(minshape[c(outline3),],lwd=4)
dev.off()

# 4.3.3 PLS with accounting for phylogenetic non-independence

# Account in PLS for phylogeny
# Order other comparative datasets by order of tree tips
pls_diet_ord <- pls_diet[ match(pruned.tree_noOrAs$tip.label, rownames(pls_diet)),]
pls_bdy_ord <- pls_bdy[,, match(pruned.tree_noOrAs$tip.label, dimnames(pls_bdy)[[3]])]
# Check if right order
all(pruned.tree_noOrAs$tip.label == rownames(pls_diet_ord)) 
all(pruned.tree_noOrAs$tip.label == dimnames(pls_bdy_ord)[[3]])
# Perform PLS analysis
fp = phylo.integration(  A = pls_diet_ord,  A2 = pls_bdy_ord,  phy = pruned.tree_noOrAs  )
# Plot PLS
Pp = plot(fp)
# Summary PLS
summary(fp)
# Check loading of the projections
fp$left.pls.vectors
# Make dataframe of scores and add tribe
f1p= data.frame( "XScores"=fp$ XScores[,1], "YScores"=fp$ YScores[,1])
f1p$Tribe = with(all_metadata_2, Tribe[match(rownames(f1p),SpeciesID)])
# R^2 for the PLS fit
summary(lm(f1p$YScores ~ f1p$XScores))

# 4.3.4 Plot PLS with accounting for phylogenetic non-independence

# Define min and max shapes
minxp <- min(Pp$plot_args$y)
maxxp <- max(Pp$plot_args$y)
predsp <- shape.predictor(Pp$A2, 
                          x = Pp$plot.args$y,
                          min = minxp, max = maxxp)
plotRefToTarget(mshape(Pp$A2), predsp$min)
plotRefToTarget(mshape(Pp$A2), predsp$max)
minshapep = predsp$min
maxshapep = predsp$max
# Plot PLS
plot_pls_phy_bdy <- ggplot(data = f1p, aes(XScores, YScores,colour=Tribe,size=3)) +
  geom_point() +
  theme_bw() +  
  scale_color_manual(values=tribe_colours) +
  ylab("body shape projection") +
  xlab("food source projection") +
  stat_smooth(method = "lm", col = "darkgray",se = FALSE,size=1, linetype=2) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text.x = element_text(color = "black", size = 11, face = "plain"),
        axis.text.y = element_text(color = "black", size = 11),  
        axis.title.x = element_text(color = "black", size = 15, hjust = .5, vjust = 0),
        axis.title.y = element_text(color = "black", size = 15, vjust = 2),
        legend.title = element_text(size=15), legend.text = element_text(size=12),
        legend.position="none",
        plot.margin = unit(c(0.2,0.2,0.2,1), "cm"))
# outlines to draw the shapes
outline1 = c(16,17,18,19,16)
outline2 = c(1,2,3,4,5,6,8,9,11,12,15,1)
outline3 = c(14,10,7)
# Plot max shape
pdf(file = "Plot_BDY_Phy_max_13012025.pdf",
    width = 7,
    height = 7) 
plot(maxshapep[,1], maxshapep[,2], xlim=c(-0.4,0.4), ylim=c(-0.4,0.4), pch=16, cex=0.3, axes=F, xlab="", ylab="", asp=1)
polygon(maxshapep[c(outline1),],lwd=4)
lines(maxshapep[c(outline2),],lwd=4)
lines(maxshapep[c(outline3),],lwd=4)
dev.off()
# Plot min shape
pdf(file = "Plot_BDY_Phy_min_13012025.pdf",
    width = 7,
    height = 7) 
plot(minshapep[,1],minshapep[,2], xlim=c(-0.4,0.4), ylim=c(-0.4,0.4), pch=16, cex=0.3, axes=F, xlab="", ylab="", asp=1)
polygon(minshapep[c(outline1),],lwd=4)
lines(minshapep[c(outline2),],lwd=4)
lines(minshapep[c(outline3),],lwd=4)
dev.off()

# 4.4 Two-block PLS and Upper Oral Jaw Shape
################################################################################

# 4.4.1 PLS without accounting for phylogenetic non-independence

# Make matrix
pls_diet_uoj <- as.matrix(t(assay(vsdata_spp)))
# Replace colnames by diet phylum
colnames(pls_diet_uoj) <- with(Tax_table_strict_2, phylum[match(colnames(pls_diet_uoj),rownames(Tax_table_strict_2))])
pls_diet_uoj <- as.matrix(pls_diet_uoj)
# Load data Ronco et al
load("sp_means_OJ.Rdata")
pls_uoj <- sp_means_OJ
# Subset
pls_uoj <- pls_uoj[ ,, dimnames(pls_uoj)[[3]] %in% rownames(pls_diet) ]
# Order in the same way
pls_diet_uoj <- pls_diet_uoj[match(dimnames(pls_uoj)[[3]], rownames(pls_diet_uoj)),]
all(dimnames(pls_uoj)[[3]] == rownames(pls_diet_uoj))
# Run 2-block PLS
f_uoj = two.b.pls(pls_diet_uoj, pls_uoj) 
# Plot PLS
P_uoj = plot(f_uoj)
# Summary
summary(f_uoj)
# check loading of the projections
f_uoj$left.pls.vectors
f_uoj$right.pls.vectors
# Get projected scores
f1_uoj= data.frame( "XScores"=f_uoj$ XScores[,1], "YScores"=f_uoj$ YScores[,1])
f1_uoj$Tribe = with(all_metadata_2, Tribe[match(rownames(f1_uoj),SpeciesID)])
# R^2 for the PLS fit
summary(lm(f1_uoj$YScores ~ f1_uoj$XScores))
# pGLS
f1_uoj$SpeciesID <- rownames(f1_uoj)
compset_uoj = comparative.data(pruned.tree, f1_uoj, SpeciesID)
print(summary(pgls(YScores ~ XScores , compset_uoj, lambda= 'ML')))

# 4.4.2 Plot PLS without accounting for phylogenetic non-independence

# Visualize shape at minimum and maximum PLS scores.
minx_uoj <- min(P_uoj$plot_args$y)
maxx_uoj <- max(P_uoj$plot_args$y)
preds_uoj <- shape.predictor(P_uoj$A2, 
                             x = P_uoj$plot.args$y,
                             min = minx_uoj, max = maxx_uoj)
plotRefToTarget(mshape(P_uoj$A2), preds_uoj$min)
plotRefToTarget(mshape(P_uoj$A2), preds_uoj$max)
minshape_uoj = preds_uoj$min
maxshape_uoj = preds_uoj$max
# Plot PLS
plot_pls_uoj <-ggplot(data = f1_uoj, aes(XScores, YScores,colour=Tribe,size=3)) +
  geom_point() +
  theme_bw() +  
  scale_color_manual(values=tribe_colours) +
  ylab("upper oral jaw shape projection") +
  xlab("food source projection") +
  stat_smooth(method = "lm", col = "darkgray",se = FALSE,size=1, linetype=2) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text.x = element_text(color = "black", size = 11, face = "plain"),
        axis.text.y = element_text(color = "black", size = 11),  
        axis.title.x = element_text(color = "black", size = 15, hjust = .5, vjust = 0),
        axis.title.y = element_text(color = "black", size = 15, vjust = 2),
        legend.title = element_text(size=15), legend.text = element_text(size=12),
        legend.position="none",
        plot.margin = unit(c(0.2,0.2,0.2,1), "cm"))
# Plot min shape
pdf(file = "Plot_UOJ_NoPhy_max_10012025.pdf",
    width = 7,
    height = 7) 
plot(maxshape_uoj[,1],maxshape_uoj[,2], xlim=c(-0.1,0.1), ylim=c(-0.1,0.1), pch=16, cex=0.3, axes=F, xlab="", ylab="", asp=1)
lines(maxshape_uoj[c(3,1,4,2),],lwd=4)
dev.off()
# Plot min shape
pdf(file = "Plot_UOJ_NoPhy_min_10012025.pdf",
    width = 7,
    height = 7) 
plot(minshape_uoj[,1],minshape_uoj[,2], xlim=c(-0.1,0.1), ylim=c(-0.1,0.1), pch=16, cex=0.3, axes=F, xlab="", ylab="", asp=1)
lines(minshape_uoj[c(3,1,4,2),],lwd=4)
dev.off()

# 4.4.3 PLS with accounting for phylogenetic non-independence

# Order by phylogenetic tree
pls_diet_ord_uoj <- pls_diet_uoj[ match(pruned.tree_noOrAs$tip.label, rownames(pls_diet_uoj)),]
pls_bdy_ord_uoj <- pls_uoj[,, match(pruned.tree_noOrAs$tip.label, dimnames(pls_uoj)[[3]])]
# Check
all(pruned.tree_noOrAs$tip.label == rownames(pls_diet_ord_uoj)) 
all(pruned.tree_noOrAs$tip.label == dimnames(pls_bdy_ord_uoj)[[3]])
# Do PLS
fp_uoj = phylo.integration(  A = pls_diet_ord_uoj,  A2 = pls_bdy_ord_uoj,  phy = pruned.tree_noOrAs  )
# Plot PLS
Pp_uoj = plot(fp_uoj)
# Summary
summary(fp_uoj)
# Check loading of the projections
fp_uoj$left.pls.vectors
fp_uoj$right.pls.vectors
# Get projected scores
f1p_uoj= data.frame( "XScores"=fp_uoj$ XScores[,1], "YScores"=fp_uoj$ YScores[,1])
f1p_uoj$Tribe = with(all_metadata_2, Tribe[match(rownames(f1p_uoj),SpeciesID)])
# R^2 for the PLS fit
summary(lm(f1p_uoj$YScores ~ f1p_uoj$XScores))

# 4.4.4 Plot PLS with accounting for phylogenetic non-independence

# Visualize shape at minimum and maximum PLS scores.
minxp_uoj <- min(Pp_uoj$plot_args$y)
maxxp_uoj <- max(Pp_uoj$plot_args$y)
predsp_uoj <- shape.predictor(Pp_uoj$A2, 
                              x = Pp_uoj$plot.args$y,
                              min = minxp_uoj, max = maxxp_uoj)
plotRefToTarget(mshape(Pp_uoj$A2), predsp_uoj$min)
plotRefToTarget(mshape(Pp_uoj$A2), predsp_uoj$max)
minshapep_uoj = predsp_uoj$min
maxshapep_uoj = predsp_uoj$max
# Plot PLS
plot_pls_phy_uoj <- ggplot(data = f1p_uoj, aes(XScores, YScores,colour=Tribe,size=3)) +
  geom_point() +
  theme_bw() +  
  scale_color_manual(values=tribe_colours) +
  ylab("upper oral jaw shape projection") +
  xlab("food source projection") +
  stat_smooth(method = "lm", col = "darkgray",se = FALSE,size=1, linetype=2) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text.x = element_text(color = "black", size = 11, face = "plain"),
        axis.text.y = element_text(color = "black", size = 11),  
        axis.title.x = element_text(color = "black", size = 15, hjust = .5, vjust = 0),
        axis.title.y = element_text(color = "black", size = 15, vjust = 2),
        legend.title = element_text(size=15), legend.text = element_text(size=12),
        legend.position="none",
        plot.margin = unit(c(0.2,0.2,0.2,1), "cm"))
# Plot min shape
pdf(file = "Plot_UOJ_Phy_max_13012025.pdf",
    width = 7,
    height = 7) 
plot(maxshapep_uoj[,1],maxshapep_uoj[,2], xlim=c(-0.1,0.1), ylim=c(-0.1,0.1), pch=16, cex=0.3, axes=F, xlab="", ylab="", asp=1)
lines(maxshapep_uoj[c(3,1,4,2),],lwd=4)
dev.off()
# Plot min shape
pdf(file = "Plot_UOJ_Phy_min_13012025.pdf",
    width = 7,
    height = 7) 
plot(minshapep_uoj[,1],minshapep_uoj[,2], xlim=c(-0.1,0.1), ylim=c(-0.1,0.1), pch=16, cex=0.3, axes=F, xlab="", ylab="", asp=1)
lines(minshapep_uoj[c(3,1,4,2),],lwd=4)
dev.off()

# 4.5 Two-block PLS and Lower Pharyngeal Jaw Morphology
################################################################################

# 4.5.1 PLS without accounting for phylogenetic non-independence

# Make matrix
pls_diet_lpj <- as.matrix(t(assay(vsdata_spp)))
# Replace colnames by diet phylum
colnames(pls_diet_lpj) <- with(Tax_table_strict_2, phylum[match(colnames(pls_diet_lpj),rownames(Tax_table_strict_2))])
pls_diet_lpj <- as.matrix(pls_diet_lpj)
# Load data Ronco et al
load("sp_means_LPJ.Rdata")
pls_lpj <- sp_means_LPJ
# Remove Astbur
pls_lpj <- pls_lpj[ ,, dimnames(pls_lpj)[[3]] %in%  rownames(pls_diet_lpj)]
# Order in the same way
pls_diet_lpj <- pls_diet_lpj[match(dimnames(pls_lpj)[[3]], rownames(pls_diet_lpj)),]
all(dimnames(pls_lpj)[[3]] == rownames(pls_diet_lpj))
# Run 2-block PLS
f_lpj = two.b.pls(pls_diet_lpj, pls_lpj) 
# Plot PLS
P_lpj = plot(f_lpj)
# Summary
summary(f_lpj)
# Check loading of the projections
f_lpj$left.pls.vectors
f_lpj$right.pls.vectors
# Get projected scores
f1_lpj= data.frame( "XScores"=f_lpj$ XScores[,1], "YScores"=f_lpj$ YScores[,1])
f1_lpj$Tribe = with(all_metadata_2, Tribe[match(rownames(f1_lpj),SpeciesID)])
# R^2 for the PLS fit
summary(lm(f1_lpj$YScores ~ f1_lpj$XScores))
# pGLS
f1_lpj$SpeciesID <- rownames(f1_lpj)
compset_lpj = comparative.data(pruned.tree, f1_lpj, SpeciesID)
print(summary(pgls(YScores ~ XScores , compset_lpj, lambda= 'ML')))

# 4.5.2 Plot PLS without accounting for phylogenetic non-independence

# Visualize shape at minimum and maximum PLS scores
minx_lpj <- min(P_lpj$plot_args$y)
maxx_lpj <- max(P_lpj$plot_args$y)
preds_lpj <- shape.predictor(P_lpj$A2, 
                             x = P_lpj$plot.args$y,
                             min = minx_lpj, max = maxx_lpj)
plotRefToTarget(mshape(P_lpj$A2), preds_lpj$min)
plotRefToTarget(mshape(P_lpj$A2), preds_lpj$max)
minshape_lpj = preds_lpj$min
maxshape_lpj = preds_lpj$max
# Plot PLS
plot_pls_lpj <-ggplot(data = f1_lpj, aes(XScores, YScores,colour=Tribe,size=3)) +
  geom_point() +
  theme_bw() +  
  scale_color_manual(values=tribe_colours) +
  ylab("lower pharyngeal jaw shape projection") +
  xlab("food source projection") +
  stat_smooth(method = "lm", col = "darkgray",se = FALSE,size=1, linetype=2) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text.x = element_text(color = "black", size = 11, face = "plain"),
        axis.text.y = element_text(color = "black", size = 11),  
        axis.title.x = element_text(color = "black", size = 15, hjust = .5, vjust = 0),
        axis.title.y = element_text(color = "black", size = 15, vjust = 2),
        legend.title = element_text(size=15), legend.text = element_text(size=12),
        legend.position="none",
        plot.margin = unit(c(0.2,0.2,0.2,1), "cm"))
# Outlines for shapes
outline1_lpj=c(1,12,11,13,10,14,9,15,8,3,23,24,25,5,6,7,22,16,21,17,20,18,19,1)
outline2_lpj=c(1,19,18,20,17,21,16,22,7,6,5,25,24,23,3,39,40,41,29,30,31,38,32,37,33,36,34,35,1)
outline3_lpj=c(1,12,11,13,10,14,9,15,8,3)
# Plot max shape
pdf(file = "Plot_LPJ_NoPhy_max_10012025.pdf",
    width = 7,
    height = 7) 
par(mar = c(0,0,0,0))
plot(maxshape_lpj[c(1:27),1], -1* maxshape_lpj[c(1:27),3], xlim=c(-0.25,0.75), ylim=c(-0.5,0.5), pch=16, type="n", axes=F, xlab="", ylab="",asp=1)
polygon(maxshape_lpj[c(outline1_lpj),1]  +0.5 , -1* maxshape_lpj[c(outline1_lpj),3]  +0.1 ,lwd=4)
segments(maxshape_lpj[15,1]  +0.5 ,-1* maxshape_lpj[15,3]  +0.1 , maxshape_lpj[26,1]  +0.5 , -1* maxshape_lpj[26,3]  +0.1 ,lwd=4)
segments(maxshape_lpj[7,1]  +0.5 ,-1* maxshape_lpj[7,3]  +0.1 , maxshape_lpj[26,1]  +0.5 ,  -1* maxshape_lpj[26,3]  +0.1 ,lwd=4)
polygon(maxshape_lpj[c(4,2,27),1]  +0.5 , -1* maxshape_lpj[c(4,2,27),3]  +0.1, col= rgb(0.51,0.5,0.5,alpha=0.2) ,lwd=4)
lines( y= maxshape_lpj[outline3_lpj,1],-1* maxshape_lpj[outline3_lpj,2], col="grey",lwd=4)
segments( y0= maxshape_lpj[15,1],-1* maxshape_lpj[15,2],  y1= maxshape_lpj[26,1],  -1* maxshape_lpj[26,2], col="grey",lwd=4)
segments( y0= maxshape_lpj[7,1],-1* maxshape_lpj[7,2],  y1=  maxshape_lpj[26,1],  -1* maxshape_lpj[26,2], col="grey",lwd=4)
segments( y0= maxshape_lpj[42,1],-1* maxshape_lpj[42,2],   y1= maxshape_lpj[15,1],  -1* maxshape_lpj[15,2], col="grey",lwd=4)
segments( y0= maxshape_lpj[42,1],-1* maxshape_lpj[42,2],   y1= maxshape_lpj[31,1],  -1* maxshape_lpj[31,2], col="grey",lwd=4)
polygon( y= maxshape_lpj[c(outline2_lpj),1], -1* maxshape_lpj[c(outline2_lpj),2],lwd=4)
polygon( y= maxshape_lpj[c(4,2,28,27),1], -1* maxshape_lpj[c(4,2,28,27),2], col= rgb(0.51,0.5,0.5,alpha=0.2),lwd=4)
polygon(maxshape_lpj[c(outline2_lpj),2]  +0.5 , -1* maxshape_lpj[c(outline2_lpj),3]  -0.1 ,lwd=4)
lines(maxshape_lpj[outline3_lpj,2]  +0.5 ,-1* maxshape_lpj[outline3_lpj,3]  -0.1 , col="black",lwd=4)
segments(maxshape_lpj[15,2]  +0.5 ,-1* maxshape_lpj[15,3]  -0.1 , maxshape_lpj[26,2]  +0.5 ,  -1* maxshape_lpj[26,3]  -0.1 , col="black",lwd=4)
segments(maxshape_lpj[7,2]  +0.5 ,-1* maxshape_lpj[7,3]  -0.1 , maxshape_lpj[26,2]  +0.5 ,  -1* maxshape_lpj[26,3]  -0.1 , col="black",lwd=4)
segments(maxshape_lpj[42,2]  +0.5 ,-1* maxshape_lpj[42,3]  -0.1 , maxshape_lpj[15,2]  +0.5 ,  -1* maxshape_lpj[15,3]  -0.1 , col="black",lwd=4)
segments(maxshape_lpj[42,2]  +0.5 ,-1* maxshape_lpj[42,3]  -0.1 , maxshape_lpj[31,2]  +0.5 ,  -1* maxshape_lpj[31,3]  -0.1 , col="black",lwd=4)
polygon(maxshape_lpj[c(4,2,28,27),2]+0.5, -1* maxshape_lpj[c(4,2,28,27),3] -0.1, col= rgb(0.51,0.5,0.5,alpha=0.2),lwd=4)
dev.off()
# Plot min shape
pdf(file = "Plot_LPJ_NoPhy_min_10012025.pdf",
    width = 7,
    height = 7) 
par(mar = c(0,0,0,0))
plot(minshape_lpj[c(1:27),1], -1* minshape_lpj[c(1:27),3], xlim=c(-0.25,0.75), ylim=c(-0.5,0.5), pch=16, col="white", axes=F, xlab="", ylab="", asp=1)
polygon(minshape_lpj[c(outline1_lpj),1]  +0.5 , -1* minshape_lpj[c(outline1_lpj),3]  +0.1 ,lwd=4)
segments(minshape_lpj[15,1]  +0.5 ,-1* minshape_lpj[15,3]  +0.1 , minshape_lpj[26,1]  +0.5 , -1* minshape_lpj[26,3]  +0.1 ,lwd=4)
segments(minshape_lpj[7,1]  +0.5 ,-1* minshape_lpj[7,3]  +0.1 , minshape_lpj[26,1]  +0.5 ,  -1* minshape_lpj[26,3]  +0.1 ,lwd=4)
polygon(minshape_lpj[c(4,2,27),1]  +0.5 , -1* minshape_lpj[c(4,2,27),3]  +0.1, col= rgb(0.51,0.5,0.5,alpha=0.2) ,lwd=4)
lines( y= minshape_lpj[outline3_lpj,1],-1* minshape_lpj[outline3_lpj,2], col="grey",lwd=4)
segments( y0= minshape_lpj[15,1],-1* minshape_lpj[15,2],  y1= minshape_lpj[26,1],  -1* minshape_lpj[26,2], col="grey",lwd=4)
segments( y0= minshape_lpj[7,1],-1* minshape_lpj[7,2],  y1=  minshape_lpj[26,1],  -1* minshape_lpj[26,2], col="grey",lwd=4)
segments( y0= minshape_lpj[42,1],-1* minshape_lpj[42,2],   y1= minshape_lpj[15,1],  -1* minshape_lpj[15,2], col="grey",lwd=4)
segments( y0= minshape_lpj[42,1],-1* minshape_lpj[42,2],   y1= minshape_lpj[31,1],  -1* minshape_lpj[31,2], col="grey",lwd=4)
polygon( y= minshape_lpj[c(outline2_lpj),1], -1* minshape_lpj[c(outline2_lpj),2],lwd=4)
polygon( y= minshape_lpj[c(4,2,28,27),1], -1* minshape_lpj[c(4,2,28,27),2], col= rgb(0.51,0.5,0.5,alpha=0.2),lwd=4)
polygon(minshape_lpj[c(outline2_lpj),2]  +0.5 , -1* minshape_lpj[c(outline2_lpj),3]  -0.1 ,lwd=4)
lines(minshape_lpj[outline3_lpj,2]  +0.5 ,-1* minshape_lpj[outline3_lpj,3]  -0.1 , col="black",lwd=4)
segments(minshape_lpj[15,2]  +0.5 ,-1* minshape_lpj[15,3]  -0.1 , minshape_lpj[26,2]  +0.5 ,  -1* minshape_lpj[26,3]  -0.1 , col="black",lwd=4)
segments(minshape_lpj[7,2]  +0.5 ,-1* minshape_lpj[7,3]  -0.1 , minshape_lpj[26,2]  +0.5 ,  -1* minshape_lpj[26,3]  -0.1 , col="black",lwd=4)
segments(minshape_lpj[42,2]  +0.5 ,-1* minshape_lpj[42,3]  -0.1 , minshape_lpj[15,2]  +0.5 ,  -1* minshape_lpj[15,3]  -0.1 , col="black",lwd=4)
segments(minshape_lpj[42,2]  +0.5 ,-1* minshape_lpj[42,3]  -0.1 , minshape_lpj[31,2]  +0.5 ,  -1* minshape_lpj[31,3]  -0.1 , col="black",lwd=4)
polygon(minshape_lpj[c(4,2,28,27),2]+0.5, -1* minshape_lpj[c(4,2,28,27),3] -0.1, col= rgb(0.51,0.5,0.5,alpha=0.2),lwd=4)
dev.off()

# 4.5.3 PLS with accounting for phylogenetic non-independence

# Order by phylogenetric tree
pls_diet_ord_lpj <- pls_diet_lpj[ match(pruned.tree_noOrAs$tip.label, rownames(pls_diet_lpj)),]
pls_bdy_ord_lpj <- pls_lpj[,, match(pruned.tree_noOrAs$tip.label, dimnames(pls_lpj)[[3]])]
# Check
all(pruned.tree_noOrAs$tip.label == rownames(pls_diet_ord_lpj)) 
all(pruned.tree_noOrAs$tip.label == dimnames(pls_bdy_ord_lpj)[[3]])
# Perform PLS
fp_lpj = phylo.integration(  A = pls_diet_ord_lpj,  A2 = pls_bdy_ord_lpj,  phy = pruned.tree_noOrAs  )
# Plot
Pp_lpj = plot(fp_lpj)
# Summary
summary(fp_lpj)
# Check loading of the SI projections
fp_lpj$left.pls.vectors
fp_lpj$right.pls.vectors
# Get projected scores
f1p_lpj= data.frame( "XScores"=fp_lpj$ XScores[,1], "YScores"=fp_lpj$ YScores[,1])
f1p_lpj$Tribe = with(all_metadata_2, Tribe[match(rownames(f1p_lpj),SpeciesID)])
# R^2 for the PLS fit
summary(lm(f1p_lpj$YScores ~ f1p_lpj$XScores))

# 4.5.4 Plot PLS with accounting for phylogenetic non-independence

# Visualize shape at minimum and maximum PLS scores.
minxp_lpj <- min(Pp_lpj$plot_args$y)
maxxp_lpj <- max(Pp_lpj$plot_args$y)
predsp_lpj <- shape.predictor(Pp_lpj$A2, 
                              x = Pp_lpj$plot.args$y,
                              min = minxp_lpj, max = maxxp_lpj)
plotRefToTarget(mshape(Pp_lpj$A2), predsp_lpj$min)
plotRefToTarget(mshape(Pp_lpj$A2), predsp_lpj$max)
minshapep_lpj = predsp_lpj$min
maxshapep_lpj = predsp_lpj$max
# Plot PLS
plot_pls_phy_lpj <- ggplot(data = f1p_lpj, aes(XScores, YScores,colour=Tribe,size=3)) +
  geom_point() +
  theme_bw() +  
  scale_color_manual(values=tribe_colours) +
  ylab("lower pharyngeal jaw shape projection") +
  xlab("food source projection") +
  stat_smooth(method = "lm", col = "darkgray",se = FALSE,size=1, linetype=2) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text.x = element_text(color = "black", size = 11, face = "plain"),
        axis.text.y = element_text(color = "black", size = 11),  
        axis.title.x = element_text(color = "black", size = 15, hjust = .5, vjust = 0),
        axis.title.y = element_text(color = "black", size = 15, vjust = 2),
        legend.title = element_text(size=15), legend.text = element_text(size=12),
        legend.position="none",
        plot.margin = unit(c(0.2,0.2,0.2,1), "cm"))
# Outlines to draw the shapes
outline1p_lpj=c(1,12,11,13,10,14,9,15,8,3,23,24,25,5,6,7,22,16,21,17,20,18,19,1)
outline2p_lpj=c(1,19,18,20,17,21,16,22,7,6,5,25,24,23,3,39,40,41,29,30,31,38,32,37,33,36,34,35,1)
outline3p_lpj=c(1,12,11,13,10,14,9,15,8,3)
# Plot max shape
pdf(file = "Plot_LPJ_Phy_max_13012025.pdf",
    width = 7,
    height = 7) 
par(mar = c(0,0,0,0))
plot(maxshapep_lpj[c(1:27),1], -1* maxshapep_lpj[c(1:27),3], xlim=c(-0.25,0.75), ylim=c(-0.5,0.5), pch=16, type="n", axes=F, xlab="", ylab="",asp=1)
polygon(maxshapep_lpj[c(outline1p_lpj),1]  +0.5 , -1* maxshapep_lpj[c(outline1p_lpj),3]  +0.1 ,lwd=4)
segments(maxshapep_lpj[15,1]  +0.5 ,-1* maxshapep_lpj[15,3]  +0.1 , maxshapep_lpj[26,1]  +0.5 , -1* maxshapep_lpj[26,3]  +0.1,lwd=4 )
segments(maxshapep_lpj[7,1]  +0.5 ,-1* maxshapep_lpj[7,3]  +0.1 , maxshapep_lpj[26,1]  +0.5 ,  -1* maxshapep_lpj[26,3]  +0.1,lwd=4 )
polygon(maxshapep_lpj[c(4,2,27),1]  +0.5 , -1* maxshapep_lpj[c(4,2,27),3]  +0.1, col= rgb(0.51,0.5,0.5,alpha=0.2) ,lwd=4)
lines( y= maxshapep_lpj[outline3p_lpj,1],-1* maxshapep_lpj[outline3p_lpj,2], col="grey",lwd=4)
segments( y0= maxshapep_lpj[15,1],-1* maxshapep_lpj[15,2],  y1= maxshapep_lpj[26,1],  -1* maxshapep_lpj[26,2], col="grey",lwd=4)
segments( y0= maxshapep_lpj[7,1],-1* maxshapep_lpj[7,2],  y1=  maxshapep_lpj[26,1],  -1* maxshapep_lpj[26,2], col="grey",lwd=4)
segments( y0= maxshapep_lpj[42,1],-1* maxshapep_lpj[42,2],   y1= maxshapep_lpj[15,1],  -1* maxshapep_lpj[15,2], col="grey",lwd=4)
segments( y0= maxshapep_lpj[42,1],-1* maxshapep_lpj[42,2],   y1= maxshapep_lpj[31,1],  -1* maxshapep_lpj[31,2], col="grey",lwd=4)
polygon( y= maxshapep_lpj[c(outline2p_lpj),1], -1* maxshapep_lpj[c(outline2p_lpj),2],lwd=4)
polygon( y= maxshapep_lpj[c(4,2,28,27),1], -1* maxshapep_lpj[c(4,2,28,27),2], col= rgb(0.51,0.5,0.5,alpha=0.2),lwd=4)
polygon(maxshapep_lpj[c(outline2p_lpj),2]  +0.5 , -1* maxshapep_lpj[c(outline2p_lpj),3]  -0.1,lwd=4 )
lines(maxshapep_lpj[outline3p_lpj,2]  +0.5 ,-1* maxshapep_lpj[outline3p_lpj,3]  -0.1 , col="black",lwd=4)
segments(maxshapep_lpj[15,2]  +0.5 ,-1* maxshapep_lpj[15,3]  -0.1 , maxshapep_lpj[26,2]  +0.5 ,  -1* maxshapep_lpj[26,3]  -0.1 , col="black",lwd=4)
segments(maxshapep_lpj[7,2]  +0.5 ,-1* maxshapep_lpj[7,3]  -0.1 , maxshapep_lpj[26,2]  +0.5 ,  -1* maxshapep_lpj[26,3]  -0.1 , col="black",lwd=4)
segments(maxshapep_lpj[42,2]  +0.5 ,-1* maxshapep_lpj[42,3]  -0.1 , maxshapep_lpj[15,2]  +0.5 ,  -1* maxshapep_lpj[15,3]  -0.1 , col="black",lwd=4)
segments(maxshapep_lpj[42,2]  +0.5 ,-1* maxshapep_lpj[42,3]  -0.1 , maxshapep_lpj[31,2]  +0.5 ,  -1* maxshapep_lpj[31,3]  -0.1 , col="black",lwd=4)
polygon(maxshapep_lpj[c(4,2,28,27),2]+0.5, -1* maxshapep_lpj[c(4,2,28,27),3] -0.1, col= rgb(0.51,0.5,0.5,alpha=0.2),lwd=4)
dev.off()
# Plot min shape
pdf(file = "Plot_LPJ_Phy_min_13012025.pdf",
    width = 7,
    height = 7) 
par(mar = c(0,0,0,0))
plot(minshapep_lpj[c(1:27),1], -1* minshapep_lpj[c(1:27),3], xlim=c(-0.25,0.75), ylim=c(-0.5,0.5), pch=16, col="white", axes=F, xlab="", ylab="",asp=1)
polygon(minshapep_lpj[c(outline1p_lpj),1]  +0.5 , -1* minshapep_lpj[c(outline1p_lpj),3]  +0.1 ,lwd=4)
segments(minshapep_lpj[15,1]  +0.5 ,-1* minshapep_lpj[15,3]  +0.1 , minshapep_lpj[26,1]  +0.5 , -1* minshapep_lpj[26,3]  +0.1,lwd=4 )
segments(minshapep_lpj[7,1]  +0.5 ,-1* minshapep_lpj[7,3]  +0.1 , minshapep_lpj[26,1]  +0.5 ,  -1* minshapep_lpj[26,3]  +0.1,lwd=4 )
polygon(minshapep_lpj[c(4,2,27),1]  +0.5 , -1* minshapep_lpj[c(4,2,27),3]  +0.1, col= rgb(0.51,0.5,0.5,alpha=0.2),lwd=4 )
lines( y= minshapep_lpj[outline3p_lpj,1],-1* minshapep_lpj[outline3p_lpj,2], col="grey",lwd=4)
segments( y0= minshapep_lpj[15,1],-1* minshapep_lpj[15,2],  y1= minshapep_lpj[26,1],  -1* minshapep_lpj[26,2], col="grey",lwd=4)
segments( y0= minshapep_lpj[7,1],-1* minshapep_lpj[7,2],  y1=  minshapep_lpj[26,1],  -1* minshapep_lpj[26,2], col="grey",lwd=4)
segments( y0= minshapep_lpj[42,1],-1* minshapep_lpj[42,2],   y1= minshapep_lpj[15,1],  -1* minshapep_lpj[15,2], col="grey",lwd=4)
segments( y0= minshapep_lpj[42,1],-1* minshapep_lpj[42,2],   y1= minshapep_lpj[31,1],  -1* minshapep_lpj[31,2], col="grey",lwd=4)
polygon( y= minshapep_lpj[c(outline2p_lpj),1], -1* minshapep_lpj[c(outline2p_lpj),2],lwd=4)
polygon( y= minshapep_lpj[c(4,2,28,27),1], -1* minshapep_lpj[c(4,2,28,27),2], col= rgb(0.51,0.5,0.5,alpha=0.2),lwd=4)
polygon(minshapep_lpj[c(outline2p_lpj),2]  +0.5 , -1* minshapep_lpj[c(outline2p_lpj),3]  -0.1 ,lwd=4)
lines(minshapep_lpj[outline3p_lpj,2]  +0.5 ,-1* minshapep_lpj[outline3p_lpj,3]  -0.1 , col="black",lwd=4)
segments(minshapep_lpj[15,2]  +0.5 ,-1* minshapep_lpj[15,3]  -0.1 , minshapep_lpj[26,2]  +0.5 ,  -1* minshapep_lpj[26,3]  -0.1 , col="black",lwd=4)
segments(minshapep_lpj[7,2]  +0.5 ,-1* minshapep_lpj[7,3]  -0.1 , minshapep_lpj[26,2]  +0.5 ,  -1* minshapep_lpj[26,3]  -0.1 , col="black",lwd=4)
segments(minshapep_lpj[42,2]  +0.5 ,-1* minshapep_lpj[42,3]  -0.1 , minshapep_lpj[15,2]  +0.5 ,  -1* minshapep_lpj[15,3]  -0.1 , col="black",lwd=4)
segments(minshapep_lpj[42,2]  +0.5 ,-1* minshapep_lpj[42,3]  -0.1 , minshapep_lpj[31,2]  +0.5 ,  -1* minshapep_lpj[31,3]  -0.1 , col="black",lwd=4)
polygon(minshapep_lpj[c(4,2,28,27),2]+0.5, -1* minshapep_lpj[c(4,2,28,27),3] -0.1, col= rgb(0.51,0.5,0.5,alpha=0.2),lwd=4)
dev.off()

# 4.6 Plot PLS figures
################################################################################

# 4.6.1 Main Figure 4

# Make plots without phylogeny
pdf(file = "Plot_PLS_NoPhyl_13012025.pdf",
    width = 20,
    height = 10) 
grid.arrange(plot_pls_iso,plot_pls_bdy,plot_pls_uoj,plot_pls_lpj,nrow=1)
dev.off()

# 4.6.2 Plot with phylogeny

# Make plots with phylogeny
pdf(file = "Plot_PLS_WithPhyl_13012025.pdf",
    width = 20,
    height = 10) 
grid.arrange(plot_pls_phy_iso,plot_pls_phy_bdy, plot_pls_phy_uoj, plot_pls_phy_lpj,nrow=1)
dev.off()


