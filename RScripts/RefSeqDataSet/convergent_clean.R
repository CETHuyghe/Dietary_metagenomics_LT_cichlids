#### Fabrizia Ronco
#### January 2025

#### set we

require(geiger)  
require(reshape2) 
require(hexbin) 



### get data
d = read.csv("~/Documents/Science/Projects/Charlotte/PCA.csv")
head(d)
dim(d)

## make species means
d2 = aggregate(d[,2:5], by=list(d$Tribe, d$ SpeciesID), mean)

## remove Astbur and Oretan
sum(d2$Group.2 %in% c("Astbur", "Oretan"))
d2 = d2[!d2$Group.2 %in% c("Astbur", "Oretan"),]

## calculate SE 
d2_se = aggregate(d[,2:5], by=list(d$Tribe, d$ SpeciesID), function(x) {sd(x) / sqrt(length(x)) } )
d2_se = d2_se[!d2$Group.2 %in% c("Astbur", "Oretan"),]

###. get tree and prune
phy= read.nexus("~/Documents/Science/Projects/CichlidX/final_stats_for_MS/01_Data/b1.tre")
phy2 = drop.tip(phy, phy$tip.label[!phy$tip.label %in% d2$Group.2])

all(phy2$tip.label %in% d2$Group.2)
all(d2$Group.2 %in% phy2$tip.label)

#####. get rates per PC axis 
PC1 = d2$PC1; names(PC1) = d2$Group.2
PC2 = d2$PC2; names(PC2) = d2$Group.2
PC3 = d2$PC3; names(PC3) = d2$Group.2
PC4 = d2$PC4; names(PC4) = d2$Group.2

PC1_se = d2_se$PC1; names(PC1_se) = d2_se$Group.2
PC2_se = d2_se$PC2; names(PC2_se) = d2_se$Group.2
PC3_se = d2_se$PC3; names(PC3_se) = d2_se$Group.2
PC4_se = d2_se$PC4; names(PC4_se) = d2_se$Group.2

BM_PC1 = fitContinuous(phy2, PC1, SE=PC1_se, model="BM")
BM_PC2 = fitContinuous(phy2, PC2, SE=PC2_se, model="BM")
BM_PC3 = fitContinuous(phy2, PC3, SE=PC3_se, model="BM")
BM_PC4 = fitContinuous(phy2, PC4, SE=PC4_se, model="BM")


PC1_z0 = BM_PC1$opt$z0
PC1_sig2 = BM_PC1$opt$sigsq 
PC2_z0 = BM_PC2$opt$z0
PC2_sig2 = BM_PC2$opt$sigsq 
PC3_z0 = BM_PC3$opt$z0
PC3_sig2 = BM_PC3$opt$sigsq 
PC4_z0 = BM_PC4$opt$z0
PC4_sig2 = BM_PC4$opt$sigsq 


#############################################
### make nsim sims for eac set of parameters
nsim=10000

sims_PC1 = fastBM(phy2, a=PC1_z0, sig2= PC1_sig2, nsim = nsim)
sims_PC2 = fastBM(phy2, a=PC2_z0, sig2= PC2_sig2, nsim = nsim)
sims_PC3 = fastBM(phy2, a=PC3_z0, sig2= PC3_sig2, nsim = nsim)
sims_PC4 = fastBM(phy2, a=PC4_z0, sig2= PC4_sig2, nsim = nsim)


#### get distances in the empirical data 
dist_2D_empirical = sqrt ( (dist(PC1)^2) + (dist(PC2)^2) + (dist(PC3)^2) + (dist(PC4)^2) )

dist_2D_empirical = as.matrix(dist_2D_empirical)
diag(dist_2D_empirical) = NA
dist_2D_empirical[upper.tri(dist_2D_empirical)]= NA

dist_2D_empirical_melt = melt(dist_2D_empirical, na.rm=T)
dist_2D_empirical_melt$cont = paste(dist_2D_empirical_melt$Var1, dist_2D_empirical_melt$Var2, sep="_")

#### get genetic distances, and attach to the other distances
gendist = read.table("~/Documents/Science/Projects/CichlidX/Genomes/genetic_dists/f5.dists_based_on_callability.095.spc.txt")
gendist = as.matrix(gendist)

## subset anbd order to match empricical trait data 
order_gendist = match(row.names(dist_2D_empirical), row.names(gendist))
gendist = gendist[order_gendist, order_gendist]
all(row.names(dist_2D_empirical) == row.names(gendist))

diag(gendist)=NA
gendist[upper.tri(gendist)]=NA
gendist_melt = melt(gendist, na.rm=T)
gendist_melt$cont = paste(gendist_melt $Var1, gendist_melt$Var2, sep="_")
all(gendist_melt$cont == dist_2D_empirical_melt$cont )

dist_2D_empirical_melt$gen_dist = gendist_melt$value

### normalise everything to 0:1
normalize = function(x){(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))}
dist_2D_empirical_melt$gen_dist = normalize(dist_2D_empirical_melt$gen_dist)
dist_2D_empirical_melt$value = normalize(dist_2D_empirical_melt$value)

#plot(dist_2D_empirical_melt$gen_dist,dist_2D_empirical_melt$value )
#text(dist_2D_empirical_melt$gen_dist,dist_2D_empirical_melt$value,  dist_2D_empirical_melt$cont ,cex=0.5 )

#write.csv(dist_2D_empirical_melt, "~/Documents/Science/Projects/Charlotte/PC1_4_dist_norm.csv", row.names=F, quote=F)

#### do the same for the simulations
#head(sims_PC1)
## reorder
sims_PC1 = sims_PC1[match(d2$Group.2, row.names(sims_PC1)),]
sims_PC2 = sims_PC2[match(d2$Group.2, row.names(sims_PC2)),]
sims_PC3 = sims_PC3[match(d2$Group.2, row.names(sims_PC3)),]
sims_PC4 = sims_PC4[match(d2$Group.2, row.names(sims_PC4)),]


simdists = as.data.frame(matrix(NA, dim(dist_2D_empirical_melt)[1], nsim))
for ( i in 1:nsim ) {
	#tmp = sqrt((dist(sims_PC1[,i]))^2 + (dist(sims_PC2[,i]))^2)
	#tmp = sqrt( (dist(sims_PC1[,i]))^2 + (dist(sims_PC2[,i]))^2 + (dist(sims_PC3[,i]))^2 )
	tmp = sqrt( (dist(sims_PC1[,i]))^2 + (dist(sims_PC2[,i]))^2 + (dist(sims_PC3[,i]))^2 + (dist(sims_PC4[,i]))^2 )
	tmp = as.matrix(tmp)
	diag(tmp)=NA
	tmp[upper.tri(tmp)]=NA
	tmp2 = melt(tmp, na.rm=T)
	tmp2$cont = paste(tmp2$Var1, tmp2$Var2, sep="_")
	if (all(tmp2$cont == dist_2D_empirical_melt$cont)) {
		simdists[,i] = tmp2$value
	}
}


#head(simdists)
simdists = apply (simdists, 2, normalize)

#i=1
#plot(dist_2D_empirical_melt$gen_dist, simdists[,i] )
#i=i+1


################ binning in hexagons
Nbins=10

## get range of all pheno data (simulations and empirical) to set bounds 
range_all = range(c(dist_2D_empirical_melt$value, simdists))

## empirical data 
empirical_bin1 = hexbin(x = dist_2D_empirical_melt$gen_dist, y = dist_2D_empirical_melt$value, xbins = Nbins, ybnds = range_all) 
empirical_bin =  cbind(empirical_bin1@cell, empirical_bin1@count)

#plot(empirical_bin1)

## simulations 
sim_bins = list()
for ( i in 1:nsim) {
	tmp_bin = hexbin(x = dist_2D_empirical_melt$gen_dist, y = simdists[,i], xbins = Nbins, ybnds = range_all) 
	sim_bins[[i]] = cbind(tmp_bin@cell, tmp_bin@count)
}

#head(simdists)
#i=1
#plot(hexbin(x = dist_2D_empirical_melt$gen_dist, y = simdists[,i], xbins = Nbins, ybnds = range_all))
#i=i+1

################
## add the 0 bins

bin_num_max_sim = max(unlist(lapply(sim_bins, function(x) max( x[,1]))))
bin_num_max_emp = max(empirical_bin[,1])
bin_num_max  = max(bin_num_max_emp, bin_num_max_sim )

## empirical data
all_bins= c(1:bin_num_max)

emptybins = all_bins[! all_bins %in% empirical_bin[,1] ]

empirical_bin = rbind(empirical_bin, cbind(emptybins, rep(0, length(emptybins))))
empirical_bin = empirical_bin[order(empirical_bin[,1]),]


## simulations 
out_all_bins = data.frame("cell"=c(1:bin_num_max))
for ( i in 1:nsim) {
	emptybins = all_bins[! all_bins %in% sim_bins[[i]][,1] ]
	sim_bins[[i]] = rbind(sim_bins[[i]], cbind(emptybins, rep(0, length(emptybins))))
	sim_bins[[i]] = sim_bins[[i]][order(sim_bins[[i]][,1]),]
	out_all_bins[,i+1] = sim_bins[[i]][,2]
}

#sim_bins

##### identifiy bins that are differnet

frac = rep(NA, dim(out_all_bins)[1])
for (k in 1:dim(out_all_bins)[1]){
	tmp_vec= unlist(out_all_bins[k,-1])
	## per bin, the proportion of simulations that have a lower denisty 
	frac[k] = mean(tmp_vec < empirical_bin[k,2])
}


#frac > 0.95
#sig < 0.05


frac2=frac
frac2[frac > 0.95] = 1
frac2[frac <= 0.95] = 0

frac3=frac
frac3[frac > 0.9] = 1
frac3[frac <= 0.9] = 0

frac4=frac
frac4[frac > 0.85] = 1
frac4[frac <= 0.85] = 0

empirical_bin2 = empirical_bin1
#str(empirical_bin2)
empirical_bin2@xcm=0
empirical_bin2@ycm=0
empirical_bin2@ncells=length(frac2)
empirical_bin2@cell=as.integer(empirical_bin[,1])
empirical_bin2@count=empirical_bin[,2]


pdf("~/Documents/Science/Projects/Charlotte/hex_dist_PC1-4_bin10_sim10000_95_90_85.pdf")
P = plot(empirical_bin2, , colorcut=seq(0, 1, length.out=max(empirical_bin2@count)))
pushHexport(P$plot.vp)
grid.hexagons(empirical_bin2, use.count=F, cell.at = frac4, style="colorscale",colramp = function(n) {colorRampPalette(c("gold"))(n)}, mincnt=1, maxcnt=1 , border=F )
grid.hexagons(empirical_bin2, use.count=F, cell.at = frac3, style="colorscale",colramp = function(n) {colorRampPalette(c("orange"))(n)}, mincnt=1, maxcnt=1 , border=F )
grid.hexagons(empirical_bin2, use.count=F, cell.at = frac2, style="colorscale",colramp = function(n) {colorRampPalette(c("darkred"))(n)}, mincnt=1, maxcnt=1 , border=F )
dev.off()





