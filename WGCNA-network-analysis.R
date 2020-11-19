# ————————————————————————————
# WGCNA-network-analysis.R
# AUTHOR: Amal Katrib
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
# OBJECTIVE:
#  Perform ene co-expression network analysis using Steve Horvath's
#  Weighted Gene Co-Expression Network Analysis (WGCNA)
#
# ———————————————————————————————
rm( list = ls (all = TRUE))
options(stringsAsFactors = F)
allowWGCNAThreads()

# load packages
library(WGCNA)
library(cluster)
library(dplyr)
library(tibble)
library(flashClust)
library(matrixStats)
library(ggplot2)
library(DESeq2)
library(XLConnect)
library(gplots)
library(reshape2)
library(dendextend)
library(pheatmap)
library(RColorBrewer)
library(grid)
library(stringr)


# ------------------------------------------------------
#  MANUAL INPUT
# ------------------------------------------------------
# set directory
dir = "x"

# ------------------------------------------------------
#  DATA INPUT
# ..............
setwd(dir)
# ------------------------------------------------------

# (a) load gene symbol conversion
geneSymbol = read.csv("GeneInfo.csv")
# (b) load log-transformed gene expression data
counts = read.csv("geneCounts.tsv", sep="\t")
# (c) load sample info / covariates
cov = read.csv("sampleInfo.txt")

# ------------------------------------------------------
#  DATA PRE-PROCESSING (optional)
# ------------------------------------------------------

# filter out low reads (median > 0)
counts = counts[rowMedians(as.matrix(counts[,-1]))>0,]
genes = counts$GeneSymbol
counts = counts[,-1]
# remove NA entries
y = cbind(genes,counts)
y = na.omit(y)
genes = y$genes
counts = y[,-1]
# plot data distribution
plot(density(as.matrix(counts[1,])))

# ------------------------------------------------------
#  NETWORK CONSTRUCTION
# ------------------------------------------------------

# sample network based on squared Euclidean distance
A=adjacency(counts,type="distance")
# calculate whole network connectivity
k=as.numeric(apply(A,2,sum))-1
# standardized connectivity
Z.k=scale(k)
# designate samples as outliers if their Z.k value is below the threshold
thresholdZ.k=-2.5 # often -2.5
# the color vector indicates outlyingness (red)
outlierColor=ifelse(Z.k<thresholdZ.k,"red","black")
# generate cluster tree using flahsClust or hclust
sampleTree = flashClust(as.dist(1-A), method = "average")
datColors=data.frame(outliers=outlierColor)
# change colors and labels (optional)
datColors[datColors=="black"] ="aliceblue"
datColors[datColors=="red"] ="dodgerblue4"
# check sample labeling
sampleTree$labels

# plot the sample dendrogram and show outliers
plotDendroAndColors(sampleTree, groupLabels=names(datColors),
cex.colorLabels = 0.8, cex.dendroLabels = 0.9,
abHeight = 0.25, abCol = "blue",
colors=datColors, colorHeight = 0.05,
main="Squared Euclidean Distance Clustering")

# remove outlier samples (optional)
remove.samples= Z.k<thresholdZ.k | is.na(Z.k)
counts = counts[,!remove.samples]
# recompute sample network among the remaining samples
A=adjacency(counts,type="distance")
# recompute Z.k values of outlyingness
k=as.numeric(apply(A,2,sum))-1
Z.k=scale(k)

# ------------------------------------------------------
#  NETWORK CONSTRUCTION
# ..............
# SCALE-FREE TOPOLOGY / SOFT-THRESHOLDING:
# to construct a weighted network (soft thresholding with the power
# adjacency matrix), we consider a vector of potential thresholds &
# investigate soft thesholding with the power adjacency function
# ------------------------------------------------------

# default = signed network. For unsigned, just add "networkType = "unsigned"
powers = c(seq(1,10,by=1),seq(12,20,by=2))
sft = pickSoftThreshold(t(counts), powerVector=powers, networkType = "signed hybrid")[[2]]

# plot scale free fit R^2 vs different soft threshold beta
par(mfrow = c(1,2))
plot(sft[,1], -sign(sft[,3])*sft[,2],
xlab=" Soft Threshold (power)", ylab="Scale Free Topology Model Fit, R^2",type="n")
text(sft[,1], -sign(sft[,3])*sft[,2], labels=powers,cex=0.7,col="red")
abline(h=0.75,col="red")
plot(sft[,1], sft[,5],xlab="Soft Threshold (power)",
ylab="Mean Connectivity", type="n")
text(sft[,1], sft[,5], labels=powers, cex=0.8,col="red")
dev.off()

# MANUAL INPUT: specify 1st power value above cutoff line for power adjacency function
beta = 10
k.data = softConnectivity(t(counts),type="signed hybrid",power=beta)-1
# create scale free topology plot
png(paste0("ScaleFreePlot_signed_",Sys.Date(),".png"),9,5,res=300,units="in")
scaleFreePlot(k.data, main=paste("Scale Free Toplogy Criterion, Power=",beta), truncated=F);
dev.off()

# ------------------------------------------------------
#  NETWORK MODULE DETECTION
# ..............
# DYNAMIC TREE CUT:
# use blockwise for automatic network construction, which calculates topological
# overlap measure & performs dynamic tree cut to identify modules
# ------------------------------------------------------

# construct network!
mergingThresh = 0.25
net = blockwiseModules(t(counts),networkType="signed hybrid",
corType="bicor",pearsonFallback = "individual",
power=beta,maxPOutliers = 0.1,
maxBlockSize=5000,minModuleSize=30,
minkmetoStay=0.7,mergeCutHeight=mergingThresh,
numericLabels=T,numericlabels=T,reassignThreshold=0,
pamRespectsDendro=FALSE,deepSplit=2,verbose=5,saveTOMs = F)

# get module labels
moduleLabels=net$colors
# convert labels to colors for plotting
moduleColors = labels2colors(moduleLabels)
# dataframe of module eigengenes
MEs=net$MEs
rownames(MEs) = colnames(counts)
geneTree = net$dendrograms[[1]]
# calculate absolute correlation between gene expression & module eigengene >= |0.7|
# this is referred to as kme - gene membership
datkme=signedkme(t(counts), MEs)

# ------------------------------------------------------
#  NETWORK VISUALIZATION
# ------------------------------------------------------

# plot the dendrogram and the module colors underneath
blocknumber=1
datColors = data.frame(moduleColors)[net$blockGenes[[blocknumber]],]
datColors = as.data.frame(datColors)
plotDendroAndColors(net$dendrograms[[blocknumber]],
  groupLabels=c("Modules"), colors=datColors[,1], cex.colorLabels = 1.1,
  dendroLabels=FALSE, hang=0.03,addGuide=TRUE,guideHang=0.05)



# ------------------------------------------------------
#  NETWORK ANALYSIS METRICS
# ..............
# TABLE WITH MODULE INFO:
# Module # ; Module Color ; # Genes in Module ;
# # Hub Genes in module ; % Hub Genes in Module
# ------------------------------------------------------

# adjust module data to keep it orderly (per module #)
kme= datkme[,order(as.integer(gsub("kme","",colnames(datkme))))]
moduleLabels = as.data.frame(moduleLabels)
rownames(moduleLabels) = rownames(kme)
moduleLabels$moduleLabels = paste0("M",moduleLabels$moduleLabels)
# create table with module labels / colors. Keep in order of increasing mod #
table = data.frame(mod = unique(moduleLabels$moduleLabels), modColor = unique(moduleColors))
table = table[order(as.integer(gsub("M","",table$mod))),]
# specify # genes in each module
modGenes = vector(mode="list",length = nrow(table))
for (i in 1:nrow(table)) {
    modGenes[[i]] = sum(table$modColor[i] == moduleColors) }
table = cbind(table, modGenes = as.integer(modGenes))
# specify % of module genes relative to all genes
modGenesRatio = as.numeric(modGenes)/nrow(counts)
table = cbind(table, modSizeRatio = round(100*modGenesRatio,2))
# use gene symbols if it's easier than IDs
genes = geneSymbol$symbol[match(rownames(moduleLabels),geneSymbol$id)]
# identify module hub genes, with |kme| >= 0.9 ... very non-elegant, but works for now
hubGenes = vector(mode="list", length=ncol(kme))
NoHubGenes = vector(mode="list",length=ncol(kme))
for (i in 1:ncol(kme)) {
    hubGenes[[i]] = genes[which(abs(kme[,i]) >= 0.9)]
    NoHubGenes[[i]] = length(hubGenes[[i]]) }
names(hubGenes) = gsub("kme","M",colnames(kme))
# specify module hub genes + their ratio relative to module size
table = cbind(table,hubGenes = as.integer(NoHubGenes))
table = table %>% mutate(hubGenesRatio = round(100*hubGenes/modGenes,2))

# ------------------------------------------------------
#  MODULE-FEATURE CORRELATION
# ------------------------------------------------------

# define # genes and samples
nGenes = ncol(t(dataset))
nSamples = nrow(t(dataset))
# recalculate MEs with color labels
MEs0 = moduleEigengenes(t(dataset),t(moduleColors))$eigengenes
MEs1 = orderMEs(MEs0)
rownames(MEs1) = rownames(MEs)
MEs = MEs1
rm(MEs1)

## —— [OPTIONAL]: change feature entries into numerical

# correlate using appropriate method
modTraitCor = stats::cor(MEs, cov, method="spearman", use = "complete.obs")
modTraitCor = stats::cor(MEs, cov, method="pearson", use = "complete.obs")
modTraitCor = naRemove(bicor(MEs, cov, use="p", robustY = F, maxPOutliers = 0.1))
# get corresponding pVal
modTraitP = corPvalueStudent(modTraitCor, nSamples)

# Generate a graphical representation to read in the table.
# Color code each association by the correlation value and display correlations as well as p-values
textMatrix = paste(signif(modTraitCor, 2), "\n(",signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

# Time to plot!!!
# Display the correlation values within a heatmap plot
