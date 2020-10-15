# created by KSB, 09.03.20
# Question: Is pathogen load or species type responsible for why we see a strong immune signature 
# detected in peripheral blood within Korowai

library(Biobase)
library(plyr)
library(ggplot2)
library(foreach)
library(xtable)
library(biomaRt)
library(GOstats)
library(cluster)
library(marray)
library(mclust)
library(RColorBrewer)
library(igraph)
library(Rgraphviz)
library(graph)
library(colorspace)
library(annotate)
library(scales)
library(gtools)
library(JADE)
library(edgeR)
library(limma)
library(MineICA)
library(ade4)
library(ComplexHeatmap)
library(circlize)
library(polycor)
library(ggpubr)
library(reshape2)

# Set paths:
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/dataPreprocessing/"
refDir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/"
outputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/"

# 1. Data preprocessing ----------------------------------------------------------------

# load expression data and OTU file
load(paste0(inputdir, "indoRNA.read_counts.TMM.filtered.Rda"))
OTUs=read.table(paste0(refDir,"scaledOTUs.txt"))
# OTUs=read.table(paste0(refDir,"OTUtable.txt"))

# transform OTU file to make it the same configuration as the expression list
OTUs=t(OTUs)
# make 'NA' values into 0's
OTUs[is.na(OTUs)] <- 0

# check to see if samplenames in the OTU file are in the same order as the expression list before appending
samplenamesOTU <- as.character(rownames(OTUs))
samplenamesOTU <- gsub("\\.","-", samplenamesOTU)
samplenamesOTU <- gsub("Batch1","", samplenamesOTU)
samplenamesOTU <- gsub("Batch3","", samplenamesOTU)
samplenamesOTU <- gsub("Batch2","", samplenamesOTU)
samplenamesOriginal <- as.character(rownames(y$samples))
samplenamesOriginal <- sapply(strsplit(samplenamesOriginal, "[_.]"), `[`, 1)
identical(samplenamesOriginal,samplenamesOTU[match(samplenamesOriginal,samplenamesOTU)])
# TRUE 

# match OTU data oreder to expression list oreder 
OTUs=OTUs[match(samplenamesOriginal,samplenamesOTU),]

# append to DGE list
for (name in colnames(OTUs)){
  y$samples[[paste0(name)]]<- OTUs[,name]
}

# Get batch-corrected data -----------------------------------

# age is a variable we use in the model and NA values are not allowed. We'll replace NA with the mean age value
y$samples$Age[which(is.na(y$samples$Age) == T)]=45

# # get batch-corrected data
lcpm <- cpm(y, log=TRUE)
design <- model.matrix(~0 + y$samples$Island)
# rename columns to be "Mentawai", "Sumba", and "Mappi"
colnames(design)=gsub("Island", "", colnames(design))
colnames(design)=gsub("y", "", colnames(design))
colnames(design)=gsub("samples", "", colnames(design))
colnames(design)=gsub("[\\$$]", "", colnames(design))
colnames(design)=gsub("West Papua", "Mappi", colnames(design))
# regress out variables influencing expression
batch.corrected.lcpm <- removeBatchEffect(lcpm, batch=y$samples$batch, covariates = cbind(y$samples$Age, y$samples$RIN, y$samples$CD8T, y$samples$CD4T, y$samples$NK, y$samples$Bcell, y$samples$Mono, y$samples$Gran), design=design)

# 2. Run ICA using JADE ---------------------------------------------------------

# in the MineICA tutorial, they use a microarray dataset, so well use our batch-corrected log-CPM data to make
# it similar. As noted in the MineICA tutorial, features are also mean-centered before ICA computation.

indoSampleSet <- t(apply(batch.corrected.lcpm,1,scale,scale=FALSE))
colnames(indoSampleSet) <- colnames(y)

# Run ICA-JADE with 5 components and 10,000 iterations
# We restrict the data to the 10,000 genes with the highest IQR
# Here, we're using Jade because it's faster than the other option, fastICA. 
# FastICA relies on random initializations and the estimated components may vary between iterations.
resJade <- runICA(X=indoSampleSet, nbComp=5, method = "JADE", maxit=10000) 

# Create a MineICAParams object -------------------------------------

## build parameters. Selected cutoff applied to the absolute feature/gene projection values to be consider as contributors is 3
# refer to line 139 for selectCutoff (3 std deviations from mean)
params <- buildMineICAParams(resPath=paste0(outputdir,"ICA/"), selCutoff=3, pvalCutoff=0.05)

# Create an IcaSet instance -------------------------------------------

# set up annotation
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
typeIDindoSampleSet <- c(geneID_annotation="SYMBOL", geneID_biomart="ensembl_gene_id", featureID_biomart="ensembl_gene_id")

## define the reference samples if any, here no normal sample is available
refSamplesindoSampleSet <- character(0)

# Run ICA using JADE 
rownames(resJade$A) = colnames(indoSampleSet)
resBuild <- buildIcaSet(params=params, A=data.frame(resJade$A), S=data.frame(resJade$S), dat=indoSampleSet, pData=y$samples, typeID= typeIDindoSampleSet, refSamples=refSamplesindoSampleSet, mart=mart)

# extract data from anIcaSetinstance
icaSetIndo <- resBuild$icaSet
params <- resBuild$params
# You can access the phenotype data using pData
annot <- pData(icaSetIndo)

# Exploratpry analysis of IcaSet -------------------------------------------------------

# When applying ICA decomposition to genomicdata, the distribution of the gene 
# projections on the ICs is expected to be super-Gaussian: a large portion of genes 
# follows a (super-)Gaussian curve centered at zero and a small portion belongs to an 
# outgrowth located on the right and/or on the left of the distribution.  
# In order to select the elements belonging to this outgrowth, I use a threshold of 3 standard deviations from the mean.  
# These can be refered to as the “contributing genes”.

pdf(paste0(outputdir,"hist.pdf"))
hist(S(icaSetIndo)[,1], breaks=50, main="Distribution of feature projection on the first component", xlab="projection values") 
abline(v=c(3,-3), col="red", lty=2)
dev.off()

# Get genes driving variation in the first 5 PCs 
contrib <- selectContrib(icaSetIndo, cutoff=3, level="genes")
## Show the first contributing genes of the first and third components
sort(abs(contrib[[1]]),decreasing=TRUE)[1:10]

###################################################
#### Run the analysis of the ICA decomposition ####
###################################################

# we're going to run the ICA anlaysis using individual functions
# there's an option to run everything at once using the 'runAn' function
# bit running everything individually gives more control

# Write description of contributing genes -----------------------

# The function writeProjByComp investigates which genes 
# are contributing to a component given the selected threshold 
# This function creates an HTML file per component containing the description 
# of the contributing genes, and a file containing the projection values of each gene across all components
resW <- writeProjByComp(icaSet=icaSetIndo, params=params, mart=mart, level='genes', selCutoffWrite=2.5)

# Get DE genes for each island
# these are all DE genes at no log fold change threshold and an adjusted p-value of 0.05
SMBvsMTW=read.table("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/DE_Island/LM_allCovarPlusBlood/topTable_SMBvsMTW.txt")
SMBvsMTW_genes=rownames(SMBvsMTW)[which(abs(SMBvsMTW$logFC) >= 1)]
SMBvsMPI=read.table("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/DE_Island/LM_allCovarPlusBlood/topTable_SMBvsMPI.txt")
SMBvsMPI_genes=rownames(SMBvsMPI)[which(abs(SMBvsMPI$logFC) >= 1)]
MTWvsMPI=read.table("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/DE_Island/LM_allCovarPlusBlood/topTable_MTWvsMPI.txt")
MTWvsMPI_genes=rownames(MTWvsMPI)[which(abs(MTWvsMPI$logFC) >= 1)]

# get barplots of proportion of contributing genes in each PCs
contributingGenes=matrix(nrow=3,ncol=5)
rownames(contributingGenes)=c("SMBvsMPI","MTWvsMPI","SMBvsMTW")
colnames(contributingGenes)=c("PC1","PC2","PC3","PC4","PC5")
for (i in 1:5){
  contributingGenes[1,i]=length(contrib[[i]][which(names(contrib[[i]]) %in% SMBvsMPI_genes)])
  #contributingGenes[1,2]=length(SMBvsMPI)
  contributingGenes[2,i]=length(contrib[[i]][which(names(contrib[[i]]) %in% MTWvsMPI_genes)])
  #contributingGenes[2,2]=length(MTWvsMPI)
  contributingGenes[3,i]=length(contrib[[i]][which(names(contrib[[i]]) %in% SMBvsMTW_genes)])
}
contributingGenes=as.data.frame(contributingGenes)
contributingGenes$type="contribGenes"
contributingGenes$Island=rownames(contributingGenes)

deGenes=do.call("cbind", replicate(5, c(length(SMBvsMPI_genes),length(MTWvsMPI_genes),length(SMBvsMTW_genes)), simplify = FALSE))
deGenes=as.data.frame(deGenes)
rownames(deGenes)=c("SMBvsMPI","MTWvsMPI","SMBvsMTW")
colnames(deGenes)=c("PC1","PC2","PC3","PC4","PC5")
deGenes$type="deGenes"
deGenes$Island=rownames(deGenes)
allGenes=rbind(contributingGenes,deGenes)
meltedGenes=melt(allGenes)
# reorder
meltedGenes$type=factor(meltedGenes$type, levels=c("deGenes","contribGenes"))

pdf(paste0(outputdir,"numberOfContributingGenesAndDEGenes.pdf"), height=5,width=14)
ggplot(a, aes(x = Island, y = value, fill = type)) + 
  geom_bar(stat = 'identity', position = 'fill') + facet_grid(~ variable) + theme_bw() + scale_fill_manual(values = c("#4477AA","#EE6677")) 
dev.off()

# make table with resW contributions to each PCA and DE information 

# TODO: finishg this section #########
# get pc1
resW[[1]][1]
#################

# Plot heatmaps of the contributing elements ----------------------

# get which pathogens contribute most to each sample
maxContributor = sapply(1:123, function(x) names(which.max(OTUs[x,])))
sort(table(maxContributor))
# only works for quantitative data
keepVar <- c("Eukaryota_Plasmodiidae","Viruses_Flaviviridae","Bacteria_Enterobacteriaceae","Bacteria_Pseudomonadaceae")

#  Plot heatmaps. The heatmap is ranked by sample contributions
resH <- plot_heatmapsOnSel(icaSet = icaSetIndo, selCutoff = 3, level = "genes", keepVar = keepVar,doSamplesDendro = TRUE, doGenesDendro = TRUE, heatmapCol = maPalette(low = "blue", high = "red", mid = "yellow", k=44),file = "heatmapWithDendro", annot2col=annot2col(params), path = paste0(outputdir,"ICA/Heatmaps/"))

# Association with sample variables ---------------------------------------

## qualitative variable analysis

# test  whether  the  groups  of  samples  formed  by  qalitative  variables  
# are differently distributed on the components in terms of contribution value 

resQual <- qualVarAnalysis(params=params, icaSet=icaSetIndo, 
                           keepVar=c("Island","Mappi"),
                           adjustBy="none", doPlot=T,
                           path="qualVarAnalysis/", typePlot="boxplot",filename="qualVar")

## Quantitative variable analysis

# Compute pearson correlations between variable numeric variables and the sample contributions
# on all components. We are interested in correlations exceeding 0.3 in absolute value, and plots will only be drawn
# for correlations exceeding this threshold.

# note on witness gene: When applying ICA to gene expression data,
# each component is typically triggered by a group of genes co-expressed on a subset of samples.
# These genes responsible for the existence of the component will typically have high projections, we
# call them the contributing genes.

resQuant <- quantVarAnalysis(params=params, icaSet=icaSetIndo, keepVar=c(colnames(OTUs)), 
                             typeCor="pearson", cutoffOn="pval",
                             cutoff=0.05, adjustBy="none",  
                             path="quantVarAnalysis/", filename="quantVar", doPlot=T)

# get Pearson correlation and reset any 'NA' values to 0
resQuant$cor[is.na(resQuant$cor)] <- 0
correlation=resQuant$cor

# make a heatmap of the correlations
pdf(paste0(outputdir,"heatmapQuantVarAnalysisCorrelations_MPI.pdf"))
Heatmap(correlation, cluster_columns=F, col=colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red")))
dev.off()

# Study the bimodality of sample contributions in matrixA -----------------------

# The distribution of the sample contributions on a component is often bimodal
# Let's plot this out

# Selection of samples associated with a component
resmix <- plotAllMix(A=A(icaSetIndo), nbMix=2, nbBreaks=50)
## ## plot the positions of the samples on the second component according to their ER status
## ## (in a file "er.pdf") 
plotPosAnnotInComp(icaSet=icaSetIndo, params=params, keepVar=c("Island","Mappi"), funClus="Mclust")


# Cluster samples --------------------------------

# cluster the samples using the mixing matrix A on teh first and second component
clus1 <- clusterSamplesByComp(params=params, icaSet=icaSetIndo[,,,1:2],
                              funClus="Mclust", clusterOn="A", nbClus=2, filename="comp1Mclust")

## The obtained clusters are written in the file "comp1Mclus.txt" of the result path.
clus1$clus[[2]][1:5]

# or perform several clusterings, using different algorithms 
clus2 <- clusterSamplesByComp_multiple(params=params, icaSet=icaSetIndo[,,1:2], 
                                     funClus="kmeans", clusterOn=c("A","S"), level="features", 
                                     nbClus=2, filename="comparKmeans")

## Access Rand index
clus2$comparClus

# get association with the qualitative variables
# perform  the  chi-square  tests  of independence to study the association between the clustering obtained on each component and the qualitative variables. 
clus2var <- clusVarAnalysis(icaSet=icaSetIndo[,,1:2], params=params, 
                            keepVar=c("Island","Mappi"), 
                            resClus=clus1$clus, funClus="Mclust", adjustBy="none", 
                            doPlot=TRUE, path="clus2var/", filename="resChitests-Mcluscomp1")


