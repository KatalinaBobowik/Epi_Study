# created by KSB, 09.03.20
# ICA analysis for variables associated with expression profiles
# using the tutorial as per: http://127.0.0.1:23118/library/MineICA/doc/MineICA.pdf

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

# Set paths:
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/dataPreprocessing/"
refDir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/"
outputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/"

# 1. load in expression data, log transform, then mean center data ----------------------------------------------------------------

# load expression data and OTU file
load(paste0(inputdir, "indoRNA.read_counts.TMM.filtered.Rda"))
OTUs=read.table(paste0(refDir,"scaledOTUs.txt"))

# transform OTU file to make it the same configuration as the other variables, then make sure sample names are 
# the same before appending to DGE list
OTUs=t(OTUs)
samplenamesOTU <- as.character(rownames(OTUs))
samplenamesOTU <- gsub("\\.","-", samplenamesOTU)
samplenamesOTU <- gsub("Batch3","", samplenamesOTU)
samplenamesOTU <- gsub("Batch2","", samplenamesOTU)
samplenamesOriginal <- as.character(rownames(y$samples))
samplenamesOriginal <- sapply(strsplit(samplenamesOriginal, "[_.]"), `[`, 1)

# check to see if samplenames are in the same order beofre appending
identical(samplenamesOriginal,samplenamesOTU[match(samplenamesOriginal,samplenamesOTU)])
# TRUE 
OTUs=OTUs[match(samplenamesOriginal,samplenamesOTU),]
OTUs2=OTUs
OTUs2[which(OTUs2>0)]=T
colnames(OTUs2)=paste(colnames(OTUs2), "Binary",sep="_")

# append to DGE list
for (name in colnames(OTUs)){
  y$samples[[paste0(name)]]<- OTUs[,name]
}

for (name in colnames(OTUs2)){
  y$samples[[paste0(name)]]<- as.factor(OTUs2[,name])
}

# get factor variable of the highest OTU within each sample
maxContributor = sapply(1:123, function(x) names(which.max(OTUs[x,])))
y$samples$maxContributor = as.factor(maxContributor)

# rename DGE list
indoSampleSet = y
indoSampleSet$samples$Age[which(is.na(y$samples$Age) == T)]=45
lcpm <- cpm(indoSampleSet, log=TRUE)

# get batch-corrected data
design <- model.matrix(~0 + indoSampleSet$samples$Island)
colnames(design)=gsub("Island", "", colnames(design))
colnames(design)=gsub("indoSampleSet", "", colnames(design))
colnames(design)=gsub("samples", "", colnames(design))
colnames(design)=gsub("[\\$$]", "", colnames(design))
colnames(design)=gsub("West Papua", "Mappi", colnames(design))
batch.corrected.lcpm <- removeBatchEffect(lcpm, batch=indoSampleSet$samples$batch, covariates = cbind(indoSampleSet$samples$Age, indoSampleSet$samples$RIN, indoSampleSet$samples$CD8T, indoSampleSet$samples$CD4T, indoSampleSet$samples$NK, indoSampleSet$samples$Bcell, indoSampleSet$samples$Mono, indoSampleSet$samples$Gran), design=design)

# 2. Run ICA using JADE ---------------------------------------------------------

# ICA setup

## Features are mean-centered before ICA computation
# indoSampleSet <- t(apply(indoSampleSet,1,scale,scale=FALSE))
indoSampleSet <- t(apply(batch.corrected.lcpm,1,scale,scale=FALSE))
colnames(indoSampleSet) <- colnames(y)

# indoSampleSet = indoSampleSet[,grep("MPI",colnames(indoSampleSet))]

## run ICA-JADE with 5 components and 10,000 iterations
resJade <- runICA(X=indoSampleSet, nbComp=5, method = "JADE", maxit=10000) 

## build parameters. Selected cutoff applied to the absolute feature/gene projection values to be consider as contributors is 3
params <- buildMineICAParams(resPath="ICA/", selCutoff=3, pvalCutoff=0.05)

## define the reference samples if any, here no normal sample is available
refSamplesindoSampleSet <- character(0)
rownames(resJade$A) = colnames(indoSampleSet)

# set up annotation
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
typeIDindoSampleSet <- c(geneID_annotation="SYMBOL", geneID_biomart="hgnc_symbol", featureID_biomart="ensembl_gene_id")

# Run ICA using JADE 
resBuild <- buildIcaSet(params=params, A=data.frame(resJade$A), S=data.frame(resJade$S), dat=indoSampleSet, pData=y$samples, typeID= typeIDindoSampleSet, refSamples=refSamplesindoSampleSet, mart=mart)
# resBuild <- buildIcaSet(params=params, A=data.frame(resJade$A), S=data.frame(resJade$S), dat=indoSampleSet, pData=y$samples[grep("MPI",colnames(y)),], typeID= typeIDindoSampleSet, refSamples=refSamplesindoSampleSet, mart=mart)

# extract data from anIcaSetinstance
icaSetIndo <- resBuild$icaSet
params <- resBuild$params
annot <- pData(icaSetIndo)

# 3. Get genes driving variation in the first 5 PCs ---------------------------------------------------------

# When applying ICA decomposition to genomicdata, for example here gene expression data, the distribution of the gene 
# projections on the ICs is expected to be super-Gaussian: a large portion of genes follows a (super-)Gaussian centered 
# at zero and a small portion belongs to an outgrowth located on the right and/or on the left of thedistribution.  
# In order to select the elements belonging to this outgrowth, I use a threshold ofe 3 standard deviations from the mean.  
# These can be refered to as the “contributing genes”.

contrib <- selectContrib(icaSetIndo, cutoff=3, level="genes")
## Show the first contributing genes of the first and third components
sort(abs(contrib[[1]]),decreasing=TRUE)[1:10]
sort(abs(contrib[[3]]),decreasing=TRUE)[1:10]
# Here is the histogram of the projection values for the first 
pdf(paste0(outputdir,"hist.pdf"))
hist(S(icaSetIndo)[,1], breaks=50, main="Distribution of feature projection on the first component", xlab="projection values") 
abline(v=c(3,-3), col="red", lty=2)
dev.off()

# One can also extract data from a specific component
# e.g., extract sample contributions and gene projections of the second component
comp2 <- getComp(icaSetIndo, level="genes", ind=2)
## access the sample contributions 
comp2$contrib[1:5]
## access the gene projections
comp2$proj[1:5]

# Get the scaled projections for each gene

# restrict the phenotype data to the variables of interest
varLabels(icaSetIndo)
keepVar <- c("Island","maxContributor")
resW <- writeProjByComp(icaSet=icaSetIndo, params=params, mart=mart, level='genes', selCutoffWrite=2.5)
head(resW$listAnnotComp[[1]])
## show the number of components a gene contributes to
head(resW$nbOccInComp)

# Which variables are driving variability in each PC? --------------------------------------------

# Association with sample variables 

# test  whether  the  groups  of  samples  formed  by  the  qualitative  variables  
# are differently distributed on the components in terms of contribution value 
# If you would like to plot densities rather than boxplots, use'typePlot=density'
OTUQual=c()
for (sample in colnames(OTUs2)){
  if(length(levels(y$samples[,sample])) > 1){
    OTUQual=c(OTUQual,sample)
  }
}

resQual <- qualVarAnalysis(params=params, icaSet=icaSetIndo, 
                           keepVar=c("Island",OTUQual),
                           adjustBy="none", doPlot=T,
                           path="qualVarAnalysis2/", typePlot="boxplot",filename="qualVar")

### Quantitative variables
## Compute pearson correlations between variable 'age' and the sample contributions
## on all components.
## We are interested in correlations exceeding 0.3 in absolute value, and plots will only be drawn
## for correlations exceeding this threshold.
resQuant <- quantVarAnalysis(params=params, icaSet=icaSetIndo, keepVar=c(colnames(OTUs),"Age"), 
                             typeCor="pearson", cutoffOn="pval",
                             cutoff=0.05, adjustBy="none",  
                             path="quantVarAnalysis/", filename="quantVar", doPlot=T)


# get Pearson correlation
resQuant$cor

# make a heatmap of the correlations
pdf(paste0(outputdir,"heatmapQuantVarAnalysisCorrelations.pdf"))
Heatmap(resQuant$cor, cluster_columns=F, col=colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red")))
dev.off()

# Study the bimodality of sample contributions in matrixA -----------------------

# The distribution of the sample contributions on a component is often bimodal
# Let's plot this out

# Selection of samples associated with a component
resmix <- plotAllMix(A=A(icaSetIndo), nbMix=2, nbBreaks=50)
## ## plot the positions of the samples on the second component according to their ER status
## ## (in a file "er.pdf") 
plotPosAnnotInComp(icaSet=icaSetIndo, params=params, keepVar=keepVar, funClus="Mclust")


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
                            keepVar=colnames(OTUs), 
                            resClus=clus1$clus, funClus="Mclust", adjustBy="none", 
                            doPlot=TRUE, path="clus2var/", filename="resChitests-Mcluscomp1")

clus2var <- clusVarAnalysis(icaSet=icaSetIndo, params=params, 
                            keepVar=colnames(OTUs), 
                            resClus=clus1$clus, funClus="Mclust", adjustBy="none", 
                            doPlot=TRUE, path="clus2var/", filename="resChitests-Mcluscomp1")


## p-values are also contained in the ouput of the function:
clus2var


##
## Run the analysis of the ICA decomposition
# runAn(params=params, icaSet=icaSetIndo, writeGenesByComp = TRUE,keepVar=keepVar, dbGOstats = "KEGG")
keepVar <- c("Island",colnames(OTUs2)[14])

#  Plot heatmaps of the contributing elements
resH <- plot_heatmapsOnSel(icaSet = icaSetIndo, selCutoff = 3, level = "genes", keepVar = keepVar,doSamplesDendro = TRUE, doGenesDendro = TRUE, heatmapCol = maPalette(low = "blue", high = "red", mid = "yellow", k=44),file = "heatmapWithDendro", annot2col=annot2col(params))
# heatmap where genes and samples are ordered by contribution values
resH <- plot_heatmapsOnSel(icaSet = icaSetIndo, selCutoff = 3, level = "genes",keepVar = keepVar,doSamplesDendro = FALSE, doGenesDendro = FALSE, heatmapCol = maPalette(low = "blue", high = "red", mid = "yellow", k=44),file = "heatmapWithoutDendro", annot2col=annot2col(params))



# ------

SMBvsMTW=read.table("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/DE_Island/LM_allCovarPlusBlood/topTable_SMBvsMTW.txt")
SMBvsMTW=rownames(SMBvsMTW)[which(abs(SMBvsMTW$logFC) >= 1)]
SMBvsMPI=read.table("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/DE_Island/LM_allCovarPlusBlood/topTable_SMBvsMPI.txt")
SMBvsMPI=rownames(SMBvsMPI)[which(abs(SMBvsMPI$logFC) >= 1)]
MTWvsMPI=read.table("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/DE_Island/LM_allCovarPlusBlood/topTable_MTWvsMPI.txt")
MTWvsMPI=rownames(MTWvsMPI)[which(abs(MTWvsMPI$logFC) >= 1)]

#The data were centered and reduced to one unit of variance prior performing the PCA analysis.
indoSampleSet <- t(apply(batch.corrected.lcpm,1,scale,scale=FALSE))
#indoSampleSet <- t(apply(lcpm,1,scale,scale=FALSE))

#apply(doubs$env,2,scale,scale=FALSE)

y$samples$Age[which(is.na(y$samples$Age) == T)]=45

test=list()
for (i in 1:100){
  env=y$samples[sample(nrow(y$samples),50),c(10,13,16:21,27:42)]
  genes=indoSampleSet[unique(c(SMBvsMTW,SMBvsMPI)),]
  genes=t(genes)
  colnames(genes)=y[colnames(genes),]$genes$SYMBOL
  colnames(genes)=make.unique(colnames(genes))
  rownames(genes) = colnames(batch.corrected.lcpm)
  genes=genes[rownames(env),]
  dudi1 <- dudi.pca(env, scale = T, scan = F, nf = 2)
  dudi2 <- dudi.pca(genes, scale = T, scan = F, nf = 2)
  coin1 <- coinertia(dudi1,dudi2, scan = F, nf = 2)
  test[[i]]=coin1$tab
}
new=t(sapply(1:24, function(x) rowMeans(test %>% map(x) %>% invoke(cbind, .))))
rownames(new)=colnames(env)
Heatmap(new, col=colorRamp2(c(-1, 0, 1), c("blue", "white", "red")))

# or 
env=y$samples[,colnames(OTUs)]
genes=indoSampleSet[unique(c(SMBvsMTW,SMBvsMPI)),]
genes=t(genes)
colnames(genes)=y[colnames(genes),]$genes$SYMBOL
colnames(genes)=make.unique(colnames(genes))
rownames(genes) = colnames(batch.corrected.lcpm)
dudi1 <- dudi.pca(env, scale = T, scan = F, nf = 2)
dudi2 <- dudi.pca(genes, scale = F, scan = F, nf = 2)
coin1 <- coinertia(dudi1,dudi2, scan = F, nf = 2)
pdf(paste0(outputdir,"correlationHeatmap.pdf"), width=15)
Heatmap(t(coin1$tab), col=colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red")))
dev.off()

# Just for MPI
# or 
env=y$samples[grep("MPI",rownames(y$samples)),colnames(OTUs)]
genes=indoSampleSet[unique(c(MTWvsMPI,SMBvsMPI)),]
#genes=indoSampleSet

genes=t(genes)
colnames(genes)=y[colnames(genes),]$genes$SYMBOL
colnames(genes)=make.unique(colnames(genes))
rownames(genes) = colnames(batch.corrected.lcpm)
genes=genes[grep("MPI",rownames(genes)),]
dudi1 <- dudi.pca(env, scale = T, scan = F, nf = 2)
dudi2 <- dudi.pca(genes, scale = T, scan = F, nf = 2)
coin1 <- coinertia(dudi1,dudi2, scan = F, nf = 2)
pdf(paste0(outputdir,"correlationHeatmap.pdf"), width=15)
Heatmap(t(coin1$tab), col=colorRamp2(c(-1, 0, 1), c("blue", "white", "red")))
dev.off()
names=names(which(abs(colSums(coin1$tab))>0))
Heatmap(t(coin1$tab[,which(abs(colSums(coin1$tab))>0)]),col=colorRamp2(c(-1, 0, 1), c("blue", "white", "red")))


# correlation between samples and Island
a=hetcor(y$samples$Island, y$samples[,names])
b=data.frame(a$correlations[,1],a$tests[,1])
b=b[order(b$a.correlations...1.),]
newOTUs=rownames(b)[which(abs(b[,1])>=0.1)]
newOTUs=newOTUs[-which(newOTUs %in% "y$samples$Island")]

# more testing - PCA
# PCA plotting function
plot.pca <- function(dataToPca, speciesCol, namesPch, sampleNames){
    pca <- prcomp(t(dataToPca), scale=T, center=T)
    pca.var <- pca$sdev^2/sum(pca$sdev^2)
    for (i in 1:9){
        pca_axis1=i
        pca_axis2=i+1
        plot(pca$x[,pca_axis1], pca$x[,pca_axis2], col=speciesCol, pch=namesPch, cex=2, xlab=paste0("PC", pca_axis1, " (", round(pca.var[pca_axis1]*100, digits=2), "% of variance)"), ylab=paste0("PC", pca_axis2, " (", round(pca.var[pca_axis2]*100, digits=2), "% of variance)", sep=""), main=name)
        # points(pca$x[,pca_axis1][which(allreplicated==T)], pca$x[,pca_axis2][which(allreplicated==T)], col="black", pch=8, cex=2)
        # text(pca$x[,pca_axis1][which(allreplicated==T)], pca$x[,pca_axis2][which(allreplicated==T)], labels=samplenames[which(allreplicated==T)], pos=3)
        #legend(legend=unique(sampleNames), pch=16, x="bottomright", col=speciesCol[unique(as.numeric(get(name)))], cex=0.6, title=name, border=F, bty="n")
        #legend(legend=unique(sampleNames), pch=16, x="bottomright", col=unique(speciesCol), cex=0.6, title=name, border=F, bty="n")
        #legend(legend=unique(sampleNames), pch=16, x="bottomright", col=village.col[unique(as.numeric(get(name)))], cex=0.6, title=name, border=F, bty="n")
        #legend(legend=unique(as.numeric(y$samples$batch)), "topright", pch=unique(as.numeric(y$samples$batch)) + 14, title="Batch", cex=0.6, border=F, bty="n")
    }

    return(pca)
}

# PCA association function
pc.assoc <- function(pca.data){
    all.pcs <- data.frame()
    for (i in 1:ncol(pca.data$x)){
        all.assoc <- vector()
        for (j in 1:ncol(all.covars.df)){
            test.assoc <- anova(lm(pca.data$x[,i] ~ all.covars.df[,j]))[1,5]
            all.assoc <- c(all.assoc, test.assoc)
        }
        single.pc <- c(i, all.assoc)
        all.pcs <- rbind(all.pcs, single.pc)
    }
    names(all.pcs) <- c("PC", colnames(all.covars.df))

    print ("Here are the relationships between PCs and some possible covariates")
    print (all.pcs)
    return (all.pcs)
}

numericVariables=newOTUs
for (name in numericVariables){
 assign(name, y$samples[[paste0(name)]])
}

# Prepare covariate matrix
all.covars.df <- y$samples[,numericVariables]

for (name in numericVariables){
    initial = .bincode(get(name), breaks=seq(min(get(name), na.rm=T), max(get(name), na.rm=T), len = 80),include.lowest = TRUE)
    bloodCol <- colorRampPalette(c("blue", "red"))(79)[initial]
    pdf(paste0(outputdir,"pcaresults_",name,".pdf"))
    pcaresults <- plot.pca(dataToPca=batch.corrected.lcpm, speciesCol=bloodCol,namesPch=as.numeric(y$samples$Island) + 14,sampleNames=get(name))
    legend(legend=c("High","Low"), pch=16, x="bottomright", col=c(bloodCol[which.max(get(name))], bloodCol[which.min(get(name))]), cex=0.6, title=name, border=F, bty="n")
    legend(legend=unique(as.numeric(y$samples$batch)), "topright", pch=unique(as.numeric(y$samples$batch)) + 14, title="Batch", cex=0.6, border=F, bty="n")
    dev.off()
}

# Get PCA associations
all.pcs <- pc.assoc(pcaresults)
all.pcs$Variance <- pcaresults$sdev^2/sum(pcaresults$sdev^2)
write.table(all.pcs, file=paste0(outputdir,"pca_covariates_blood_RNASeqDeconCell.txt"), col.names=T, row.names=F, quote=F, sep="\t")

