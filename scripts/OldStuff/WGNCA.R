# created by KSB, 09.03.20
# WGCNA analysis for variables associated with expression profiles
# using the tutorial as per: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html

### R code from vignette source 'MineICA.Rnw'

###################################################
### code chunk number 1: morelib
###################################################
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


###################################################
### code chunk number 2: lib
###################################################
library(MineICA)


###################################################
### code chunk number 3: helpIcaSet (eval = FALSE)
###################################################
## help(IcaSet)


###################################################
### code chunk number 4: loadMainz
###################################################
## load Mainz expression data and sample annotations.
library(breastCancerMAINZ)
data(mainz)
show(mainz)
## we restrict the data to the 10,000 probe sets with the highest IQR
mainz <- selectFeatures_IQR(mainz,10000)


###################################################
### code chunk number 5: runJade
###################################################
library(JADE)
## Features are mean-centered before ICA computation
exprs(mainz) <- t(apply(exprs(mainz),1,scale,scale=FALSE))
colnames(exprs(mainz)) <- sampleNames(mainz)
## run ICA-JADE
resJade <- runICA(X=exprs(mainz), nbComp=5, method = "JADE", maxit=10000) 


###################################################
### code chunk number 6: fastica (eval = FALSE)
###################################################
## library(fastICA)
## ## Random initializations are used for each iteration of FastICA
## ## Estimates are clustered using hierarchical clustering with average linkage
## res <- clusterFastICARuns(X=exprs(mainz), nbComp=5, alg.type="deflation", nbIt=10, 
##                           funClus="hclust", method="average")


###################################################
### code chunk number 7: buildParams
###################################################
## build params
params <- buildMineICAParams(resPath="mainz/", selCutoff=3, pvalCutoff=0.05)


###################################################
### code chunk number 8: libannot
###################################################
## load annotation package
library(hgu133a.db)


###################################################
### code chunk number 9: lspackage
###################################################
ls("package:hgu133a.db")


###################################################
### code chunk number 10: mart
###################################################
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")


###################################################
### code chunk number 11: martfilters
###################################################
listFilters(mart)[120:125,]


###################################################
### code chunk number 12: martattr
###################################################
listAttributes(mart)[grep(x=listAttributes(mart)[,1],pattern="affy")[1:5],]


###################################################
### code chunk number 13: buildIcaSet
###################################################
## Define typeID, Mainz data originate from affymetrix HG-U133a  microarray 
## and are indexed by probe sets.
## The probe sets are annotated into Gene Symbols
typeIDmainz <-  c(geneID_annotation="SYMBOL", geneID_biomart="hgnc_symbol", 
                    featureID_biomart="ensembl_gene_id")

## define the reference samples if any, here no normal sample is available
refSamplesMainz <- character(0)

resBuild <- buildIcaSet(params=params, A=data.frame(resJade$A), S=data.frame(resJade$S),
                        dat=mainz, pData=pData(mainz), refSamples=refSamplesMainz,
                        annotation="hgu133a.db", typeID= typeIDmainz,
                        chipManu = "affymetrix", mart=mart)

icaSetMainz <- resBuild$icaSet
params <- resBuild$params


###################################################
### code chunk number 14: showIcaSet
###################################################
icaSetMainz


###################################################
### code chunk number 15: defannot
###################################################
annot <- pData(icaSetMainz)


###################################################
### code chunk number 16: lookvar
###################################################
varLabels(icaSetMainz)[1:5]
icaSetMainz$grade[1:5]


###################################################
### code chunk number 17: lookfeatures
###################################################
featureNames(icaSetMainz)[1:5] # probe set ids
geneNames(icaSetMainz)[1:5] #gene symbols
sampleNames(icaSetMainz)[1:5] 


###################################################
### code chunk number 18: lookdat (eval = FALSE)
###################################################
## head(dat(icaSetMainz)) #probe set level
## head(datByGene(icaSetMainz)) #gene level


###################################################
### code chunk number 19: lookAS (eval = FALSE)
###################################################
## A(icaSetMainz)
## S(icaSetMainz)
## SByGene(icaSetMainz)


###################################################
### code chunk number 20: lookcomp (eval = FALSE)
###################################################
## nbComp(icaSetMainz)
## compNames(icaSetMainz)
## indComp(icaSetMainz)


###################################################
### code chunk number 21: witG
###################################################
witGenes(icaSetMainz)[1:5]
## We can for example modify the second contributing gene 
witGenes(icaSetMainz)[2] <- "KRT16"


###################################################
### code chunk number 22: MineICA.Rnw:362-363 (eval = FALSE)
###################################################
## compNames(icaSetMainz) <- paste("IC",1:nbComp(icaSetMainz),sep="")


###################################################
### code chunk number 23: MineICA.Rnw:371-379 (eval = FALSE)
###################################################
## ## select tumor samples of grade 3
## keepSamples <- sampleNames(icaSetMainz)[icaSetMainz$grade=="3"]
## ## Subset icaSetMainz to the grade-3 samples 
## icaSetMainz[,keepSamples]
## ## Subset icaSetMainz to the grade-3 samples and the first five components
## icaSetMainz[,keepSamples,1:5]
## ## Subset icaSetMainz to the first 10 features
## icaSetMainz[featureNames(icaSetMainz)[1:10],keepSamples]


###################################################
### code chunk number 24: histProj6
###################################################
hist(S(icaSetMainz)[,1], breaks=50, main="Distribution of feature projection on the first component", 
     xlab="projection values")
abline(v=c(3,-3), col="red", lty=2)


###################################################
### code chunk number 25: selectContribGenes
###################################################
## Extract the contributing genes
contrib <- selectContrib(icaSetMainz, cutoff=3, level="genes")
## Show the first contributing genes of the first and third components
sort(abs(contrib[[1]]),decreasing=TRUE)[1:10]
sort(abs(contrib[[3]]),decreasing=TRUE)[1:10]


###################################################
### code chunk number 26: selectContribGenes (eval = FALSE)
###################################################
## ## One can also want to apply different cutoffs depending on the components
## ## for example using the first 4 components:
## contrib <- selectContrib(icaSetMainz[,,1:4], cutoff=c(4,4,4,3), level="genes")


###################################################
### code chunk number 27: extractComp
###################################################
## extract sample contributions and gene projections of the second component
comp2 <- getComp(icaSetMainz, level="genes", ind=2)
## access the sample contributions 
comp2$contrib[1:5]
## access the gene projections
comp2$proj[1:5]


###################################################
### code chunk number 28: runAn
###################################################
## select the annotations of interest
varLabels(icaSetMainz)
# restrict the phenotype data to the variables of interest
keepVar <- c("age","er","grade")
# specify the variables that should be treated as character
icaSetMainz$er <- c("0"="ER-","1"="ER+")[as.character(icaSetMainz$er)]
icaSetMainz$grade <- as.character(icaSetMainz$grade)


###################################################
### code chunk number 29: runAn (eval = FALSE)
###################################################
## ## Run the analysis of the ICA decomposition
## # only enrichment in KEGG gene sets are tested
## runAn(params=params, icaSet=icaSetMainz, writeGenesByComp = TRUE, 
##       keepVar=keepVar, dbGOstats = "KEGG")


###################################################
### code chunk number 30: runAn
###################################################
resPath(params)


###################################################
### code chunk number 31: writeProj
###################################################
resW <- writeProjByComp(icaSet=icaSetMainz, params=params, mart=mart, 
                        level='genes', selCutoffWrite=2.5)
## the description of the contributing genes of each component is contained 
## in res$listAnnotComp which contains the gene id, its projection value, the number and 
## the indices of the components on which it exceeds the threshold, and its description.
head(resW$listAnnotComp[[1]])
## The number of components a gene contributes to is available 
## in res$nbOccInComp
head(resW$nbOccInComp)
## The output HTML files are located in the path:
genesPath(params)


###################################################
### code chunk number 32: plotHeatmaps (eval = FALSE)
###################################################
# ## selection of the variables we want to display on the heatmap
# keepVar <- c("er","grade")
# ## For the second component, select contributing genes using a threshold of 3 
# ## on the absolute projection values,
# ## heatmap with dendrogram
# resH <- plot_heatmapsOnSel(icaSet = icaSetMainz, selCutoff = 3, level = "genes", 
#                            keepVar = keepVar,
#                            doSamplesDendro = TRUE, doGenesDendro = TRUE, keepComp = 2,
#                            heatmapCol = maPalette(low = "blue", high = "red", mid = "yellow", k=44),
#                            file = "heatmapWithDendro", annot2col=annot2col(params))
# 
# 
# ## heatmap where genes and samples are ordered by contribution values
# resH <- plot_heatmapsOnSel(icaSet = icaSetMainz, selCutoff = 3, level = "genes", 
#                            keepVar = keepVar,
#                            doSamplesDendro = FALSE, doGenesDendro = FALSE, keepComp = 2,
#                            heatmapCol = maPalette(low = "blue", high = "red", mid = "yellow", k=44),
#                            file = "heatmapWithoutDendro", annot2col=annot2col(params))
# 


###################################################
### code chunk number 33: runEnrich (eval = FALSE)
###################################################
## ## run enrichment analysis on the first three components of icaSetMainz, 
## ## using gene sets from the ontology 'Biological Process' (BP) of Gene Ontology (GO)
## resEnrich <- runEnrich(params=params,icaSet=icaSetMainz[,,1:3],
##                        dbs=c("GO"), ontos="BP")


###################################################
### code chunk number 34: runEnrich3 (eval = FALSE)
###################################################
## ## Access results obtained for GO/BP for the first three components 
## # First component, when gene selection was based on the negative projection values
## head(resEnrich$GO$BP[[1]]$left)


###################################################
### code chunk number 35: runEnrich3bis
###################################################
structure(list(GOBPID = c("GO:0006955", "GO:0002694", "GO:0050867", 
"GO:0002429", "GO:0050863", "GO:0051251"), Pvalue = c(4.13398630676443e-16, 
3.53565704068228e-14, 1.0598364083037e-11, 2.10474926262594e-11, 
2.96408182272612e-11, 3.41481775787203e-11), OddsRatio = c(16.1139281129653, 
10.2399932157395, 10.34765625, 13.1975308641975, 14.5422535211268, 
14.3646536754775), ExpCount = c(2.18081918081918, 3.24691805656273, 
2.3821609862219, 1.56635242929659, 1.33711359063472, 1.35043827611395
), Count = c(21L, 23L, 18L, 15L, 14L, 14L), Size = c(185L, 199L, 
146L, 96L, 85L, 85L), Term = c("immune response", "regulation of leukocyte activation", 
"positive regulation of cell activation", "immune response-activating cell surface receptor signaling pathway", 
"regulation of T cell activation", "positive regulation of lymphocyte activation"
), In_geneSymbols = c("CD7,MS4A1,CD27,CTSW,GZMA,HLA-DOB,IGHD,IGHM,IGJ,IL2RG,CXCL10,LTB,LY9,CXCL9,CCL18,CXCL11,CST7,FCGR2C,IL32,ADAMDEC1,ICOS", 
"AIF1,CD2,CD3D,CD3G,CD247,CD27,CD37,CD38,HLA-DQB1,LCK,PTPRC,CCL5,CCL19,XCL1,EBI3,LILRB1,PTPN22,SIT1,TRBC1,TRAC,ICOS,MZB1,SLAMF7", 
"AIF1,CD2,CD3D,CD3G,CD247,CD27,CD38,HLA-DQB1,LCK,PTPRC,CCL5,CCL19,XCL1,EBI3,LILRB1,TRBC1,TRAC,ICOS", 
"CD3D,CD3G,CD247,CD38,HLA-DQB1,IGHG1,IGKC,IGLC1,LCK,PRKCB,PTPRC,PTPN22,TRBC1,TRAC,TRAT1", 
"AIF1,CD3D,CD3G,CD247,HLA-DQB1,PTPRC,CCL5,XCL1,EBI3,PTPN22,SIT1,TRBC1,TRAC,ICOS", 
"AIF1,CD3D,CD3G,CD247,CD38,HLA-DQB1,PTPRC,CCL5,XCL1,EBI3,LILRB1,TRBC1,TRAC,ICOS"
)), .Names = c("GOBPID", "Pvalue", "OddsRatio", "ExpCount", "Count", 
"Size", "Term", "In_geneSymbols"), row.names = c(NA, 6L), class = "data.frame")


###################################################
### code chunk number 36: runEnrich4 (eval = FALSE)
###################################################
## # Second component
## head(resEnrich$GO$BP[[2]]$both, n=5)


###################################################
### code chunk number 37: runEnrich4bis
###################################################
structure(list(GOBPID = c("GO:0045104", "GO:0031581", "GO:0030318", 
"GO:0070488", "GO:0072602", "GO:0034329"), Pvalue = c(2.16044773513962e-05, 
6.04616151867683e-05, 0.000394387232705895, 0.000461592979000511, 
0.000461592979000511, 0.00107959102524193), OddsRatio = c(19.6820175438597, 
26.7826086956522, 14.4053511705686, Inf, Inf, 4.29841077032088
), ExpCount = c(0.366751269035533, 0.237309644670051, 0.366751269035533, 
0.0431472081218274, 0.0431472081218274, 2.09263959390863), Count = c(5L, 
4L, 4L, 2L, 2L, 8L), Size = c(17L, 11L, 17L, 2L, 2L, 97L), Term = c("intermediate filament cytoskeleton organization", 
"hemidesmosome assembly", "melanocyte differentiation", "neutrophil aggregation", 
"interleukin-4 secretion", "cell junction assembly"), In_geneSymbols = c("DST,KRT14,KRT16,PKP1,SYNM", 
"DST,COL17A1,KRT5,KRT14", "EDN3,KIT,SOX10,MLPH", "S100A8,S100A9", 
"GATA3,VTCN1", "DST,CDH3,COL17A1,GPM6B,KRT5,KRT14,SFRP1,UGT8"
)), .Names = c("GOBPID", "Pvalue", "OddsRatio", "ExpCount", "Count", 
"Size", "Term", "In_geneSymbols"), row.names = c(NA, 6L), class = "data.frame")


###################################################
### code chunk number 38: runEnrich2 (eval = FALSE)
###################################################
## # Third component, when gene selection was based  on the absolute projection values
## head(resEnrich$GO$BP[[3]]$both)


###################################################
### code chunk number 39: runEnrich2bis
###################################################
structure(list(GOBPID = c("GO:0048285", "GO:0051301", "GO:0007067", 
"GO:0007076", "GO:0000086", "GO:0006271"), Pvalue = c(2.1806594780167e-24, 
7.64765765268202e-16, 4.85530963990747e-12, 9.30475956146815e-06, 
1.95701394548966e-05, 5.65649365383205e-05), OddsRatio = c(16.8486540378863, 
14.5290697674419, 17.0874279123414, Inf, 8.18195718654434, 27.2667509481669
), ExpCount = c(3.1816533720087, 2.14748665070889, 1.18268700606506, 
0.063633067440174, 1.18781725888325, 0.233321247280638), Count = c(32L, 
21L, 14L, 3L, 8L, 4L), Size = c(150L, 107L, 65L, 3L, 56L, 11L
), Term = c("organelle fission", "cell division", "mitosis", 
"mitotic chromosome condensation", "G2/M transition of mitotic cell cycle", 
"DNA strand elongation involved in DNA replication"), In_geneSymbols = c("BIRC5,BUB1,CCNA2,CDK1,CDC20,CDC25A,CENPE,IGF1,KIFC1,MYBL2,NEK2,AURKA,CCNB2,KIF23,DLGAP5,NDC80,UBE2C,TPX2,NCAPH,UBE2S,NUSAP1,ERCC6L,CDCA8,CEP55,CENPN,PBK,NCAPG,DSCC1,CDCA3,KIF18B,SKA1,ASPM", 
"BUB1,CCNA2,CDK1,CDC20,CDC25A,CENPE,KIFC1,NEK2,AURKA,TOP2A,CCNB2,NDC80,UBE2C,TPX2,NCAPH,UBE2S,ERCC6L,CDCA8,NCAPG,CDCA3,SKA1", 
"CCNA2,CDK1,CDC25A,CCNB2,TPX2,ERCC6L,CDCA8,CEP55,CENPN,PBK,CDCA3,KIF18B,SKA1,ASPM", 
"NCAPH,NUSAP1,NCAPG", "BIRC5,CCNA2,CDK1,CDC25A,FOXM1,NEK2,CCNB2,MELK", 
"MCM5,CDC45,GINS1,GINS2")), .Names = c("GOBPID", "Pvalue", "OddsRatio", 
"ExpCount", "Count", "Size", "Term", "In_geneSymbols"), row.names = c(NA, 
6L), class = "data.frame")


###################################################
### code chunk number 40: stagebis
###################################################
icaSetMainz$grade <- factor(as.character(icaSetMainz$grade), levels=c("1","2","3"))


###################################################
### code chunk number 41: varAn
###################################################
### Qualitative variables
## Compute Wilcoxon and Kruskall-Wallis tests to compare the distribution 
## of the samples according to their grade and ER status on all components.
resQual <- qualVarAnalysis(params=params, icaSet=icaSetMainz, 
                           keepVar=c("er","grade"),
                           adjustBy="none", typePlot="boxplot", 
                           path="qualVarAnalysis/", filename="qualVar")


###################################################
### code chunk number 42: quantVarAn
###################################################
### Quantitative variables
## Compute pearson correlations between variable 'age' and the sample contributions
## on all components.
## We are interested in correlations exceeding 0.3 in absolute value, and plots will only be drawn
## for correlations exceeding this threshold.
resQuant <- quantVarAnalysis(params=params, icaSet=icaSetMainz, keepVar="age", 
                             typeCor="pearson", cutoffOn="cor",
                             cutoff=0.3, adjustBy="none",  
                             path="quantVarAnalysis/", filename="quantVar")


###################################################
### code chunk number 43: quantVarAncor
###################################################
resQuant$cor[2]


###################################################
### code chunk number 44: plotMix
###################################################
resmix <- plotAllMix(A=A(icaSetMainz), nbMix=2, nbBreaks=50)


###################################################
### code chunk number 45: stageOnHist (eval = FALSE)
###################################################
## ## plot the positions of the samples on the second component according to their ER status
## ## (in a file "er.pdf") 
## plotPosAnnotInComp(icaSet=icaSetMainz, params=params, keepVar=c("er"), keepComp=2,  
##                    funClus="Mclust")


###################################################
### code chunk number 46: cluster1
###################################################
## clustering of the samples in 1-dim using the vector 
## of sample contributions of the two first components 
## and Gaussian mixture modeling (Mclust)
clus1 <- clusterSamplesByComp(params=params, icaSet=icaSetMainz[,,,1:2],
                              funClus="Mclust", clusterOn="A", nbClus=2, filename="comp1Mclust")

## The obtained clusters are written in the file "comp1Mclus.txt" of the result path.
clus1$clus[[2]][1:5]


###################################################
### code chunk number 47: cluster2
###################################################
clus2 <- clusterSamplesByComp_multiple(params=params, icaSet=icaSetMainz[,,1:2], 
                                     funClus="kmeans", clusterOn=c("A","S"), level="features", 
                                     nbClus=2, filename="comparKmeans")

## The obtained clusters and their comparison with adjusted Rand indices are written 
## in file "comparKmeans.txt" of the result path.

## Both clustering results are stored in a common data.frame
head(clus2$clus)
## Access Rand index
clus2$comparClus


###################################################
### code chunk number 48: clus2var (eval = FALSE)
###################################################
## ## Test the association between the clustering obtained by Mclust for the first
## ## component and the variables:  
clus2var <- clusVarAnalysis(icaSet=icaSetMainz[,,1:2], params=params, 
                            keepVar=c("er","grade"), 
                            resClus=clus1$clus, funClus="Mclust", adjustBy="none", 
                            doPlot=TRUE, path="clus2var/", filename="resChitests-Mcluscomp1")
## Look at the filename which contains p-values and links to the barplots
## p-values are also contained in the ouput of the function:
clus2var
## 
## 


###################################################
### code chunk number 49: clus2var (eval = FALSE)
###################################################
## structure(list(`1` = c(1.657e-06, 2.315e-08), `2` = c(6.775e-07, 
## 0.0001899)), .Names = c("1", "2"), row.names = c("er", "grade"
## ), class = "data.frame")


###################################################
### code chunk number 50: tuut (eval = FALSE)
###################################################
## load three other breast cancer datasets also based on  Affymetrix HG-U133a microarray
library(breastCancerUPP)
library(breastCancerTRANSBIG)
library(breastCancerVDX)
data(upp)
data(transbig)
data(vdx)
## function to build IcaSet instances from these three datasets
treat <- function(es, annot="hgu133a.db") {
    es <- selectFeatures_IQR(es,10000)
    exprs(es) <- t(apply(exprs(es),1,scale,scale=FALSE))
    colnames(exprs(es)) <- sampleNames(es)
    resJade <- runICA(X=exprs(es), nbComp=5, method = "JADE", maxit=10000)
    resBuild <- buildIcaSet(params=buildMineICAParams(), A=data.frame(resJade$A), S=data.frame(resJade$S),
                            dat=exprs(es), pData=pData(es), refSamples=character(0),
                            annotation=annot, typeID= typeIDmainz,
                            chipManu = "affymetrix", mart=mart)
    icaSet <- resBuild$icaSet
}
icaSetUpp <- treat(upp, annot="hgu133plus2.db")
icaSetVdx <- treat(vdx)
icaSetTransbig <- treat(transbig)


###################################################
### code chunk number 51: compareAn (eval = FALSE)
###################################################
## resGraph <- runCompareIcaSets(icaSets=list(icaSetMainz, icaSetUpp, 
##                               icaSetTransbig, icaSetVdx), 
##                               labAn=c("Mainz", "Upp","Transbig","Vdx"), 
##                               type.corr="pearson", level="genes", 
##                               cutoff_zval=0, fileNodeDescr="nodeDescr.txt", 
##                               fileDataGraph="dataGraph.txt", tkplot=TRUE)


###################################################
### code chunk number 52: compareAnColors (eval = FALSE)
###################################################
## barplot(names.arg=unique(resGraph[[2]]$labAn),height=rep(1,4), 
##         col=unique(resGraph[[2]]$col))


###################################################
### code chunk number 53: inter1 (eval = FALSE)
###################################################
## ## comparison of four components included in the clique of the correlation-based graph
## # that includes the second component of Mainz.
## inter <- compareGenes(keepCompByIcaSet =   c(2,2,2,2),  
##                       icaSets = list(icaSetMainz, icaSetTransbig, icaSetUpp, icaSetVdx),
##                       lab=c("Mainz", "Transbig", "Upp", "Vdx"), cutoff=3, 
##                       type="intersection",  annotate=F)
## head(inter)


###################################################
### code chunk number 54: inter1bis
###################################################
structure(list(min_rank = structure(c(1L, 12L, 22L, 12L, 31L, 
31L), .Label = c("1", "10", "101", "11", "12", "13", "14", "16", 
"17", "18", "19", "2", "20", "21", "22", "23", "24", "25", "26", 
"27", "29", "3", "30", "31", "33", "34", "35", "36", "37", "38", 
"4", "40", "41", "43", "46", "47", "5", "51", "52", "53", "54", 
"58", "6", "63", "69", "7", "70", "74", "8", "9"), class = "factor"), 
    median_rank = c(1, 2, 3, 5, 5.5, 6), ranks = c("1,1,1,1", 
    "2,2,3,2", "3,3,4,3", "5,8,2,5", "4,7,9,4", "6,4,6,6"), scaled_proj = c("-8.9,-9.2,9.2,-9.3", 
    "-8.4,-8.7,7.9,-8.9", "-7,-8.3,7.7,-8.8", "-6.4,-6.1,8.4,-7.6", 
    "-6.6,-6.5,6.7,-7.7", "-6.3,-7.8,7.1,-7.6")), .Names = c("min_rank", 
"median_rank", "ranks", "scaled_proj"), row.names = c("IGL@", 
"IGKV4-1", "IGLV2-23", "NKG7", "IGHM", "TNFRSF17"), class = "data.frame")


###################################################
### code chunk number 55: inter2 (eval = FALSE)
###################################################
## ## comparison of four components included in the clique of the correlation-based graph
## # that includes the third component of Mainz.
## inter <- compareGenes(keepCompByIcaSet = c(3,3,3,1), 
##                       icaSets = list(icaSetMainz, icaSetTransbig, icaSetUpp, icaSetVdx),
##                       lab=c("Mainz", "Transbig", "Upp", "Vdx"), cutoff=3, 
##                       type="intersection",  annotate=F)
## head(inter)


###################################################
### code chunk number 56: inter2bis
###################################################
structure(list(min_rank = structure(c(1L, 7L, 21L, 7L, 7L, 16L
), .Label = c("1", "10", "11", "16", "17", "18", "2", "23", "24", 
"25", "26", "3", "32", "38", "39", "4", "40", "42", "43", "49", 
"5", "51", "58", "6", "7", "8", "9"), class = "factor"), median_rank = c(1, 
4, 6.5, 8.5, 9.5, 9.5), ranks = c("1,1,1,3", "5,5,3,2", "7,7,6,5", 
"13,2,4,21", "2,15,15,4", "4,31,12,7"), scaled_proj = c("7.7,8.7,8.1,7.2", 
"6.3,7.2,7.3,7.3", "6.2,7.1,5.9,6.5", "5.7,7.3,7,5.1", "7.3,5.8,5.4,7", 
"6.3,4.6,5.6,6.2")), .Names = c("min_rank", "median_rank", "ranks", 
"scaled_proj"), row.names = c("GABRP", "MIA", "SERPINB5", "KRT14", 
"KRT16", "KRT81"), class = "data.frame")




