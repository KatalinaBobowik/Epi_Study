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
library("factoextra")

# Set paths:
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/dataPreprocessing/"
refDir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/"
outputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/"

# 1. load in expression data, log transform, then mean center data ----------------------------------------------------------------

# Fucking rad!!!
# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/119-pca-in-r-using-ade4-quick-scripts/

# load expression data and OTU file
load(paste0(inputdir, "indoRNA.read_counts.TMM.filtered.Rda"))
OTUs=read.table(paste0(refDir,"scaledOTUs.txt"))
#OTUs=read.table(paste0(refDir,"OTUtable.txt"))

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

# Co-inertia analysis between Mappi samples gene expression levels and pathogen load -------------------

# these are all DE genes at no log fold change threshold and an adjusted p-value of 0.05
SMBvsMTW=read.table("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/DE_Island/LM_allCovarPlusBlood/topTable_SMBvsMTW.txt")
SMBvsMTW=rownames(SMBvsMTW)[which(abs(SMBvsMTW$logFC) >= 1)]
SMBvsMPI=read.table("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/DE_Island/LM_allCovarPlusBlood/topTable_SMBvsMPI.txt")
SMBvsMPI=rownames(SMBvsMPI)[which(abs(SMBvsMPI$logFC) >= 1)]
MTWvsMPI=read.table("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/DE_Island/LM_allCovarPlusBlood/topTable_MTWvsMPI.txt")
MTWvsMPI=rownames(MTWvsMPI)[which(abs(MTWvsMPI$logFC) >= 1)]

# center and reduce to one unit of variance prior to performing the PCA analysis, as suggested by Fave et al.
indoSampleSet <- t(apply(batch.corrected.lcpm,1,scale,scale=F))

# First assign environmental variables and genes
env=y$samples[,colnames(OTUs)]
#env=env[grep("MPI", rownames(env)),]
genes=indoSampleSet[unique(c(SMBvsMPI,MTWvsMPI)),]
#genes=indoSampleSet
genes=t(genes)
colnames(genes)=y[colnames(genes),]$genes$SYMBOL
colnames(genes)=make.unique(colnames(genes))
rownames(genes) = colnames(lcpm)
genes=genes[,grep("RNF182|MDGA1|TMTC1|TUBB2A|SIGLEC14|MARCO|MYOM2|KCNMA1|GRB14|RAP1GAP|UTS2|SOX5|OCLN",colnames(genes))]
#genes=genes[grep("MPI", rownames(genes)),]
# Set up Coinertia analysis -
# dudi1 and 2 do not need to be scaled since they are all in the same units 
# i.e., the environmental variables are all in the same units and the gene expression
# data is all in the same unit
dudi1 <- dudi.pca(env, scale = F, scan = F, nf = 2)
dudi2 <- dudi.pca(genes, scale = F, scan = F, nf = 3)
pdf(paste0(outputdir,"DudiPCA_ENV.pdf"), width=15)
Heatmap(dudi1$tab)
dev.off()
coin1 <- coinertia(dudi1,dudi2, scan = F, nf = 2)
# plot coinartia analysis
pdf(paste0(outputdir,"coinertiaANalysis.pdf"), width=15)
plot(coin1)
dev.off()

rv1 <- RV.rtest(dudi1$tab, dudi2$tab, 99)
plot(rv1)
pdf(paste0(outputdir,"correlationHeatmap.pdf"), width=15)
Heatmap(t(coin1$tab), col=colorRamp2(c(-0.2, 0, 0.2), c("blue", "white", "red")))
dev.off()
names=names(which(abs(colSums(coin1$tab))>0))
Heatmap(t(coin1$tab[,which(abs(colSums(coin1$tab))>0)]),col=colorRamp2(c(-1, 0, 1), c("blue", "white", "red")))

# factoextra pipeine

pdf(paste0(outputdir,"indDudiPCA.pdf"))
fviz_pca_ind(dudi1,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = F     # Avoid text overlapping
             )
dev.off()

pdf(paste0(outputdir,"envDudiPCA.pdf"))
fviz_pca_var(dudi1,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = F     # Avoid text overlapping
             )
dev.off()

pdf(paste0(outputdir,"envAndindDudiPCA.pdf"))
fviz_pca_biplot(dudi1, repel = F,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
                )

dev.off()


pdf(paste0(outputdir,"Dudi2_indDudiPCA.pdf"))
fviz_pca_ind(dudi2,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = F     # Avoid text overlapping
             )
dev.off()

pdf(paste0(outputdir,"Dudi2_envDudiPCA.pdf"))
fviz_pca_var(dudi2,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = F     # Avoid text overlapping
             )
dev.off()

pdf(paste0(outputdir,"Dudi2_envAndindDudiPCA.pdf"))
fviz_pca_biplot(dudi2, repel = F,
                col.var = "contrib", # Variables color
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                col.ind = "#696969"  # Individuals color
                )
dev.off()
# PCA analysis ----------------------------------------

# Is there structure/variance between MPI and other islands?

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

# get variable names to make associations
for (name in colnames(OTUs)){
 assign(name, y$samples[[paste0(name)]])
}

# Prepare covariate matrix
all.covars.df <- y$samples[,c("Island","Mappi","Age","RIN",colnames(OTUs))]

for (name in colnames(OTUs)){
    initial = .bincode(get(name), breaks=seq(min(get(name), na.rm=T), max(get(name), na.rm=T), len = 80),include.lowest = TRUE)
    bloodCol <- colorRampPalette(c("blue", "red"))(79)[initial]
    pdf(paste0(outputdir,"pcaresults_",name,".pdf"))
    pcaresults <- plot.pca(dataToPca=batch.corrected.lcpm, speciesCol=bloodCol,namesPch=as.numeric(y$samples$Island)+14,sampleNames=get(name))
    #legend(legend=c("High","Low"), pch=16, x="bottomright", col=c(bloodCol[which.max(get(name))], bloodCol[which.min(get(name))]), cex=0.6, title=name, border=F, bty="n")
    #legend(legend=unique(as.numeric(y$samples$batch)), "topright", pch=unique(as.numeric(y$samples$batch)) + 14, title="Batch", cex=0.6, border=F, bty="n")
    dev.off()
}

# Get PCA associations
all.pcs <- pc.assoc(pcaresults)
all.pcs$Variance <- pcaresults$sdev^2/sum(pcaresults$sdev^2)
write.table(all.pcs, file=paste0(outputdir,"pca_covariates_blood_RNASeqDeconCell.txt"), col.names=T, row.names=F, quote=F, sep="\t")

