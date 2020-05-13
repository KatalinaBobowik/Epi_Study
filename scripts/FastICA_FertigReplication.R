# created by KSB, 09.03.20
# WGCNA analysis for variables associated with expression profiles
# using the tutorial as per: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html

### R code from vignette source 'MineICA.Rnw'
dana <- batch.corrected.lcpm

##############
eigs <- princomp(cor(dana))$sdev^2
eigs <- eigs/sum(eigs)

barplot(eigs[1:10], ylab='Proportion of variance')

sum(eigs[1:2])
## [1]  0.7103975
par(mfrow=c(1,2))
pca2d(princomp(cor(dana)),col=as.numeric(as.factor(y$samples$Island)))
title('PCA colored by Island')

# pca2d(princomp(cor(dana)),col=as.numeric(as.factor(pd.GTExGeneRPKM.SUB10Brain$SubjectID)))
# title('PCA colored by GTeX subject')

###########

dana <- dana[apply(dana,1,sd)>0,]

ica <- list()
for (p in 1:10) {
  ica[[paste0('p',p)]] <- fastICA(dana,n.comp = p)
}

ica <- fastICA(dana,n.comp = 10)

par(mfrow=c(2,5))
for (i in 1:10) {
  plot(as.numeric(as.factor(y$samples$Eukaryota_Plasmodiidae)),ica$A[i,], 
       col=as.numeric(as.factor(y$samples$Island)),axes=F,pch=19)
  axis(2)
  axis(1,at=1:length(levels(as.factor(y$samples$Eukaryota_Plasmodiidae))),
       labels=levels(as.factor(y$samples$Eukaryota_Plasmodiidae)),las=2)
}

# get variable names to make associations
for (name in colnames(OTUs)){
 assign(name, y$samples[[paste0(name)]])
}

plot.ica <- function(dataToICA, speciesCol, namesPch, sampleNames){
	ics <- fastICA(dataToICA,n.comp=10)
	for (i in 1:9){
        ica_axis1=i
        ica_axis2=i+1
        plot(ics$A[,ica_axis1], ics$A[,ica_axis2], col=speciesCol, pch=namesPch, cex=2, xlab=paste0("IC", ica_axis1), ylab=paste0("IC",ica_axis1), main=name)
     }

    return(ica)
}

indoSampleSet <- t(apply(batch.corrected.lcpm,1,scale,scale=F))
for (name in colnames(OTUs)){
	initial = .bincode(get(name), breaks=seq(min(get(name), na.rm=T), max(get(name), na.rm=T), len = 80),include.lowest = TRUE)
    bloodCol <- colorRampPalette(c("blue", "red"))(79)[initial]
    #pcaresults <- plot.pca(dataToPca=batch.corrected.lcpm, speciesCol=bloodCol,namesPch=as.numeric(y$samples$Island)+14,sampleNames=get(name))
	#plot(ics$A[1,],ics$A[2,],col=bloodCol,xlab="IC1",ylab="IC2")
	pdf(paste0(outputdir,"icaresults_",name,".pdf"))
	plot.ica(dataToICA=indoSampleSet, speciesCol=bloodCol, namesPch=19, sampleNames=name)
	dev.off()
}

map=c(SMBvsMPI,MTWvsMPI)
#map=y$genes[unique(c(SMBvsMPI,MTWvsMPI)),]$SYMBOL)
length(contrib[[4]][which(names(contrib[[4]]) %in% map)])/length(map)

#ment=make.unique(y$genes[unique(c(SMBvsMTW,MTWvsMPI)),]$SYMBOL)
ment=c(SMBvsMTW,MTWvsMPI)
length(contrib[[4]][which(names(contrib[[4]]) %in% ment)])/length(ment)

#sumb=make.unique(y$genes[unique(c(SMBvsMPI,SMBvsMTW)),]$SYMBOL)
sumb=c(SMBvsMPI,SMBvsMTW)
length(contrib[[4]][which(names(contrib[[4]]) %in% sumb)])/length(sumb)

####
contributingGenes=matrix(nrow=3,ncol=5)
rownames(contributingGenes)=c("SMBvsMPI","MTWvsMPI","SMBvsMTW")
colnames(contributingGenes)=c("PC1","PC2","PC3","PC4","PC5")
for (i in 1:5){
	contributingGenes[1,i]=length(contrib[[i]][which(names(contrib[[i]]) %in% SMBvsMPI)])
	#contributingGenes[1,2]=length(SMBvsMPI)
	contributingGenes[2,i]=length(contrib[[i]][which(names(contrib[[i]]) %in% MTWvsMPI)])
	#contributingGenes[2,2]=length(MTWvsMPI)
	contributingGenes[3,i]=length(contrib[[i]][which(names(contrib[[i]]) %in% SMBvsMTW)])
}
contributingGenes=as.data.frame(contributingGenes)
contributingGenes$type="contribGenes"
contributingGenes$island=rownames(contributingGenes)

c=do.call("cbind", replicate(5, c(length(SMBvsMPI),length(MTWvsMPI),length(SMBvsMTW)), simplify = FALSE))
c=as.data.frame(c)
rownames(c)=c("SMBvsMPI","MTWvsMPI","SMBvsMTW")
colnames(c)=c("PC1","PC2","PC3","PC4","PC5")
c$type="deGenes"
c$island=rownames(c)
allGenes=rbind(contributingGenes,c)
meltedGenes=melt(allGenes)

ggplot(meltedGenes, aes(x = island, y = value, fill = type)) + 
  geom_bar(stat = 'identity', position = 'stack') + facet_grid(~ variable)

# g <- ggplot(x, aes(islandComp)) + geom_bar(aes(fill = value))
# ggplot(y, aes(x = islandComp, y = value, fill = variable)) + geom_bar(stat = "identity")
