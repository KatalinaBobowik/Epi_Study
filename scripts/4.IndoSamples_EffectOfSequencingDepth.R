# created by KSB, 20.01.19
# Script created to inspect where reads are mapping to. 

# load packages
require(ggplot2)
require(RColorBrewer)
library(dplyr)
library(plyr)
library(phyloseq)
library(reshape2)
library(ggpubr)
library(vegan)
library(metacoder)

# set ggplot colour theme to white
theme_set(theme_bw())

# set up directories
AllReadsdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/All_Reads/"
Indo250Kdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/Indo_250K/"
outputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/Rarefaction_Singletons/"

#################
### All Reads ###
#################

# load in Indonesian phyloseq object
load(paste0(AllReadsdir,"allReads_Counts_physeq.Rda"))

# First, decide on how many OTUs should be filtered out
allReads <- otu_table(allReads_Counts_physeq)

# Try out different filtering methods:
for (i in 1:5){
	# remove OTUs at or under thresholding of 1:5
	allReads[allReads<=i]<-0
	# rarefaction curve to check the effect of sampling depth on species found
	# here, we'll do it with labels. Note, you need to make the DF into an otu table, since it's
	# a phyloseq object
	pdf(paste0(outputdir,"IndoRarefactionCurve_AllReads_NoSpeciesFiltering_RemovingCounts_",i,"_andBelow_Counts.pdf"), width=10)
	rarecurve(t(otu_table(allReads, taxa_are_rows = TRUE)), step=20, label=T, xlab="Counts",ylab="Number Species", main=paste("Removing Reads",i,"And Below",sep=" "))
	dev.off()

	# rarefaction curve to check the effect of sampling depth on species found
	# without labels
	pdf(paste0(outputdir,"IndoRarefactionCurve_AllReads_NoSpeciesFiltering_RemovingCounts_",i,"_andBelow_Counts_noLabels.pdf"), width=10)
	rarecurve(t(otu_table(allReads, taxa_are_rows = TRUE)), step=20, label=F, xlab="Counts",ylab="Number Species", main=paste("Removing Reads",i,"And Below",sep=" "))
	dev.off()
}

#######################
### 250K SubSampled ###
#######################

# load in 250K Indonesian phyloseq object
load(paste0(Indo250Kdir,"subSampled_250K_Indo_Counts_physeq.Rda"))

# Separate species' abundances and taxonomy columns
subSamp_250K <- otu_table(subSampled_250K_Indo_Counts_physeq)

# Try out different filtering methods:
for (i in 1:5){
	# remove OTUs at or under thresholding of 1:5
	# get rid of singletons
	subSamp_250K[subSamp_250K<=i]<-0

	# rarefaction curve to check the effect of sampling depth on species found
	# here, we'll do it with labels. Note, you need to make the DF into an otu table, since it's
	# a phyloseq object
	pdf(paste0(outputdir,"IndoRarefactionCurve_250K_NoSpeciesFiltering_RemovingCounts_",i,"_andBelow_Counts.pdf"), width=10)
	rarecurve(t(otu_table(subSamp_250K, taxa_are_rows = TRUE)), step=20, label=T, xlab="Counts",ylab="Number Species",main=paste("Removing Reads",i,"And Below",sep=" "))
	dev.off()

	# rarefaction curve to check the effect of sampling depth on species found
	# without labels
	pdf(paste0(outputdir,"IndoRarefactionCurve_250K_NoSpeciesFiltering_RemovingCounts_",i,"_andBelow_Counts_noLabels.pdf"), width=10)
	rarecurve(t(otu_table(subSamp_250K, taxa_are_rows = TRUE)), step=20, label=F, xlab="Counts",ylab="Number Species",main=paste("Removing Reads",i,"And Below",sep=" "))
	dev.off()
}

#########################
###### Mantel Test ######
#########################

# For 250K Reads -----------------------
subSamp_250K <- otu_table(subSampled_250K_Indo_Counts_physeq)
# copy the all reads dataframe, as you will be altering it and need to compare to the original 
subSamp_250K_noFilter=subSamp_250K

# make an empty matrix to save data to
mantel_result_matrix=matrix(nrow=5,ncol=2)
colnames(mantel_result_matrix)=c("Spearman_Rho","Significance")
rownames(mantel_result_matrix)=c(1:5)
for (i in 1:5){
	# get rid of anything with 5 counts or less
	subSamp_250K[subSamp_250K<=i]<-0
	dist.250K=vegdist(subSamp_250K, method = "bray")
	dist.250K.Nofilter=vegdist(subSamp_250K_noFilter, method = "bray")
	mantel_result=mantel(dist.250K, dist.250K.Nofilter, method = "spearman", permutations = 999, na.rm = TRUE)
	mantel_result_matrix[i,]=c(mantel_result$statistic, mantel_result$signif)
}
mantel_result_matrix <- cbind(mantel_result_matrix, "nFiltered"=1:nrow(mantel_result_matrix)) 
write.table(mantel_result_matrix,file=paste0(outputdir,"250K_Mantel_result_matrix_BraysDissimilarity_SpearmanRho.txt"),sep="\t")

# For All Reads -----------------------------
allReads <- otu_table(allReads_Counts_physeq)
# copy the all reads dataframe, as you will be altering it and need to compare to the original 
allReads_noFilter=allReads

# make an empty matrix to save data to
mantel_result_matrix=matrix(nrow=5,ncol=2)
colnames(mantel_result_matrix)=c("Spearman_Rho","Significance")
rownames(mantel_result_matrix)=c(1:5)
for (i in 1:5){
	# get rid of anything with 5 counts or less
	allReads[allReads<=i]<-0
	dist.allReads=vegdist(allReads, method = "bray")
	dist.allReads.Nofilter=vegdist(allReads_noFilter, method = "bray")
	mantel_result=mantel(dist.allReads, dist.allReads.Nofilter, method = "spearman", permutations = 999, na.rm = TRUE)
	mantel_result_matrix[i,]=c(mantel_result$statistic, mantel_result$signif)
}
mantel_result_matrix <- cbind(mantel_result_matrix, "nFiltered"=1:nrow(mantel_result_matrix)) 
write.table(mantel_result_matrix,file=paste0(outputdir,"Mantel_result_matrix_BraysDissimilarity_SpearmanRho.txt"),sep="\t")

#########################
### Exploratory Plots ###
#########################

# Filtered data ------------------

# Make melted dataframe of number of reads and species for all reads

# remove empty rows
allReads=allReads[rowSums(allReads)>0, ]

depth=colSums(allReads)
nSpecies=colSums(allReads!=0)
speciesSummary=data.frame(depth,nSpecies)
speciesSummary$librarySize="AllReads"
speciesSummary$samples=rownames(speciesSummary)

# Make melted dataframe of number of reads and species for subsampled data
# remove empty rows
subSamp_250K=subSamp_250K[rowSums(subSamp_250K)>0, ]

depth_subSampled=colSums(subSamp_250K)
nSpecies_subSampled=colSums(subSamp_250K!=0)
speciesSummary_subSampled=data.frame(depth_subSampled,nSpecies_subSampled)
speciesSummary_subSampled$librarySize="SubSampled"
colnames(speciesSummary_subSampled)=c("depth","nSpecies","librarySize")
speciesSummary_subSampled$samples=rownames(speciesSummary_subSampled)

# bind dataframes together
summaryTable=rbind(speciesSummary,speciesSummary_subSampled)

# correlation of sequencing depth on number of species found
pdf(paste0(outputdir,"DepthAndnSpecies_CorrelationPlots.pdf"), width=10)
ggplot(data=summaryTable, aes(x=depth, y=nSpecies,group=librarySize)) +
  geom_point(aes(color=librarySize)) + geom_smooth(method="lm") + stat_cor()
dev.off()

# histogram of the number of species by sequencing depth
pdf(paste0(outputdir,"HistogramLibrarySize_DepthComparison.pdf"), width=10)
ggplot(data=summaryTable, aes(x=depth)) +
	geom_histogram(aes(color = librarySize, fill= librarySize), bins = 80, alpha=0.3, position="identity")
dev.off()

pdf(paste0(outputdir,"HistogramLibrarySize_nSpeciesComparison.pdf"), width=10)
ggplot(data=summaryTable, aes(x=nSpecies)) +
	geom_histogram(aes(color = librarySize, fill= librarySize), bins = 10, alpha=0.3, position="identity")
dev.off()


