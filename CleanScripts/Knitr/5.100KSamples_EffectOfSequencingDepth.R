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
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/100KSamples/"
outputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/Rarefaction_Singletons_100K/"

#################
### Indo Data ###
#################

# load in 100K Indonesian phyloseq object
load(paste0(inputdir,"subSampled_100K_Indo_Counts_physeq.Rda"))

# First, decide on how many OTUs should be filtered out
Indonesian <- otu_table(subSampled_100K_Indo_Counts_physeq)

# Try out different filtering methods:
for (i in 1:5){
	# remove OTUs at or under thresholding of 1:5
	Indonesian[Indonesian<=i]<-0
	# rarefaction curve to check the effect of sampling depth on species found
	# here, we'll do it with labels. Note, you need to make the DF into an otu table, since it's
	# a phyloseq object
	pdf(paste0(outputdir,"IndoRarefactionCurve_100K_NoSpeciesFiltering_RemovingCounts_",i,"_andBelow_Counts.pdf"), width=10)
	rarecurve(t(otu_table(Indonesian, taxa_are_rows = TRUE)), step=20, label=T, xlab="Counts",ylab="Number Species", main=paste("Removing Reads",i,"And Below",sep=" "))
	dev.off()

	# rarefaction curve to check the effect of sampling depth on species found
	# without labels
	pdf(paste0(outputdir,"IndoRarefactionCurve_100K_NoSpeciesFiltering_RemovingCounts_",i,"_andBelow_Counts_noLabels.pdf"), width=10)
	rarecurve(t(otu_table(Indonesian, taxa_are_rows = TRUE)), step=20, label=F, xlab="Counts",ylab="Number Species", main=paste("Removing Reads",i,"And Below",sep=" "))
	dev.off()
}

####################
### Control Data ###
####################

# load in 100K control phyloseq object
load(paste0(inputdir,"control_100K_Counts_physeq.Rda"))

# Separate species' abundances and taxonomy columns
control <- otu_table(control_100K_Counts_physeq)

# Try out different filtering methods:
for (i in 1:5){
	# remove OTUs at or under thresholding of 1:5
	# get rid of singletons
	control[control<=i]<-0

	# rarefaction curve to check the effect of sampling depth on species found
	# here, we'll do it with labels. Note, you need to make the DF into an otu table, since it's
	# a phyloseq object
	pdf(paste0(outputdir,"controlRarefactionCurve_100K_NoSpeciesFiltering_RemovingCounts_",i,"_andBelow_Counts.pdf"), width=10)
	rarecurve(t(otu_table(control, taxa_are_rows = TRUE)), step=20, label=T, xlab="Counts",ylab="Number Species",main=paste("Removing Reads",i,"And Below",sep=" "))
	dev.off()

	# rarefaction curve to check the effect of sampling depth on species found
	# without labels
	pdf(paste0(outputdir,"controlRarefactionCurve_100K_NoSpeciesFiltering_RemovingCounts_",i,"_andBelow_Counts_noLabels.pdf"), width=10)
	rarecurve(t(otu_table(control, taxa_are_rows = TRUE)), step=20, label=F, xlab="Counts",ylab="Number Species",main=paste("Removing Reads",i,"And Below",sep=" "))
	dev.off()
}

#########################
###### Mantel Test ######
#########################

# For Indonesian data -------------------------------
Indonesian <- otu_table(subSampled_100K_Indo_Counts_physeq)
# copy the all reads dataframe, as you will be altering it and need to compare to the original 
Indonesian_noFilter=Indonesian

# make an empty matrix to save data to
mantel_result_matrix=matrix(nrow=5,ncol=2)
colnames(mantel_result_matrix)=c("Spearman_Rho","Significance")
rownames(mantel_result_matrix)=c(1:5)
for (i in 1:5){
	# get rid of anything with 5 counts or less
	Indonesian[Indonesian<=i]<-0
	dist.allReads=vegdist(Indonesian, method = "bray")
	dist.allReads.Nofilter=vegdist(Indonesian_noFilter, method = "bray")
	mantel_result=mantel(dist.allReads, dist.allReads.Nofilter, method = "spearman", permutations = 99, na.rm = TRUE)
	mantel_result_matrix[i,]=c(mantel_result$statistic, mantel_result$signif)
}
mantel_result_matrix <- cbind(mantel_result_matrix, "nFiltered"=1:nrow(mantel_result_matrix)) 
write.table(mantel_result_matrix,file=paste0(outputdir,"Indonesian100K_Mantel_result_matrix_BraysDissimilarity_SpearmanRho.txt"),sep="\t")

# For control data -------------------------------
control <- otu_table(control_100K_Counts_physeq)
# copy the all reads dataframe, as you will be altering it and need to compare to the original 
control_noFilter=control

# make an empty matrix to save data to
mantel_result_matrix=matrix(nrow=5,ncol=2)
colnames(mantel_result_matrix)=c("Spearman_Rho","Significance")
rownames(mantel_result_matrix)=c(1:5)
for (i in 1:5){
	# get rid of anything with 5 counts or less
	control[control<=i]<-0
	dist.allReads=vegdist(control, method = "bray")
	dist.allReads.Nofilter=vegdist(control_noFilter, method = "bray")
	mantel_result=mantel(dist.allReads, dist.allReads.Nofilter, method = "spearman", permutations = 99, na.rm = TRUE)
	mantel_result_matrix[i,]=c(mantel_result$statistic, mantel_result$signif)
}
mantel_result_matrix <- cbind(mantel_result_matrix, "nFiltered"=1:nrow(mantel_result_matrix)) 
write.table(mantel_result_matrix,file=paste0(outputdir,"Control100K_Mantel_result_matrix_BraysDissimilarity_SpearmanRho.txt"),sep="\t")

