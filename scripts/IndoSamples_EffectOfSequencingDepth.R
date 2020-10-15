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
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/EpiStudy/"
outputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/"

# no filtering of data (i.e., human reads remain)
raw_CCMetagen_data <- read.csv(paste0(inputdir,"Indo_AllReads_Counts.csv"),check.names=FALSE)

# Data preprocessing ------------------------------------

# add a new column containing family names and superkingdom
raw_CCMetagen_data$SuperKFamily <- paste(raw_CCMetagen_data$Superkingdom, raw_CCMetagen_data$Species, sep="_")

# Separate species' abundances and taxonomy columns
abund_raw_noSingletons <- as.matrix(raw_CCMetagen_data[,c(1:123)])

# get rid of anything with 5 counts or less
abund_raw_noSingletons[abund_raw_noSingletons==1]<-0
abund_raw_noSingletons[abund_raw_noSingletons==2]<-0
abund_raw_noSingletons[abund_raw_noSingletons==3]<-0
abund_raw_noSingletons[abund_raw_noSingletons==4]<-0
abund_raw_noSingletons[abund_raw_noSingletons==5]<-0

# rarefaction curve to check the effect of sampling depth on species found
# with labels
pdf(paste0(outputdir,"rarefactionCurve_IndoSamples_NoSpeciesFiltering_RemovingCounts5andBelow_Counts.pdf"), width=10)
rarecurve(t(otu_table(abund_raw_noSingletons, taxa_are_rows = TRUE)), step=20, label=T, xlab="Counts",ylab="Number Species")
dev.off()

# rarefaction curve to check the effect of sampling depth on species found
# without labels
pdf(paste0(outputdir,"rarefactionCurve_IndoSamples_NoSpeciesFiltering_RemovingCounts5andBelow_Counts_noLabels.pdf"), width=10)
rarecurve(t(otu_table(abund_raw_noSingletons, taxa_are_rows = TRUE)), step=20, label=F, xlab="Counts",ylab="Number Species")
dev.off()

# rarefaction curves in subsampled data
# now for subsampled data
SubSampled_noFilter_counts =  read.csv(paste0(inputdir,"Indo_noFiltering_Subsampled_Counts.csv"),check.names=FALSE)
# add a new column containing family names and superkingdom
SubSampled_noFilter_counts$SuperKFamily <- paste(SubSampled_noFilter_counts$Superkingdom, SubSampled_noFilter_counts$Species, sep="_")

# Separate species' abundances and taxonomy columns
abund_raw_subsampled_noSingletons <- as.matrix(SubSampled_noFilter_counts[,c(1:123)])

# get rid of anything with 5 counts or less
abund_raw_subsampled_noSingletons[abund_raw_subsampled_noSingletons==1]<-0
abund_raw_subsampled_noSingletons[abund_raw_subsampled_noSingletons==2]<-0
abund_raw_subsampled_noSingletons[abund_raw_subsampled_noSingletons==3]<-0
abund_raw_subsampled_noSingletons[abund_raw_subsampled_noSingletons==4]<-0
abund_raw_subsampled_noSingletons[abund_raw_subsampled_noSingletons==5]<-0

# rarefaction curve to check the effect of sampling depth on species found
# with labels
pdf(paste0(outputdir,"rarefactionCurve_IndoSamplesSubSampled_NoSpeciesFiltering_RemovingCounts5andBelow_Counts.pdf"), width=10)
rarecurve(t(otu_table(abund_raw_subsampled_noSingletons, taxa_are_rows = TRUE)), step=20, label=T, xlab="Counts",ylab="Number Species")
dev.off()

# rarefaction curve to check the effect of sampling depth on species found
# without labels
pdf(paste0(outputdir,"rarefactionCurve_IndoSamplesSubSampled_NoSpeciesFiltering_RemovingCounts5andBelow_Counts_noLabels.pdf"), width=10)
rarecurve(t(otu_table(abund_raw_subsampled_noSingletons, taxa_are_rows = TRUE)), step=20, label=F, xlab="Counts",ylab="Number Species")
dev.off()

# Compare subsampled data to non-subsampled data --------------------

# first, do this without filtering singletons
abund_raw <- as.matrix(raw_CCMetagen_data[,c(1:123)])
abund_raw=abund_raw!=0

abund_raw_subsampled=as.matrix(SubSampled_noFilter_counts[,c(1:123)])
abund_raw_subsampled=abund_raw_subsampled!=0
subsampleMelted=melt(abund_raw_subsampled)
colnames(subsampleMelted)[1]="ID"
subsampleMelted[1]="Subsampled"
meltedAllReads=melt(abund_raw)
colnames(meltedAllReads)[1]="ID"
meltedAllReads[1]="AllReads"
allReadsVsSubsample=rbind(subsampleMelted,meltedAllReads)

# plot side by side 
ggplot(allReadsVsSubsample, aes(fill=ID, y=as.numeric(value), x=Var2)) + 
    stat_summary(fun = sum, geom = "bar", position = "fill")

# now do this when removing reads 5 and less ------------------------
abund_raw_noSingletons=abund_raw_noSingletons!=0
abund_raw_subsampled_noSingletons=abund_raw_subsampled_noSingletons!=0

abund_raw_subsampled_SingletonsRemoved=melt(abund_raw_subsampled_noSingletons)
colnames(abund_raw_subsampled_SingletonsRemoved)[1]="ID"
abund_raw_subsampled_SingletonsRemoved[1]="Subsampled"
allReads_SingletonsRemoved=melt(abund_raw_noSingletons)
colnames(allReads_SingletonsRemoved)[1]="ID"
allReads_SingletonsRemoved[1]="AllReads"
allReadsVsSubsample=rbind(abund_raw_subsampled_SingletonsRemoved,allReads_SingletonsRemoved)

# plot side by side
ggplot(allReadsVsSubsample, aes(fill=ID, y=as.numeric(value), x=Var2)) + 
    stat_summary(fun = sum, geom = "bar", position = "fill", stat="identity") +
    theme(axis.text.x = element_text(angle = 90))

# plot depth by species found ------------------

depth <- as.matrix(raw_CCMetagen_data[,c(1:123)])
depth=colSums(depth)
nSpecies=abund_raw2!=0
nSpecies=colSums(nSpecies)
speciesSummary=data.frame(depth,nSpecies)
speciesSummary$librarySize="AllReads"
speciesSummary$samples=rownames(speciesSummary)

# for sub-sampled data
depth_subSampled=as.matrix(SubSampled_noFilter_counts[,c(1:123)])
depth_subSampled=colSums(depth_subSampled)
nSpecies_subSampled=abund_raw_subsampled!=0
nSpecies_subSampled=colSums(nSpecies_subSampled)
speciesSummary_subSampled=data.frame(depth_subSampled,nSpecies_subSampled)
speciesSummary_subSampled$librarySize="SubSampled"
colnames(speciesSummary_subSampled)=c("depth","nSpecies","librarySize")
speciesSummary_subSampled$samples=rownames(speciesSummary_subSampled)

# bind dataframes together
summaryTable=rbind(speciesSummary,speciesSummary_subSampled)

# exploratory plots 

# correlation of sequencing depth on number of species found
ggplot(data=summaryTable, aes(x=depth, y=nSpecies,group=librarySize)) +
  geom_point(aes(color=librarySize)) + geom_smooth(method="lm")

ggplot(data=summaryTable, aes(x=depth, y=nSpecies,fill=librarySize)) +
	geom_bar(position="dodge", stat="identity", width = 500)

ggplot(data=summaryTable, aes(x=depth, fill=librarySize)) +
	geom_density(alpha=0.4)


# ----------------

rarefy_obs(abund_raw_noSingletons, "otu_counts", other_cols = TRUE)

ps.rarefied = rarefy_even_depth(otu_table(abund_raw_noSingletons), rngseed=1, sample.size=500, replace=F)

