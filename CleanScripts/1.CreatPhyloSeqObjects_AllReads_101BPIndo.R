# created by KSB, 20.01.19
# Create a PhyloSeq object for the sub-sampled, 250K Indonesian data

# load packages
require(ggplot2)
require(RColorBrewer)
library(dplyr)
library(plyr)
library(phyloseq)
library(reshape2)
library(ggpubr)
library(vegan)

# set up directories
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/EpiStudy/"
outputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/All_Reads/"

####################
##### RPM Data #####
####################

# Load in data
allReads_RPM <- read.csv(paste0(inputdir,"Indo_AllReads_NoFilteringHumanIncluded_RPM.csv"),check.names=FALSE)

# Remove Metazoa and Viridiplantae
# allReads_RPM=allReads_RPM[-which(allReads_RPM$Kingdom %in% c("Viridiplantae","Metazoa")),]

# Convert to PhyloSeq object ------------------------------------------

# Separate species' abundances and taxonomy columns
taxa_raw <- as.matrix(allReads_RPM[,c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species")])
abund_raw <- as.matrix(allReads_RPM[,-which(colnames(allReads_RPM) %in% c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species"))])

# convert to Phyloseq object
tax = tax_table(taxa_raw)
taxa = otu_table(abund_raw, taxa_are_rows = TRUE)
allReads_RPM_physeq = phyloseq(taxa, tax)
# add in sample information, starting with Island
samplenames <- colnames(otu_table(allReads_RPM_physeq))
island <- sapply(strsplit(samplenames, "[-.]"), `[`, 1)
# make this into a df and add to the Phloseq object
samples_df=data.frame(SampleName=colnames(otu_table(allReads_RPM_physeq)), SamplePop=island)
samples = sample_data(samples_df)
rownames(samples)=samples$SampleName
sample_data(allReads_RPM_physeq) <- samples

# save phyloseq object
save(allReads_RPM_physeq, file = paste0(outputdir, "allReads_RPM_physeq.Rda"))

####################
#### Count Data ####
####################

# Load in data
allReads_Counts <- read.csv(paste0(inputdir,"Indo_AllReads_Counts.csv"),check.names=FALSE)

# Remove Metazoa and Viridiplantae
# allReads_Counts=allReads_Counts[-which(allReads_Counts$Kingdom %in% c("Viridiplantae","Metazoa")),]

# Convert to PhyloSeq object ------------------------------------------

# Separate species' abundances and taxonomy columns
taxa_raw <- as.matrix(allReads_Counts[,c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species")])
abund_raw <- as.matrix(allReads_Counts[,-which(colnames(allReads_Counts) %in% c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species"))])

# convert to Phyloseq object
tax = tax_table(taxa_raw)
taxa = otu_table(abund_raw, taxa_are_rows = TRUE)
allReads_Counts_physeq = phyloseq(taxa, tax)
# add in sample information, starting with Island
samplenames <- colnames(otu_table(allReads_Counts_physeq))
island <- sapply(strsplit(samplenames, "[-.]"), `[`, 1)
# make this into a df and add to the Phloseq object
samples_df=data.frame(SampleName=colnames(otu_table(allReads_Counts_physeq)), SamplePop=island)
samples = sample_data(samples_df)
rownames(samples)=samples$SampleName
sample_data(allReads_Counts_physeq) <- samples

# make sure counts and RPM phyloseq objects are the same
allReads_RPM_physeq
# otu_table()   OTU Table:         [ 151 taxa and 123 samples ]
# sample_data() Sample Data:       [ 123 samples by 2 sample variables ]
# tax_table()   Taxonomy Table:    [ 151 taxa by 8 taxonomic ranks ]
allReads_Counts_physeq
# otu_table()   OTU Table:         [ 151 taxa and 123 samples ]
# sample_data() Sample Data:       [ 123 samples by 2 sample variables ]
# tax_table()   Taxonomy Table:    [ 151 taxa by 8 taxonomic ranks ]

# save phyloseq object
save(allReads_Counts_physeq, file = paste0(outputdir, "allReads_Counts_physeq.Rda"))

