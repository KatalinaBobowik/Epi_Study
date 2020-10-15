# created by KSB, 20.01.19
# Create a PhyloSeq object for the sub-sampled, AllREadsSE Indonesian data

# load packages
require(ggplot2)
require(RColorBrewer)
library(dplyr)
library(plyr)
library(phyloseq)
library(reshape2)
library(ggpubr)
library(vegan)
library(ggfortify)

# set up directories
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/EpiStudy/"
outputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/100KSamples/"

#########################
##### Indo RPM Data #####
#########################

# Load in data
AllREadsSE_Indo_RPM <- read.csv(paste0(inputdir,"RPM_Indo_AllReads_SE.csv"),check.names=FALSE)

# Remove Metazoa and Viridiplantae
# AllREadsSE_Indo_RPM=subSampled_100K_Indo_RPM[-which(AllREadsSE_Indo_RPM$Kingdom %in% c("Viridiplantae","Metazoa")),]

# Convert to PhyloSeq object ------------------------------------------

# Separate species' abundances and taxonomy columns
taxa_raw <- as.matrix(AllREadsSE_Indo_RPM[,c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species")])
abund_raw <- as.matrix(AllREadsSE_Indo_RPM[,-which(colnames(AllREadsSE_Indo_RPM) %in% c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species"))])

# convert to Phyloseq object
tax = tax_table(taxa_raw)
taxa = otu_table(abund_raw, taxa_are_rows = TRUE)
AllREadsSE_Indo_RPM_physeq = phyloseq(taxa, tax)
# add in sample information, starting with Island
samplenames <- colnames(otu_table(AllREadsSE_Indo_RPM_physeq))
pop <- rep("Indonesia",ncol(otu_table(AllREadsSE_Indo_RPM_physeq)))
# make this into a df and add to the Phloseq object
samples_df=data.frame(SampleName=colnames(otu_table(AllREadsSE_Indo_RPM_physeq)), SamplePop=pop)
samples = sample_data(samples_df)
rownames(samples)=samples$SampleName
sample_data(AllREadsSE_Indo_RPM_physeq) <- samples

# save phyloseq object
save(AllREadsSE_Indo_RPM_physeq, file = paste0(outputdir, "AllREadsSE_Indo_RPM_physeq.Rda"))

#########################
#### Indo Count Data ####
#########################

# Load in data
AllREadsSE_Indo_Counts <- read.csv(paste0(inputdir,"Counts_Indo_AllReads_SE.csv"),check.names=FALSE)

# Remove Metazoa and Viridiplantae
# subSampled_AllREadsSE_Indo_Counts=subSampled_AllREadsSE_Indo_Counts[-which(subSampled_AllREadsSE_Indo_Counts$Kingdom %in% c("Viridiplantae","Metazoa")),]

# Convert to PhyloSeq object ------------------------------------------

# Separate species' abundances and taxonomy columns
taxa_raw <- as.matrix(AllREadsSE_Indo_Counts[,c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species")])
abund_raw <- as.matrix(AllREadsSE_Indo_Counts[,-which(colnames(AllREadsSE_Indo_Counts) %in% c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species"))])

# convert to Phyloseq object
tax = tax_table(taxa_raw)
taxa = otu_table(abund_raw, taxa_are_rows = TRUE)
AllREadsSE_Indo_Counts_physeq = phyloseq(taxa, tax)
# add in sample information, starting with Island
samplenames <- colnames(otu_table(AllREadsSE_Indo_Counts_physeq))
pop <- rep("Indonesia",ncol(otu_table(AllREadsSE_Indo_Counts_physeq)))
# make this into a df and add to the Phloseq object
samples_df=data.frame(SampleName=colnames(otu_table(AllREadsSE_Indo_Counts_physeq)), SamplePop=pop)
samples = sample_data(samples_df)
rownames(samples)=samples$SampleName
sample_data(AllREadsSE_Indo_Counts_physeq) <- samples

# make sure counts and RPM phyloseq objects are the same
AllREadsSE_Indo_RPM_physeq
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 2863 taxa and 123 samples ]
# sample_data() Sample Data:       [ 123 samples by 2 sample variables ]
# tax_table()   Taxonomy Table:    [ 2863 taxa by 8 taxonomic ranks ]
AllREadsSE_Indo_Counts_physeq
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 2863 taxa and 123 samples ]
# sample_data() Sample Data:       [ 123 samples by 2 sample variables ]
# tax_table()   Taxonomy Table:    [ 2863 taxa by 8 taxonomic ranks ]

# save phyloseq object
save(AllREadsSE_Indo_Counts_physeq, file = paste0(outputdir, "AllREadsSE_Indo_Counts_physeq.Rda"))

