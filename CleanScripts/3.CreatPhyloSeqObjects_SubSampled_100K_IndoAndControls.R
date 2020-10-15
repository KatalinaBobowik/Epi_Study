# created by KSB, 20.01.19
# Create a PhyloSeq object for the sub-sampled, 100K Indonesian data

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
controldir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/Controls/"

#########################
##### Indo RPM Data #####
#########################

# Load in data
subSampled_100K_Indo_RPM <- read.csv(paste0(inputdir,"RPM_Indo_SubSampled_100K.csv"),check.names=FALSE)

# Remove Metazoa and Viridiplantae
# subSampled_100K_Indo_RPM=subSampled_100K_Indo_RPM[-which(subSampled_100K_Indo_RPM$Kingdom %in% c("Viridiplantae","Metazoa")),]

# Convert to PhyloSeq object ------------------------------------------

# Separate species' abundances and taxonomy columns
taxa_raw <- as.matrix(subSampled_100K_Indo_RPM[,c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species")])
abund_raw <- as.matrix(subSampled_100K_Indo_RPM[,-which(colnames(subSampled_100K_Indo_RPM) %in% c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species"))])

# convert to Phyloseq object
tax = tax_table(taxa_raw)
taxa = otu_table(abund_raw, taxa_are_rows = TRUE)
subSampled_100K_Indo_RPM_physeq = phyloseq(taxa, tax)
# add in sample information, starting with Island
samplenames <- colnames(otu_table(subSampled_100K_Indo_RPM_physeq))
pop <- rep("Indonesia",ncol(otu_table(subSampled_100K_Indo_RPM_physeq)))
# make this into a df and add to the Phloseq object
samples_df=data.frame(SampleName=colnames(otu_table(subSampled_100K_Indo_RPM_physeq)), SamplePop=pop)
samples = sample_data(samples_df)
rownames(samples)=samples$SampleName
sample_data(subSampled_100K_Indo_RPM_physeq) <- samples

# save phyloseq object
save(subSampled_100K_Indo_RPM_physeq, file = paste0(outputdir, "subSampled_100K_Indo_RPM_physeq.Rda"))

#########################
#### Indo Count Data ####
#########################

# Load in data
subSampled_100K_Indo_Counts <- read.csv(paste0(inputdir,"Indo_noFiltering_Subsampled_Counts_100K.csv"),check.names=FALSE)

# Remove Metazoa and Viridiplantae
# subSampled_100K_Indo_Counts=subSampled_100K_Indo_Counts[-which(subSampled_100K_Indo_Counts$Kingdom %in% c("Viridiplantae","Metazoa")),]

# Convert to PhyloSeq object ------------------------------------------

# Separate species' abundances and taxonomy columns
taxa_raw <- as.matrix(subSampled_100K_Indo_Counts[,c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species")])
abund_raw <- as.matrix(subSampled_100K_Indo_Counts[,-which(colnames(subSampled_100K_Indo_Counts) %in% c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species"))])

# convert to Phyloseq object
tax = tax_table(taxa_raw)
taxa = otu_table(abund_raw, taxa_are_rows = TRUE)
subSampled_100K_Indo_Counts_physeq = phyloseq(taxa, tax)
# add in sample information, starting with Island
samplenames <- colnames(otu_table(subSampled_100K_Indo_Counts_physeq))
pop <- rep("Indonesia",ncol(otu_table(subSampled_100K_Indo_Counts_physeq)))
# make this into a df and add to the Phloseq object
samples_df=data.frame(SampleName=colnames(otu_table(subSampled_100K_Indo_Counts_physeq)), SamplePop=pop)
samples = sample_data(samples_df)
rownames(samples)=samples$SampleName
sample_data(subSampled_100K_Indo_Counts_physeq) <- samples

# make sure counts and RPM phyloseq objects are the same
subSampled_100K_Indo_RPM_physeq
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 1496 taxa and 123 samples ]
# sample_data() Sample Data:       [ 123 samples by 2 sample variables ]
# tax_table()   Taxonomy Table:    [ 1496 taxa by 8 taxonomic ranks ]
subSampled_100K_Indo_Counts_physeq
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 1496 taxa and 123 samples ]
# sample_data() Sample Data:       [ 123 samples by 2 sample variables ]
# tax_table()   Taxonomy Table:    [ 1496 taxa by 8 taxonomic ranks ]

# save phyloseq object
save(subSampled_100K_Indo_Counts_physeq, file = paste0(outputdir, "subSampled_100K_Indo_Counts_physeq.Rda"))

############################
##### Control RPM Data #####
############################

# Load in data
control_100K_RPM <- read.csv(paste0(inputdir,"Control_100K_NoFiltering_RPM.csv"),check.names=FALSE)

# Remove Metazoa and Viridiplantae
# control_100K_RPM=control_100K_RPM[-which(control_100K_RPM$Kingdom %in% c("Viridiplantae","Metazoa")),]

# Convert to PhyloSeq object ------------------------------------------

# Separate species' abundances and taxonomy columns
taxa_raw <- as.matrix(control_100K_RPM[,c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species")])
abund_raw <- as.matrix(control_100K_RPM[,-which(colnames(control_100K_RPM) %in% c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species"))])

# convert to Phyloseq object
tax = tax_table(taxa_raw)
taxa = otu_table(abund_raw, taxa_are_rows = TRUE)
control_100K_RPM_physeq = phyloseq(taxa, tax)

# only get control samples
diseaseStatus=read.table(paste0(inputdir,"ControlSamples_DiseaseStatus.txt"))
colnames(diseaseStatus)=c("Samples","diseaseStatus")
controlSamples=as.character(diseaseStatus[grep("CSF", diseaseStatus$diseaseStatus),"Samples"])
# prune out samples we don't want
control_100K_RPM_physeq=prune_samples(controlSamples,control_100K_RPM_physeq)

# add in sample information, starting with Island
samplenames <- colnames(otu_table(control_100K_RPM_physeq))
pop <- rep("Netherlands",ncol(otu_table(control_100K_RPM_physeq)))
# make this into a df and add to the Phloseq object
samples_df=data.frame(SampleName=colnames(otu_table(control_100K_RPM_physeq)), SamplePop=pop)
samples = sample_data(samples_df)
rownames(samples)=samples$SampleName
sample_data(control_100K_RPM_physeq) <- samples

# save phyloseq object
save(control_100K_RPM_physeq, file = paste0(outputdir, "control_100K_RPM_physeq.Rda"))

###############################
##### Control Counts Data #####
###############################

# Load in data
control_100K_Counts <- read.csv(paste0(inputdir,"Controls_NoFiltering_Counts.csv"),check.names=FALSE)

# Remove Metazoa and Viridiplantae
# control_100K_Counts=control_100K_Counts[-which(control_100K_Counts$Kingdom %in% c("Viridiplantae","Metazoa")),]

# Convert to PhyloSeq object ------------------------------------------

# Separate species' abundances and taxonomy columns
taxa_raw <- as.matrix(control_100K_Counts[,c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species")])
abund_raw <- as.matrix(control_100K_Counts[,-which(colnames(control_100K_Counts) %in% c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species"))])

# convert to Phyloseq object
tax = tax_table(taxa_raw)
taxa = otu_table(abund_raw, taxa_are_rows = TRUE)
control_100K_Counts_physeq = phyloseq(taxa, tax)

# prune out samples we don't want
control_100K_Counts_physeq=prune_samples(controlSamples,control_100K_Counts_physeq)

# add in sample information, starting with Island
samplenames <- colnames(otu_table(control_100K_Counts_physeq))
pop <- rep("Netherlands",ncol(otu_table(control_100K_Counts_physeq)))
# make this into a df and add to the Phloseq object
samples_df=data.frame(SampleName=colnames(otu_table(control_100K_Counts_physeq)), SamplePop=pop)
samples = sample_data(samples_df)
rownames(samples)=samples$SampleName
sample_data(control_100K_Counts_physeq) <- samples

# make sure counts and RPM phyloseq objects are the same
control_100K_RPM_physeq
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 971 taxa and 49 samples ]
# sample_data() Sample Data:       [ 49 samples by 2 sample variables ]
# tax_table()   Taxonomy Table:    [ 971 taxa by 8 taxonomic ranks ]
control_100K_Counts_physeq
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 971 taxa and 49 samples ]
# sample_data() Sample Data:       [ 49 samples by 2 sample variables ]
# tax_table()   Taxonomy Table:    [ 971 taxa by 8 taxonomic ranks ]

# save phyloseq object
save(control_100K_Counts_physeq, file = paste0(outputdir, "control_100K_Counts_physeq.Rda"))



