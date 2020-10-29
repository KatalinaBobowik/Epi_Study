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
refdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/EpiStudy/"
outputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/Indo_250K/"


# Load in count data
AllREadsSE_Indo_Counts <- read.csv(paste0(refdir,"Counts_Indo_AllReads_SE.csv"),check.names=FALSE)

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

# get information on phyloseq objects
AllREadsSE_Indo_Counts_physeq

# save phyloseq object
save(AllREadsSE_Indo_Counts_physeq, file = paste0(outputdir, "AllREadsSE_Indo_Counts_physeq.Rda"))

