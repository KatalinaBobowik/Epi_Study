# created by KSB, 20.01.19
# Script created to inspect where reads are mapping to. 

# load packages
require(ggplot2)
require(RColorBrewer)
library(dplyr)
library(plyr)
library(reshape2)
library(ggpubr)
library(vegan)
library(ggfortify)
library(microbiome)
library(microbiomeutilities)
library(viridis)
library(tibble)
library(knitr)

# set ggplot colour theme to white
theme_set(theme_bw())

# set up directories
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/All_Reads/"
outputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/All_Reads/"

# load in allReads Indonesian phyloseq counts object
load(paste0(inputdir,"allReads_Counts_physeq.Rda"))
# load in allreads Indonesian phyloseq RPM object
load(paste0(inputdir,"allReads_RPM_physeq.Rda"))

# filter reads 1 and below
otu_table(allReads_Counts_physeq)[otu_table(allReads_Counts_physeq)<=1]<-0
otu_table(allReads_RPM_physeq)[otu_table(allReads_RPM_physeq)<=1]<-0

# Get a Summary of the read depth before any filtering
SeqDepth_preFilter = colSums(otu_table(allReads_Counts_physeq))
sample_data(allReads_Counts_physeq)$SeqDepth_preFilter = SeqDepth_preFilter

# investigate where reads are mapping to
p=plot_bar(allReads_Counts_physeq, fill = "Phylum")
fig <- p + geom_bar(stat="identity", position="stack")
fig

# Some human reads are getting through. We also see a low level of plants
# from Streptophya sequences
Viridiplantae=subset_taxa(allReads_Counts_physeq, (Kingdom=="Viridiplantae"))
p=plot_bar(Viridiplantae, fill = "Phylum")
fig <- p + geom_bar(stat="identity", position="stack")
fig

# Remove Viridiplantae and get number of sequences after removal
allReads_Counts_physeq=subset_taxa(allReads_Counts_physeq, (Kingdom!="Viridiplantae"))
allReads_RPM_physeq=subset_taxa(allReads_RPM_physeq, (Kingdom!="Viridiplantae"))

# Get a Summary of the read depth after filtering
SeqDepth_noPlants = colSums(otu_table(allReads_Counts_physeq))
sample_data(allReads_Counts_physeq)$SeqDepth_noPlants = SeqDepth_noPlants

# Now let's investigate the Metazoa reads
Metazoa=subset_taxa(allReads_Counts_physeq, (Kingdom=="Metazoa"))
p=plot_bar(Metazoa, fill = "Phylum")
fig <- p + geom_bar(stat="identity", position="stack")
fig

# we see that most of the metazoa reads aru humans that got through.
# Let's remove them and plot again
Metazoa=subset_taxa(Metazoa, (Phylum!="Chordata"))
# filter out any 0's that remain in the dataframe
Metazoa <- prune_taxa(taxa_sums(Metazoa) > 0, Metazoa)
p=plot_bar(Metazoa, fill = "Phylum")
fig <- p + geom_bar(stat="identity", position="stack")
fig

# We can see that arthropods (probably caused by mismapping), nematodes,
# and platyhelminths are found in the data. 

# Blasted sequences for Spirometra erinaceieuropaei (found in MTW-MDB-015), and Caenorhabditis remanei,
# found in MTW-TLL-006. These sequences both have reads that map equally well to human
# sequences. Thus, I will remove these in particular since helminths are known pathogens
# and I don't want the results to look misleading.

# Remove Metazoa and get number of sequences after removal
allReads_Counts_physeq=subset_taxa(allReads_Counts_physeq, (Kingdom!="Metazoa"))
allReads_RPM_physeq=subset_taxa(allReads_RPM_physeq, (Kingdom!="Metazoa"))

# Get a Summary of the read depth after filtering
SeqDepth_noAnimals = colSums(otu_table(allReads_Counts_physeq))
sample_data(allReads_Counts_physeq)$SeqDepth_noAnimals = SeqDepth_noAnimals

