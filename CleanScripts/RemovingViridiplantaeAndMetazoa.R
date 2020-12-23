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
inputdir_100BP = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/SEIndonesianSamples/"
inputdir75BP = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/ControlSampleComparison_TB/"
outputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/SEIndonesianSamples/"

# load in SE Indonesian phyloseq counts object
load(paste0(inputdir_100BP,"AllREadsSE_Indo_Counts_physeq_noSingletons.Rda"))

# investigate where reads are mapping to
p=plot_bar(AllREadsSE_Indo_Counts_physeq, fill = "Kingdom")
fig <- p + geom_bar(stat="identity", position="stack")
fig

# Some human reads are getting through. We also see a low level of plants
# from Streptophya sequences
Viridiplantae=subset_taxa(AllREadsSE_Indo_Counts_physeq, (Kingdom=="Viridiplantae"))
p=plot_bar(Viridiplantae, fill = "Phylum")
fig <- p + geom_bar(stat="identity", position="stack")
fig

# Remove Viridiplantae and get number of sequences after removal
AllREadsSE_Indo_Counts_physeq=subset_taxa(AllREadsSE_Indo_Counts_physeq, (Kingdom!="Viridiplantae"))
AllREadsSE_Indo_Counts_physeq=subset_taxa(AllREadsSE_Indo_Counts_physeq, (Kingdom!="Viridiplantae"))

# Get a Summary of the read depth after filtering
SeqDepth_noPlants = colSums(otu_table(AllREadsSE_Indo_Counts_physeq))
sample_data(AllREadsSE_Indo_Counts_physeq)$SeqDepth_noPlants = SeqDepth_noPlants

# Now let's investigate the Metazoa reads
Metazoa=subset_taxa(AllREadsSE_Indo_Counts_physeq, (Kingdom=="Metazoa"))
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

##############
# 75 BP data # 
##############

# load in SE phyloseq counts object from all pops
load(paste0(inputdir75BP,"merged_phylo_counts_noSingletons.Rda"))

# investigate where reads are mapping to
p=plot_bar(merged_phylo_counts, fill = "Phylum")
fig <- p + geom_bar(stat="identity", position="stack")
fig

# Some human reads are getting through. We also see a low level of plants
# from Streptophya sequences
Viridiplantae=subset_taxa(merged_phylo_counts, (Kingdom=="Viridiplantae"))
p=plot_bar(Viridiplantae, fill = "Phylum")
fig <- p + geom_bar(stat="identity", position="stack") + facet_wrap(~ SamplePop, scales = "free")
fig

# Remove Viridiplantae and get number of sequences after removal
merged_phylo_counts=subset_taxa(merged_phylo_counts, (Kingdom!="Viridiplantae"))

# Now let's investigate the Metazoa reads
Metazoa=subset_taxa(merged_phylo_counts, (Kingdom=="Metazoa"))
p=plot_bar(Metazoa, fill = "Phylum")
fig <- p + geom_bar(stat="identity", position="stack") + facet_wrap(~ SamplePop, scales = "free")
fig

# we see that most of the metazoa reads aru humans that got through.
# Let's remove them and plot again
Metazoa=subset_taxa(Metazoa, (Phylum!="Chordata"))
# filter out any 0's that remain in the dataframe
Metazoa <- prune_taxa(taxa_sums(Metazoa) > 0, Metazoa)
p=plot_bar(Metazoa, fill = "Phylum")
fig <- p + geom_bar(stat="identity", position="stack") + facet_wrap(~ SamplePop, scales = "free")
fig

Metazoa=subset_taxa(Metazoa, (Phylum!="Mollusca"))
# filter out any 0's that remain in the dataframe
Metazoa <- prune_taxa(taxa_sums(Metazoa) > 0, Metazoa)
p=plot_bar(Metazoa, fill = "Phylum")
fig <- p + geom_bar(stat="identity", position="stack") + facet_wrap(~ SamplePop, scales = "free")
fig

Metazoa=subset_taxa(Metazoa, (Phylum!="Arthropoda"))
# filter out any 0's that remain in the dataframe
Metazoa <- prune_taxa(taxa_sums(Metazoa) > 0, Metazoa)
p=plot_bar(Metazoa, fill = "Phylum")
fig <- p + geom_bar(stat="identity", position="stack") + facet_wrap(~ SamplePop, scales = "free")
fig



