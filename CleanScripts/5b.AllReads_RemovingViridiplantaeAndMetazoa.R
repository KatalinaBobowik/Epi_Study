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
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/Indo_250K/"
outputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/Indo_250K/"

# load in allReads Indonesian phyloseq counts object
load(paste0(inputdir,"subSampled_250K_Indo_Counts_physeq.Rda"))
# load in allreads Indonesian phyloseq RPM object
load(paste0(inputdir,"subSampled_250K_Indo_RPM_physeq.Rda"))

# filter reads 1 and below
otu_table(subSampled_250K_Indo_Counts_physeq)[otu_table(subSampled_250K_Indo_Counts_physeq)<=1]<-0
otu_table(subSampled_250K_Indo_RPM_physeq)[otu_table(subSampled_250K_Indo_RPM_physeq)<=1]<-0

# Get a Summary of the read depth before any filtering
SeqDepth_preFilter = colSums(otu_table(subSampled_250K_Indo_Counts_physeq))
sample_data(subSampled_250K_Indo_Counts_physeq)$SeqDepth_preFilter = SeqDepth_preFilter

# investigate where reads are mapping to
p=plot_bar(subSampled_250K_Indo_RPM_physeq, fill = "Kingdom")
fig <- p + geom_bar(stat="identity", position="stack")
pdf(paste0(outputdir,"AllReads_IncludingHumanandPlants_KingdomLevel_RPM.pdf"), width=10)
fig
dev.off()

# Some human reads are getting through. We also see a low level of plants
# from Streptophya sequences
Viridiplantae=subset_taxa(subSampled_250K_Indo_RPM_physeq, (Kingdom=="Viridiplantae"))
p=plot_bar(Viridiplantae, fill = "Phylum")
fig <- p + geom_bar(stat="identity", position="stack")
pdf(paste0(outputdir,"Plants_RPM.pdf"), width=10)
fig
dev.off()

# Remove Viridiplantae and get number of sequences after removal
subSampled_250K_Indo_Counts_physeq=subset_taxa(subSampled_250K_Indo_Counts_physeq, (Kingdom!="Viridiplantae"))
subSampled_250K_Indo_RPM_physeq=subset_taxa(subSampled_250K_Indo_RPM_physeq, (Kingdom!="Viridiplantae"))

# Get a Summary of the read depth after filtering
SeqDepth_noPlants = colSums(otu_table(subSampled_250K_Indo_Counts_physeq))
sample_data(subSampled_250K_Indo_Counts_physeq)$SeqDepth_noPlants = SeqDepth_noPlants

# Now let's investigate the Metazoa reads
Metazoa=subset_taxa(subSampled_250K_Indo_RPM_physeq, (Kingdom=="Metazoa"))
p=plot_bar(Metazoa, fill = "Phylum")
fig <- p + geom_bar(stat="identity", position="stack")
pdf(paste0(outputdir,"Metazoa_RPM.pdf"), width=10)
fig
dev.off()

# we see that most of the metazoa reads aru humans that got through.
# Let's remove them and plot again
Metazoa=subset_taxa(Metazoa, (Phylum!="Chordata"))
# filter out any 0's that remain in the dataframe
Metazoa <- prune_taxa(taxa_sums(Metazoa) > 0, Metazoa)
p=plot_bar(Metazoa, fill = "Phylum")
fig <- p + geom_bar(stat="identity", position="stack")
pdf(paste0(outputdir,"Metazoa_NoHumans_RPM.pdf"), width=10)
fig
dev.off()

# We can see that arthropods (probably caused by mismapping), nematodes,
# and platyhelminths are found in the data. 

# Blasted sequences for Spirometra erinaceieuropaei (found in MTW-MDB-015), and Caenorhabditis remanei,
# found in MTW-TLL-006. These sequences both have reads that map equally well to human
# sequences. Thus, I will remove these in particular since helminths are known pathogens
# and I don't want the results to look misleading.

# Remove Metazoa and get number of sequences after removal
subSampled_250K_Indo_Counts_physeq=subset_taxa(subSampled_250K_Indo_Counts_physeq, (Kingdom!="Metazoa"))
subSampled_250K_Indo_RPM_physeq=subset_taxa(subSampled_250K_Indo_RPM_physeq, (Kingdom!="Metazoa"))

# Get a Summary of the read depth after filtering
SeqDepth_noAnimals = colSums(otu_table(subSampled_250K_Indo_Counts_physeq))
sample_data(subSampled_250K_Indo_Counts_physeq)$SeqDepth_noAnimals = SeqDepth_noAnimals

# save table
write.table(sample_data(subSampled_250K_Indo_Counts_physeq)[,c("SeqDepth_preFilter","SeqDepth_noPlants","SeqDepth_noAnimals")], file=paste0(outputdir, "taxaFilteringSummary.txt"))
