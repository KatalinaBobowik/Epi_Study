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
refdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/EpiStudy/"

# Read in the data
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

# filter reads 1 and below
otu_table(AllREadsSE_Indo_Counts_physeq)[otu_table(AllREadsSE_Indo_Counts_physeq)<=1]<-0

# Get a Summary of the read depth before any filtering
SeqDepth_preFilter = colSums(otu_table(AllREadsSE_Indo_Counts_physeq))
sample_data(AllREadsSE_Indo_Counts_physeq)$SeqDepth_preFilter = SeqDepth_preFilter

# investigate where reads are mapping to
p=plot_bar(AllREadsSE_Indo_Counts_physeq, fill = "Phylum")
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

# Remove Metazoa and get number of sequences after removal
AllREadsSE_Indo_Counts_physeq=subset_taxa(AllREadsSE_Indo_Counts_physeq, (Kingdom!="Metazoa"))

# Get a Summary of the read depth after filtering
SeqDepth_noAnimals = colSums(otu_table(AllREadsSE_Indo_Counts_physeq))
sample_data(AllREadsSE_Indo_Counts_physeq)$SeqDepth_noAnimals = SeqDepth_noAnimals

#####################
# For Mali Samples #
#####################

Mali_Counts <- read.csv(paste0(refdir,"MaliControls_Counts_NoFilter.csv"),check.names=FALSE)

# Separate species' abundances and taxonomy columns
taxa_raw <- as.matrix(Mali_Counts[,c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species")])
abund_raw <- as.matrix(Mali_Counts[,-which(colnames(Mali_Counts) %in% c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species"))])
# convert to Phyloseq object
tax = tax_table(taxa_raw)
taxa = otu_table(abund_raw, taxa_are_rows = TRUE)
Mali_Counts_physeq = phyloseq(taxa, tax)
# add in sample information, i.e., the sample names and population they're from
samplenames <- colnames(otu_table(Mali_Counts_physeq))
pop <- rep("Mali",ncol(otu_table(Mali_Counts_physeq)))

# make this into a df and add to the Phloseq object
samples_df=data.frame(SampleName=colnames(otu_table(Mali_Counts_physeq)), SamplePop=pop)
samples = sample_data(samples_df)
rownames(samples)=samples$SampleName
sample_data(Mali_Counts_physeq) <- samples

# filter reads 1 and below
otu_table(Mali_Counts_physeq)[otu_table(Mali_Counts_physeq)<=1]<-0

# Get a Summary of the read depth before any filtering
SeqDepth_preFilter = colSums(otu_table(Mali_Counts_physeq))
sample_data(Mali_Counts_physeq)$SeqDepth_preFilter = SeqDepth_preFilter

# investigate where reads are mapping to
p=plot_bar(Mali_Counts_physeq, fill = "Phylum")
fig <- p + geom_bar(stat="identity", position="stack")
fig

# Some human reads are getting through. We also see a low level of plants
# from Streptophya sequences
Viridiplantae=subset_taxa(Mali_Counts_physeq, (Kingdom=="Viridiplantae"))
p=plot_bar(Viridiplantae, fill = "Phylum")
fig <- p + geom_bar(stat="identity", position="stack")
fig

# Remove Viridiplantae and get number of sequences after removal
Mali_Counts_physeq=subset_taxa(Mali_Counts_physeq, (Kingdom!="Viridiplantae"))

# Get a Summary of the read depth after filtering
SeqDepth_noPlants = colSums(otu_table(Mali_Counts_physeq))
sample_data(Mali_Counts_physeq)$SeqDepth_noPlants = SeqDepth_noPlants

# Now let's investigate the Metazoa reads
Metazoa=subset_taxa(Mali_Counts_physeq, (Kingdom=="Metazoa"))
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
Mali_Counts_physeq=subset_taxa(Mali_Counts_physeq, (Kingdom!="Metazoa"))

# Get a Summary of the read depth after filtering
SeqDepth_noAnimals = colSums(otu_table(Mali_Counts_physeq))
sample_data(Mali_Counts_physeq)$SeqDepth_noAnimals = SeqDepth_noAnimals

################
# USA Samples #
################

Americans_Counts <- read.csv(paste0(refdir,"DengueSamples_Counts_NoFiltering.csv"),check.names=FALSE)

# Separate species' abundances and taxonomy columns
taxa_raw <- as.matrix(Americans_Counts[,c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species")])
abund_raw <- as.matrix(Americans_Counts[,-which(colnames(Americans_Counts) %in% c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species"))])
# convert to Phyloseq object
tax = tax_table(taxa_raw)
taxa = otu_table(abund_raw, taxa_are_rows = TRUE)
Americans_Counts_physeq = phyloseq(taxa, tax)

# add in sample information, i.e., the sample names and population they're from
samplenames <- colnames(otu_table(Americans_Counts_physeq))
pop <- rep("USA",ncol(otu_table(Americans_Counts_physeq)))

# make this into a df and add to the Phloseq object
samples_df=data.frame(SampleName=colnames(otu_table(Americans_Counts_physeq)), SamplePop=pop)
samples = sample_data(samples_df)
rownames(samples)=samples$SampleName
sample_data(Americans_Counts_physeq) <- samples

# filter reads 1 and below
otu_table(Americans_Counts_physeq)[otu_table(Americans_Counts_physeq)<=1]<-0

# Get a Summary of the read depth before any filtering
SeqDepth_preFilter = colSums(otu_table(Americans_Counts_physeq))
sample_data(Americans_Counts_physeq)$SeqDepth_preFilter = SeqDepth_preFilter

# investigate where reads are mapping to
p=plot_bar(Americans_Counts_physeq, fill = "Phylum")
fig <- p + geom_bar(stat="identity", position="stack")
fig

# Some human reads are getting through. We also see a low level of plants
# from Streptophya sequences
Viridiplantae=subset_taxa(Americans_Counts_physeq, (Kingdom=="Viridiplantae"))
p=plot_bar(Viridiplantae, fill = "Phylum")
fig <- p + geom_bar(stat="identity", position="stack")
fig

# Remove Viridiplantae and get number of sequences after removal
Americans_Counts_physeq=subset_taxa(Americans_Counts_physeq, (Kingdom!="Viridiplantae"))

# Get a Summary of the read depth after filtering
SeqDepth_noPlants = colSums(otu_table(Americans_Counts_physeq))
sample_data(Americans_Counts_physeq)$SeqDepth_noPlants = SeqDepth_noPlants

# Now let's investigate the Metazoa reads
Metazoa=subset_taxa(Americans_Counts_physeq, (Kingdom=="Metazoa"))
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
Americans_Counts_physeq=subset_taxa(Americans_Counts_physeq, (Kingdom!="Metazoa"))

# Get a Summary of the read depth after filtering
SeqDepth_noAnimals = colSums(otu_table(Americans_Counts_physeq))
sample_data(Americans_Counts_physeq)$SeqDepth_noAnimals = SeqDepth_noAnimals