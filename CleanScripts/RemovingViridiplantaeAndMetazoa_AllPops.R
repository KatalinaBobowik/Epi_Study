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
batchInfodir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/dataPreprocessing/"

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

# add in batch information
load(paste0(batchInfodir, "indoRNA.read_counts.TMM.filtered.Rda"))
batch_df = data.frame(rownames(y$samples),y$samples$batch)
colnames(batch_df)=c("Sample","Batch")
batch_df=batch_df[order(batch_df$Sample),]
batch=batch_df[,"Batch"]

# make this into a df and add to the Phloseq object
samples_df=data.frame(SampleName=colnames(otu_table(AllREadsSE_Indo_Counts_physeq)), SamplePop=pop, batch=batch)
samples = sample_data(samples_df)
rownames(samples)=samples$SampleName
sample_data(AllREadsSE_Indo_Counts_physeq) <- samples

# For the Indonesian samples, we can see that we have a phyloseq object consisting of 123 samples, 6 of which are replicates. We'll take the replicates with the highest library depth and then remove the rest.

# get replicates
trimmed_samplenames = gsub("Batch1",'',samplenames) %>% gsub("Batch2",'', .) %>% gsub("Batch3",'', .)
trimmed_samplenames = sub("([A-Z]{3})([0-9]{3})", "\\1-\\2", trimmed_samplenames)
replicate_index = which(duplicated(trimmed_samplenames) | duplicated(trimmed_samplenames, fromLast = TRUE))
replicates = samplenames[replicate_index]

# add sequencing depth information to the Physeq object in order to filter replicates by seqDepth
SeqDepth = colSums(otu_table(AllREadsSE_Indo_Counts_physeq))
sample_data(AllREadsSE_Indo_Counts_physeq)$SeqDepth = SeqDepth

# find out which replicates have the highest sequencing depth
sample_data(AllREadsSE_Indo_Counts_physeq)[replicates,]
replicateDF=as.data.frame(sample_data(AllREadsSE_Indo_Counts_physeq)[replicates,])
replicateDF$SampleName = sub("([A-Z]{3})([0-9]{3})", "\\1-\\2", replicateDF$SampleName)
replicateDF$SampleName = gsub("Batch1",'',replicateDF$SampleName) %>% gsub("Batch2",'', .) %>% gsub("Batch3",'', .)
replicateDF=replicateDF[with(replicateDF, order(-SeqDepth)), ]
removeReplicates=rownames(replicateDF[which(duplicated(replicateDF$SampleName)),])
keepReplicates=rownames(sample_data(AllREadsSE_Indo_Counts_physeq))[-which(rownames(sample_data(AllREadsSE_Indo_Counts_physeq)) %in% removeReplicates)]

# prune these out
AllREadsSE_Indo_Counts_physeq=prune_samples(keepReplicates,AllREadsSE_Indo_Counts_physeq)

# remove taxa with only 0's in the phyloseq object
any(taxa_sums(AllREadsSE_Indo_Counts_physeq) == 0)
AllREadsSE_Indo_Counts_physeq=prune_taxa(taxa_sums(AllREadsSE_Indo_Counts_physeq) > 0, AllREadsSE_Indo_Counts_physeq)

# filter reads 1 and below
otu_table(AllREadsSE_Indo_Counts_physeq)[otu_table(AllREadsSE_Indo_Counts_physeq)<=1]<-0

# Get a Summary of the read depth before any filtering
SeqDepth_preFilter = colSums(otu_table(AllREadsSE_Indo_Counts_physeq))
sample_data(AllREadsSE_Indo_Counts_physeq)$SeqDepth_preFilter = SeqDepth_preFilter

# we now have a phyloseq object with 117 samples
AllREadsSE_Indo_Counts_physeq

# investigate where reads are mapping to
# p=plot_bar(AllREadsSE_Indo_Counts_physeq, fill = "Phylum")
# fig <- p + geom_bar(stat="identity", position="stack")
# fig

# Some human reads are getting through. We also see a low level of plants
# from Streptophya sequences
Viridiplantae=subset_taxa(AllREadsSE_Indo_Counts_physeq, (Kingdom=="Viridiplantae"))
p=plot_bar(Viridiplantae, fill = "Phylum")
pdf(paste0(outputdir,"Viridiplantae_Indonesians.pdf"))
fig <- p + geom_bar(stat="identity", position="stack") + theme_bw(base_size = 15) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
	ggtitle("Viridiplantae", subtitle = "Indonesia") + ylab("Counts") + xlab("Samples")
fig
dev.off()

# Remove Viridiplantae and get number of sequences after removal
AllREadsSE_Indo_Counts_physeq=subset_taxa(AllREadsSE_Indo_Counts_physeq, (Kingdom!="Viridiplantae"))

# Get a Summary of the read depth after filtering
SeqDepth_noPlants = colSums(otu_table(AllREadsSE_Indo_Counts_physeq))
sample_data(AllREadsSE_Indo_Counts_physeq)$SeqDepth_noPlants = SeqDepth_noPlants

# Now let's investigate the Metazoa reads
Metazoa=subset_taxa(AllREadsSE_Indo_Counts_physeq, (Kingdom=="Metazoa"))
p=plot_bar(Metazoa, fill = "Phylum")
pdf(paste0(outputdir,"Metazoa_Indonesians.pdf"))
fig <- p + geom_bar(stat="identity", position="stack") + theme_bw(base_size = 15) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
	ggtitle("Metazoa", subtitle = "Indonesia") + ylab("Counts") + xlab("Samples")
fig
dev.off()

# we see that most of the metazoa reads aru humans that got through.
# Let's remove them and plot again
pdf(paste0(outputdir,"Metazoa_noChordata_Indonesians.pdf"))
Metazoa=subset_taxa(Metazoa, (Phylum!="Chordata"))
# filter out any 0's that remain in the dataframe
Metazoa <- prune_taxa(taxa_sums(Metazoa) > 0, Metazoa)
p=plot_bar(Metazoa, fill = "Phylum")
fig <- p + geom_bar(stat="identity", position="stack") + theme_bw(base_size = 15) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
	ggtitle("Metazoa, no Chordata", subtitle = "Indonesia") + ylab("Counts") + xlab("Samples")
fig
dev.off()

# only platyhelminths and nematodes
pdf(paste0(outputdir,"Metazoa_Platyhelminths_Nematodes_Indonesians.pdf"))
Worms=subset_taxa(Metazoa, Phylum=="Platyhelminthes" | Phylum=="Nematoda" | Phylum=="Annelida")
# filter out any 0's that remain in the dataframe
Worms <- prune_taxa(taxa_sums(Worms) > 0, Worms)
p=plot_bar(Worms, fill = "Phylum")
fig <- p + geom_bar(stat="identity", position="stack") + theme_bw(base_size = 15) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
	ggtitle("Helminths", subtitle = "Indonesia") + ylab("Counts") + xlab("Samples")
fig
dev.off()


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
# p=plot_bar(Mali_Counts_physeq, fill = "Phylum")
# fig <- p + geom_bar(stat="identity", position="stack")
# fig

# Some human reads are getting through. We also see a low level of plants
# from Streptophya sequences
Viridiplantae=subset_taxa(Mali_Counts_physeq, (Kingdom=="Viridiplantae"))
p=plot_bar(Viridiplantae, fill = "Phylum")
pdf(paste0(outputdir,"Viridiplantae_Mali.pdf"))
fig <- p + geom_bar(stat="identity", position="stack") + theme_bw(base_size = 15) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
	ggtitle("Viridiplantae", subtitle = "Mali") + ylab("Counts") + xlab("Samples")
fig
dev.off()

# Remove Viridiplantae and get number of sequences after removal
Mali_Counts_physeq=subset_taxa(Mali_Counts_physeq, (Kingdom!="Viridiplantae"))

# Get a Summary of the read depth after filtering
SeqDepth_noPlants = colSums(otu_table(Mali_Counts_physeq))
sample_data(Mali_Counts_physeq)$SeqDepth_noPlants = SeqDepth_noPlants

# Now let's investigate the Metazoa reads
Metazoa=subset_taxa(Mali_Counts_physeq, (Kingdom=="Metazoa"))
p=plot_bar(Metazoa, fill = "Phylum")
pdf(paste0(outputdir,"Metazoa_Mali.pdf"))
fig <- p + geom_bar(stat="identity", position="stack") + theme_bw(base_size = 15) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
	ggtitle("Metazoa", subtitle = "Mali") + ylab("Counts") + xlab("Samples")
fig
dev.off()

# we see that most of the metazoa reads aru humans that got through.
# Let's remove them and plot again
Metazoa=subset_taxa(Metazoa, (Phylum!="Chordata"))
# filter out any 0's that remain in the dataframe
Metazoa <- prune_taxa(taxa_sums(Metazoa) > 0, Metazoa)
p=plot_bar(Metazoa, fill = "Phylum")
pdf(paste0(outputdir,"Metazoa_noChordata_Mali.pdf"))
fig <- p + geom_bar(stat="identity", position="stack") + theme_bw(base_size = 15) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
	ggtitle("Metazoa, no Chordata", subtitle = "Mali") + ylab("Counts") + xlab("Samples")
fig
dev.off()

# only platyhelminths and nematodes
pdf(paste0(outputdir,"Metazoa_Platyhelminths_Nematodes_Mali.pdf"))
Worms=subset_taxa(Metazoa, Phylum=="Platyhelminthes" | Phylum=="Nematoda" | Phylum=="Annelida")
# filter out any 0's that remain in the dataframe
Worms <- prune_taxa(taxa_sums(Worms) > 0, Worms)
p=plot_bar(Worms, fill = "Phylum")
fig <- p + geom_bar(stat="identity", position="stack") + theme_bw(base_size = 15) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
	ggtitle("Helminths", subtitle = "Mali") + ylab("Counts") + xlab("Samples")
fig
dev.off()

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
# UK Samples #
################

European_Counts <- read.csv(paste0(refdir,"TBControls_NoFiltering_Counts.csv"),check.names=FALSE)

# Separate species' abundances and taxonomy columns
taxa_raw <- as.matrix(European_Counts[,c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species")])
abund_raw <- as.matrix(European_Counts[,-which(colnames(European_Counts) %in% c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species"))])
# convert to Phyloseq object
tax = tax_table(taxa_raw)
taxa = otu_table(abund_raw, taxa_are_rows = TRUE)
European_Counts_physeq = phyloseq(taxa, tax)

# add in sample information, i.e., the sample names and population they're from
samplenames <- colnames(otu_table(European_Counts_physeq))
pop <- rep("UK",ncol(otu_table(European_Counts_physeq)))

# make this into a df and add to the Phloseq object
samples_df=data.frame(SampleName=colnames(otu_table(European_Counts_physeq)), SamplePop=pop)
samples = sample_data(samples_df)
rownames(samples)=samples$SampleName
sample_data(European_Counts_physeq) <- samples

# filter reads 1 and below
otu_table(European_Counts_physeq)[otu_table(European_Counts_physeq)<=1]<-0

# Get a Summary of the read depth before any filtering
SeqDepth_preFilter = colSums(otu_table(European_Counts_physeq))
sample_data(European_Counts_physeq)$SeqDepth_preFilter = SeqDepth_preFilter

# investigate where reads are mapping to
# p=plot_bar(European_Counts_physeq, fill = "Phylum")
# fig <- p + geom_bar(stat="identity", position="stack")
# fig

# Some human reads are getting through. We also see a low level of plants
# from Streptophya sequences
Viridiplantae=subset_taxa(European_Counts_physeq, (Kingdom=="Viridiplantae"))
p=plot_bar(Viridiplantae, fill = "Phylum")
pdf(paste0(outputdir,"Viridiplantae_UK.pdf"))
fig <- p + geom_bar(stat="identity", position="stack") + theme_bw(base_size = 15) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
	ggtitle("Viridiplantae", subtitle = "UK") + ylab("Counts") + xlab("Samples")
fig
dev.off()

# Remove Viridiplantae and get number of sequences after removal
European_Counts_physeq=subset_taxa(European_Counts_physeq, (Kingdom!="Viridiplantae"))

# Get a Summary of the read depth after filtering
SeqDepth_noPlants = colSums(otu_table(European_Counts_physeq))
sample_data(European_Counts_physeq)$SeqDepth_noPlants = SeqDepth_noPlants

# Now let's investigate the Metazoa reads
Metazoa=subset_taxa(European_Counts_physeq, (Kingdom=="Metazoa"))
p=plot_bar(Metazoa, fill = "Phylum")
pdf(paste0(outputdir,"Metazoa_UK.pdf"))
fig <- p + geom_bar(stat="identity", position="stack") + theme_bw(base_size = 15) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
	ggtitle("Metazoa", subtitle = "UK") + ylab("Counts") + xlab("Samples")
fig
dev.off()

# we see that most of the metazoa reads aru humans that got through.
# Let's remove them and plot again
Metazoa=subset_taxa(Metazoa, (Phylum!="Chordata"))
# filter out any 0's that remain in the dataframe
Metazoa <- prune_taxa(taxa_sums(Metazoa) > 0, Metazoa)
p=plot_bar(Metazoa, fill = "Phylum")
pdf(paste0(outputdir,"Metazoa_noChordata_UK.pdf"))
fig <- p + geom_bar(stat="identity", position="stack") + theme_bw(base_size = 15) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
	ggtitle("Metazoa, no Chordata", subtitle = "UK") + ylab("Counts") + xlab("Samples")
fig
dev.off()

# only platyhelminths and nematodes
pdf(paste0(outputdir,"Metazoa_Platyhelminths_Nematodes_UK.pdf"))
Worms=subset_taxa(Metazoa, Phylum=="Platyhelminthes" | Phylum=="Nematoda" | Phylum=="Annelida")
# filter out any 0's that remain in the dataframe
Worms <- prune_taxa(taxa_sums(Worms) > 0, Worms)
p=plot_bar(Worms, fill = "Phylum")
fig <- p + geom_bar(stat="identity", position="stack") + theme_bw(base_size = 15) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
	ggtitle("Helminths", subtitle = "UK") + ylab("Counts") + xlab("Samples")
fig
dev.off()

# We can see that arthropods (probably caused by mismapping), nematodes,
# and platyhelminths are found in the data. 

# Blasted sequences for Spirometra erinaceieuropaei (found in MTW-MDB-015), and Caenorhabditis remanei,
# found in MTW-TLL-006. These sequences both have reads that map equally well to human
# sequences. Thus, I will remove these in particular since helminths are known pathogens
# and I don't want the results to look misleading.

# Remove Metazoa and get number of sequences after removal
European_Counts_physeq=subset_taxa(European_Counts_physeq, (Kingdom!="Metazoa"))

# Get a Summary of the read depth after filtering
SeqDepth_noAnimals = colSums(otu_table(European_Counts_physeq))
sample_data(European_Counts_physeq)$SeqDepth_noAnimals = SeqDepth_noAnimals