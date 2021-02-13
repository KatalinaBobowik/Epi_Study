# 29.01.2021 by KSB

# The aim of this script is to fill in empty Phylum values using Family values
# Since the NCBI database is constantly updated, this script is intended to vbe run once, then
# the results from this script will be used in my downstream analyses.
# This script accessed the NCBI datase on January 29th, 2021


# Loading packages and colour setup

library(phyloseq) 
library(tidyverse)
library(taxize)

# set up directories
refdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/EpiStudy/"
filteringDir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/ControlSampleComparison/"
outputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/SEIndonesianSamples/"
batchInfodir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/dataPreprocessing/"


##################################
# Reading in the Indonesian data #
##################################

# Many microbiome studies use the package [phyloseq](https://github.com/joey711/phyloseq) to analyse data due to its comprehensive packages. The data structures in phyloseq (taxa data, otu data, and sample data) are also contained in a single object, which makes it easy to keep everything together.
# Let's read in our Indonesian single-ended, 101bp count data and convert taxa and read abundance information into a phyloseq object. We'll then assign sample information to the phyloseq object.

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
pop <- sapply(strsplit(samplenames, "[-.]"), `[`, 1)

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

# get information on phyloseq objects
AllREadsSE_Indo_Counts_physeq

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
# save replicates file
keepReplicates=rownames(sample_data(AllREadsSE_Indo_Counts_physeq))[-which(rownames(sample_data(AllREadsSE_Indo_Counts_physeq)) %in% removeReplicates)]

# prune these out
AllREadsSE_Indo_Counts_physeq=prune_samples(keepReplicates,AllREadsSE_Indo_Counts_physeq)

# remove taxa with only 0's in the phyloseq object
any(taxa_sums(AllREadsSE_Indo_Counts_physeq) == 0)
AllREadsSE_Indo_Counts_physeq=prune_taxa(taxa_sums(AllREadsSE_Indo_Counts_physeq) > 0, AllREadsSE_Indo_Counts_physeq)

# add sequencing depth information
SeqDepth_Prefilter = colSums(otu_table(AllREadsSE_Indo_Counts_physeq))
sample_data(AllREadsSE_Indo_Counts_physeq)$SeqDepth_Prefilter = SeqDepth_Prefilter

AllREadsSE_Indo_Counts_physeq

# We now have a phyloseq object of 117 samples and 2,763 taxa.

###################
# Data processing #
###################

## Removing singletons from the data ---------

# assign unfiltered phyloseq object to a new object for downstream use
AllREadsSE_Indo_Counts_physeq_withSingletons <- AllREadsSE_Indo_Counts_physeq
# Filter out singletons
otu_table(AllREadsSE_Indo_Counts_physeq)[otu_table(AllREadsSE_Indo_Counts_physeq)<=1]<-0
AllREadsSE_Indo_Counts_physeq <- prune_taxa(taxa_sums(AllREadsSE_Indo_Counts_physeq) > 0, AllREadsSE_Indo_Counts_physeq)
# add sequencing depth information after filtering out singletons
SeqDepth_noSingletons = colSums(otu_table(AllREadsSE_Indo_Counts_physeq))
sample_data(AllREadsSE_Indo_Counts_physeq)$SeqDepth_noSingletons = SeqDepth_noSingletons

## Removing humans and plants ---------------

# From the script 'RemovingViridiplantaeAndMetazoa.R', we saw that human reads and viridiplantae are not of interest and we want to filter these out. When we do this, we'll also keep track of hom many reads we're filtering out after each step.

# Filter out Viridiplantae 
AllREadsSE_Indo_Counts_physeq <- subset_taxa(AllREadsSE_Indo_Counts_physeq, (Kingdom!="Viridiplantae"))
AllREadsSE_Indo_Counts_physeq <- prune_taxa(taxa_sums(AllREadsSE_Indo_Counts_physeq) > 0, AllREadsSE_Indo_Counts_physeq)
# add sequencing depth information after filtering out plants
SeqDepth_noViridiplantae = colSums(otu_table(AllREadsSE_Indo_Counts_physeq))
sample_data(AllREadsSE_Indo_Counts_physeq)$SeqDepth_noViridiplantae = SeqDepth_noViridiplantae

# Filter out Chordata
AllREadsSE_Indo_Counts_physeq <- subset_taxa(AllREadsSE_Indo_Counts_physeq, (Phylum!="Chordata"))
AllREadsSE_Indo_Counts_physeq <- prune_taxa(taxa_sums(AllREadsSE_Indo_Counts_physeq) > 0, AllREadsSE_Indo_Counts_physeq)
# add sequencing depth information after filtering out Metazoa
SeqDepth_noChordata = colSums(otu_table(AllREadsSE_Indo_Counts_physeq))
sample_data(AllREadsSE_Indo_Counts_physeq)$SeqDepth_noChordata = SeqDepth_noChordata

# Filter out Metazoa
AllREadsSE_Indo_Counts_physeq <- subset_taxa(AllREadsSE_Indo_Counts_physeq, (Kingdom!="Metazoa"))
AllREadsSE_Indo_Counts_physeq <- prune_taxa(taxa_sums(AllREadsSE_Indo_Counts_physeq) > 0, AllREadsSE_Indo_Counts_physeq)
# add sequencing depth information after filtering out Metazoa
SeqDepth_noMetazoa = colSums(otu_table(AllREadsSE_Indo_Counts_physeq))
sample_data(AllREadsSE_Indo_Counts_physeq)$SeqDepth_noMetazoa = SeqDepth_noMetazoa

# Filter out taxa which are unassigned at the Kingdom level
AllREadsSE_Indo_Counts_physeq <- subset_taxa(AllREadsSE_Indo_Counts_physeq, (Superkingdom!="unk_sk"))
AllREadsSE_Indo_Counts_physeq <- prune_taxa(taxa_sums(AllREadsSE_Indo_Counts_physeq) > 0, AllREadsSE_Indo_Counts_physeq)
# add sequencing depth information after filtering out Metazoa
SeqDepth_noUnkSk = colSums(otu_table(AllREadsSE_Indo_Counts_physeq))
sample_data(AllREadsSE_Indo_Counts_physeq)$SeqDepth_noUnkSk = SeqDepth_noUnkSk


## Summarising the data ---------------

# Now that we've done all of the filtering, we can plot the final library sizes.

# Finally, fill in missing Phlya names using the package 'taxize'
index=which(tax_table(AllREadsSE_Indo_Counts_physeq)[,"Phylum"]=="unk_p")
specieslist <- c(unname(tax_table(AllREadsSE_Indo_Counts_physeq)[index,"Family"]))
taxize_query=tax_name(query = c(specieslist), get = c("phylum"), db = "ncbi")
taxize_query[is.na(taxize_query)]="unk_p"
tax_table(AllREadsSE_Indo_Counts_physeq)[index,"Phylum"]=taxize_query[,"phylum"]

save(taxize_query, file=paste0(outputdir,"taxize_query.Rda"))

# do this for data with singletons

# Remove Viridiplantae and Metazoa
AllREadsSE_Indo_Counts_physeq_withSingletons <- subset_taxa(AllREadsSE_Indo_Counts_physeq_withSingletons, (Kingdom!="Viridiplantae"))
AllREadsSE_Indo_Counts_physeq_withSingletons <- subset_taxa(AllREadsSE_Indo_Counts_physeq_withSingletons, (Kingdom!="Metazoa"))
AllREadsSE_Indo_Counts_physeq_withSingletons <- subset_taxa(AllREadsSE_Indo_Counts_physeq_withSingletons, (Superkingdom!="unk_sk"))
# remove any empty rows
AllREadsSE_Indo_Counts_physeq_withSingletons <- prune_taxa(taxa_sums(AllREadsSE_Indo_Counts_physeq_withSingletons) > 0, AllREadsSE_Indo_Counts_physeq_withSingletons)

# add in missing taxa using the package Taxize (a few are missing when using the results from CCMetagen)
index=which(tax_table(AllREadsSE_Indo_Counts_physeq_withSingletons)[,"Phylum"]=="unk_p")
specieslist <- c(unname(tax_table(AllREadsSE_Indo_Counts_physeq_withSingletons)[index,"Family"]))
taxize_query_singletons=tax_name(query = c(specieslist), get = c("phylum"), db = "ncbi")
taxize_query_singletons[is.na(taxize_query_singletons)]="unk_p"
tax_table(AllREadsSE_Indo_Counts_physeq_withSingletons)[index,"Phylum"]=taxize_query_singletons[,"phylum"]

# save the data
save(taxize_query_singletons, file=paste0(outputdir,"taxize_query_singletons.Rda"))


