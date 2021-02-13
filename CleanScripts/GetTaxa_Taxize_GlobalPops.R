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
outputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/ControlSampleComparison_TB/"

# Reading in the Indonesian data

# First, we'll read in our Indonesian single-ended count data and convert taxa and read abundance information into a phyloseq object. We'll then assign sample information to the phyloseq object so that we can compare it to the control datasets.
AllREadsSE_Indo_Counts <- read.csv(paste0(refdir,"Counts_Indo_75BP_SE.csv"),check.names=FALSE)

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

# add sequencing depth information
SeqDepth_Prefilter = colSums(otu_table(AllREadsSE_Indo_Counts_physeq))
sample_data(AllREadsSE_Indo_Counts_physeq)$SeqDepth_Prefilter = SeqDepth_Prefilter

AllREadsSE_Indo_Counts_physeq

# We now have a phyloseq object of 117 samples and 2,400 taxa.

####################
# control datasets #
####################

# The first control dataset we'll read in is from a study conducted by [Tran et al](https://www.nature.com/articles/srep31291#Tab1), conducted in 2016. This study looked at transcriptomic differences between individuals pre and post (natural) infection with P. falciparum. We'll only be using the healthy samples pre-infection, of which there are 54.
# One of the reasons why this dataset is interesting is becase it comes from individuals living in areas that have a high pathogen load, similar to our study. Therefore, we can compare how our data looks like to another population that lives in an area with a high rate of endemic malaria and high pathogen load.
# Let's read in the data, make it into a PhyloSeq object, and add in our sample information. 

## Malian Samples ---------------

Mali_Counts <- read.csv(paste0(refdir,"MaliControls_Counts_NoFilter_75BP.csv"),check.names=FALSE)

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

# add sequencing depth information 
SeqDepth_Prefilter = colSums(otu_table(Mali_Counts_physeq))
sample_data(Mali_Counts_physeq)$SeqDepth_Prefilter = SeqDepth_Prefilter

# summarise data
Mali_Counts_physeq

# We now have a phyloseq object of 54 samples by 3355 taxa.

## European samples -------------

# The next set of control samples is from a study conducted by [Singhania et al](https://www.nature.com/articles/s41467-018-04579-w) in 2018. This study was conducted on whole blood of tubercolosis (TB) and latent-(TB) patients to distinguish between active and asymptomatic TB infection, however the samples used for this analysis are all healthy controls.
# As with the other studies, we'll first convert the data to a phyloseq object.

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

# add sequencing depth information 
SeqDepth_Prefilter = colSums(otu_table(European_Counts_physeq))
sample_data(European_Counts_physeq)$SeqDepth_Prefilter = SeqDepth_Prefilter

# get phyloseq summary information
European_Counts_physeq

# We now have a phyloseq object of 12 samples with 1100 taxa.

####################
# Merging the data #
####################

# The next step is to merge the Indonesian and control datasets together. Phyloseq makes this relatively easy with its function merge_phyloseq(). 
# HOWEVER, there is a caveat to this: although the [merge_phyloseq function](https://www.rdocumentation.org/packages/phyloseq/versions/1.16.2/topics/merge_phyloseq) says that it merges by 'first separating higher-order objects into a list of their component objects; then, merging any component objects of the same class into one object', this is not the case (for my data, anyhow)!! I'm not the only one who [had this problem](https://github.com/joey711/phyloseq/issues/574), so in the end, the solution is to give each taxa name a unique ID, which in this case, is the species name (the most specific name). We'll then make sure that our phyloseq objects from above (before merging) are the same as after we merge.

# assign unique taxa names to all phyloseq objects
taxa_names(AllREadsSE_Indo_Counts_physeq) <- paste(tax_table(AllREadsSE_Indo_Counts_physeq)[,"Superkingdom"], tax_table(AllREadsSE_Indo_Counts_physeq)[,"Kingdom"], tax_table(AllREadsSE_Indo_Counts_physeq)[,"Phylum"], tax_table(AllREadsSE_Indo_Counts_physeq)[,"Class"], tax_table(AllREadsSE_Indo_Counts_physeq)[,"Order"], tax_table(AllREadsSE_Indo_Counts_physeq)[,"Family"], tax_table(AllREadsSE_Indo_Counts_physeq)[,"Genus"], tax_table(AllREadsSE_Indo_Counts_physeq)[,"Species"], sep="_")
taxa_names(Mali_Counts_physeq) <- make.unique(paste(tax_table(Mali_Counts_physeq)[,"Superkingdom"], tax_table(Mali_Counts_physeq)[,"Kingdom"], tax_table(Mali_Counts_physeq)[,"Phylum"], tax_table(Mali_Counts_physeq)[,"Class"], tax_table(Mali_Counts_physeq)[,"Order"], tax_table(Mali_Counts_physeq)[,"Family"], tax_table(Mali_Counts_physeq)[,"Genus"], tax_table(Mali_Counts_physeq)[,"Species"], sep="_"))
taxa_names(European_Counts_physeq) <- make.unique(paste(tax_table(European_Counts_physeq)[,"Superkingdom"], tax_table(European_Counts_physeq)[,"Kingdom"], tax_table(European_Counts_physeq)[,"Phylum"], tax_table(European_Counts_physeq)[,"Class"], tax_table(European_Counts_physeq)[,"Order"], tax_table(European_Counts_physeq)[,"Family"], tax_table(European_Counts_physeq)[,"Genus"], tax_table(European_Counts_physeq)[,"Species"], sep="_"))
merged_phylo_counts=merge_phyloseq(AllREadsSE_Indo_Counts_physeq, Mali_Counts_physeq, European_Counts_physeq)

# subset populations
Indonesian_subset <- phyloseq::subset_samples(merged_phylo_counts, SamplePop == "Indonesia")
Mali_subset <- phyloseq::subset_samples(merged_phylo_counts, SamplePop == "Mali")
UK_subset <- phyloseq::subset_samples(merged_phylo_counts, SamplePop == "UK")

# define function to check if both objects are identical
is.identical <- function(pop_to_subset, original_phyloseq_obj){
	pruned_pop_subset <- prune_taxa(taxa_sums(pop_to_subset) > 0, pop_to_subset)
	merged_df=merge(as.data.frame(otu_table(pruned_pop_subset)), as.data.frame(otu_table(original_phyloseq_obj)),by=0)
	nsamples=ncol(otu_table(pruned_pop_subset))
	merged_df[,"Row.names"] <- NULL
	colnames(merged_df)=rep("Rep",ncol(merged_df))
	is.identical <- identical(merged_df[,1:nsamples],merged_df[,(nsamples+1):ncol(merged_df)])
	return(is.identical)
}

is.identical(pop_to_subset=Indonesian_subset, original_phyloseq_obj=AllREadsSE_Indo_Counts_physeq)
# TRUE
is.identical(pop_to_subset=Mali_subset, original_phyloseq_obj=Mali_Counts_physeq)
# TRUE
is.identical(pop_to_subset=UK_subset, original_phyloseq_obj=European_Counts_physeq)
# TRUE

# They're all the same, so we can move forward.

###################
# Data processing #
###################

# Removing singletons from the data -----------

# The easiest way to get rid of some error in your data is to throw out any count information below some threshold. Oddly, in microbiomics, there's no set thresholding for this. In the end, it's really a compromise between accuracy and keeping rare taxa. What *is* decided is that filtering out at least singletons is standard, since these are regarded as sources of error or contamination. Some resources I found that were helpful on this can be found [here](http://drive5.com/usearch/manual/singletons.html) and [here](https://forum.qiime2.org/t/do-you-guys-still-remove-singletons-or-doubletons-these-days/7138/2).  
# In most of the datasets, we can see that there's a high number of low counts within the data.

# Because our starting read depth is small, we will stick with removing singletons. We will also add the sequencing depth information to the phyloseq object to keep track of oir library size after filtering.

# Filter out singletons
otu_table(merged_phylo_counts)[otu_table(merged_phylo_counts)<=1]<-0
merged_phylo_counts <- prune_taxa(taxa_sums(merged_phylo_counts) > 0, merged_phylo_counts)
# add sequencing depth information after filtering out singletons
SeqDepth_noSingletons = colSums(otu_table(merged_phylo_counts))
sample_data(merged_phylo_counts)$SeqDepth_noSingletons = SeqDepth_noSingletons

## Removing humans and plants ------------

# From the script 'RemovingViridiplantaeAndMetazoa.R', we saw that human reads and viridiplantae are not of interest and we want to filter these out. When we do this, we'll also keep track of hom many reads we're filtering out after each step.

# Filter out Viridiplantae 
merged_phylo_counts <- subset_taxa(merged_phylo_counts, (Kingdom!="Viridiplantae"))
merged_phylo_counts <- prune_taxa(taxa_sums(merged_phylo_counts) > 0, merged_phylo_counts)
# add sequencing depth information after filtering out plants
SeqDepth_noViridiplantae = colSums(otu_table(merged_phylo_counts))
sample_data(merged_phylo_counts)$SeqDepth_noViridiplantae = SeqDepth_noViridiplantae

# Filter out Chordata
merged_phylo_counts <- subset_taxa(merged_phylo_counts, (Phylum!="Chordata"))
merged_phylo_counts <- prune_taxa(taxa_sums(merged_phylo_counts) > 0, merged_phylo_counts)
# add sequencing depth information after filtering out Metazoa
SeqDepth_noChordata = colSums(otu_table(merged_phylo_counts))
sample_data(merged_phylo_counts)$SeqDepth_noChordata = SeqDepth_noChordata

# Filter out Metazoa
merged_phylo_counts <- subset_taxa(merged_phylo_counts, (Kingdom!="Metazoa"))
merged_phylo_counts <- prune_taxa(taxa_sums(merged_phylo_counts) > 0, merged_phylo_counts)
# add sequencing depth information after filtering out Metazoa
SeqDepth_noMetazoa = colSums(otu_table(merged_phylo_counts))
sample_data(merged_phylo_counts)$SeqDepth_noMetazoa = SeqDepth_noMetazoa

# Filter out taxa which are unassigned at the Kingdom level
merged_phylo_counts <- subset_taxa(merged_phylo_counts, (Superkingdom!="unk_sk"))
merged_phylo_counts <- prune_taxa(taxa_sums(merged_phylo_counts) > 0, merged_phylo_counts)
# add sequencing depth information after filtering out Metazoa
SeqDepth_noUnkSk = colSums(otu_table(merged_phylo_counts))
sample_data(merged_phylo_counts)$SeqDepth_noUnkSk = SeqDepth_noUnkSk

## Summarising the data -------------

# We can see that in the Indonesian and Malian dataset, most of the reads are removed when removing Chordates, however for the UK dataset, most reads are removed when removing Metazoa. For the UK dataset, this is due to a high number of reads mapping to molluscs, as I discuss in my script 'RemovingViridiplantaeAndMetazoa.R'.

# Finally, fill in taxa that CCMetagen couldn't identify with taxize
index=which(tax_table(merged_phylo_counts)[,"Phylum"]=="unk_p")
specieslist <- c(unname(tax_table(merged_phylo_counts)[index,"Family"]))
taxize_query=tax_name(query = c(specieslist), get = c("phylum"), db = "ncbi")
taxize_query[is.na(taxize_query)]="unk_p"
tax_table(merged_phylo_counts)[index,"Phylum"]=taxize_query[,"phylum"]

# save the data
save(taxize_query, file=paste0(outputdir,"taxize_query.Rda"))

# do this for data with singletons

# Get phyloseq object without singletons by merging original phyloseq objects
merged_phylo_counts_withSingletons <- merge_phyloseq(AllREadsSE_Indo_Counts_physeq, Mali_Counts_physeq, European_Counts_physeq)
# Remove Viridiplantae and Metazoa
merged_phylo_counts_withSingletons <- subset_taxa(merged_phylo_counts_withSingletons, (Kingdom!="Viridiplantae"))
merged_phylo_counts_withSingletons <- subset_taxa(merged_phylo_counts_withSingletons, (Kingdom!="Metazoa"))
merged_phylo_counts_withSingletons <- subset_taxa(merged_phylo_counts_withSingletons, (Superkingdom!="unk_sk"))
# remove any empty rows
merged_phylo_counts_withSingletons <- prune_taxa(taxa_sums(merged_phylo_counts_withSingletons) > 0, merged_phylo_counts_withSingletons)

# fill in taxa that CCMetagen couldn't identify with taxize
index=which(tax_table(merged_phylo_counts_withSingletons)[,"Phylum"]=="unk_p")
specieslist <- c(unname(tax_table(merged_phylo_counts_withSingletons)[index,"Family"]))
taxize_query_singletons=tax_name(query = c(specieslist), get = c("phylum"), db = "ncbi")
taxize_query_singletons[is.na(taxize_query_singletons)]="unk_p"
tax_table(merged_phylo_counts_withSingletons)[index,"Phylum"]=taxize_query_singletons[,"phylum"]

# save the data
save(taxize_query_singletons, file=paste0(outputdir,"taxize_query_singletons.Rda"))





