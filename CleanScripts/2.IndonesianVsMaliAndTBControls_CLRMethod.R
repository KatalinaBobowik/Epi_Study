# Comparing Indonesian Samples Vs. Controls
# 2020-09-13 by KSB


# This part of the study will analyse differences between our Indonesian samples and compare them with healthy controls. From my previous analysis of healthy Indonesians, I found that Plasmodium, Flavivirus, and Bacteria make up the majority of the taxa found within unmapped reads. However, I want to see how different this is to populations in more sterile environments. 
# The aim of this analysis is therefore to test whether the pathogen signature identified in our Indonesian samples is unique. We will test this by looking at the sample clustering, relative abundance of taxa, differential abundance testing, and diversity estimates. 


# Loading packages and colour setup

require(ggplot2)
require(RColorBrewer)
library(dplyr)
library(plyr)
library(reshape2)
library(ggpubr)
library(metacoder)
library(tidyverse)             
library(phyloseq)                   
library(DESeq2)                       
library(microbiome)               
library(vegan)                         
library(picante)                     
library(ALDEx2)                      
library(metagenomeSeq)          
library(HMP)                             
library(dendextend)               
library(selbal)                       
library(rms)
library(breakaway)        
library(microbiomeutilities)
library(mixOmics)
library(SRS)
library(ggrepel)
library(DivNet)
library(taxize)

# set up directories
refdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/EpiStudy/"
filteringDir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/ControlSampleComparison/"
outputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/ControlSampleComparison_TB/"

# set ggplot colour theme to white
theme_set(theme_bw())

# Set up colour schemes
IndonesiaCol="#332288"
MaliCol="#999933"
UKCol="#CC6677"

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
save(Mali_Counts_physeq, file=paste0(outputdir,"Mali_Counts_physeq.Rda"))

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
save(European_Counts_physeq, file=paste0(outputdir,"European_Counts_physeq.Rda"))

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

# histogram of data
pdf(paste0(outputdir,"seqDepthHistpgram.pdf"))
ggplot(meta(merged_phylo_counts)) + geom_histogram(aes(x = SeqDepth_Prefilter), alpha= 0.6, bins=100) + facet_wrap(~SamplePop)
dev.off()

# See how rarefaction looks like when we remove singletons and when we remove counts 5 and below.

# # Separate species' abundances and taxonomy columns
# rarecurve_counts <- otu_table(merged_phylo_counts)
# col <- c(rep(IndonesiaCol,sum(sample_data(merged_phylo_counts)[,"SamplePop"]=="Indonesia")),rep(MaliCol,sum(sample_data(merged_phylo_counts)[,"SamplePop"]=="Mali")),rep(UKCol,sum(sample_data(merged_phylo_counts)[,"SamplePop"]=="UK")))
# # Try with different filtering thresholds:
# pdf(paste0(outputdir,"rarefaction.pdf"))
# for (i in c(1,5)){
#  	rarecurve_counts[rarecurve_counts<=i]<-0
# 	rarecurve(t(otu_table(rarecurve_counts, taxa_are_rows = TRUE)), step=200, col=col,label=F, xlab="Counts",ylab="Number Species",main=paste("Removing Reads",i,"And Below",sep=" "),xlim=c(0,50000))
# }
# dev.off()

# # Rarefaction curves for the Indonesian dataset only (since it is lower sequencing depth and hard to see)
# Indonesian_subset <- prune_taxa(taxa_sums(Indonesian_subset) > 0, Indonesian_subset)
# rarecurve_counts <- otu_table(Indonesian_subset)
# col <- IndonesiaCol
# pdf(paste0(outputdir,"rarefaction_IndoOnly.pdf"))
# for (i in c(1,5)){
#  	rarecurve_counts[rarecurve_counts<=i]<-0
# 	rarecurve(t(otu_table(rarecurve_counts, taxa_are_rows = TRUE)), step=200, col=col,label=F, xlab="Counts",ylab="Number Species",main=paste("Removing Reads",i,"And Below \nIndonesian Samples",sep=" "),xlim=c(0,10000))
# }
# dev.off()

# Because our starting read depth is small, we will stick with removing singletons. We will also add the sequencing depth information to the phyloseq object to keep track of oir library size after filtering.

# Filter out singletons
otu_table(merged_phylo_counts)[otu_table(merged_phylo_counts)<=1]<-0
merged_phylo_counts <- prune_taxa(taxa_sums(merged_phylo_counts) > 0, merged_phylo_counts)
# add sequencing depth information after filtering out singletons
SeqDepth_noSingletons = colSums(otu_table(merged_phylo_counts))
sample_data(merged_phylo_counts)$SeqDepth_noSingletons = SeqDepth_noSingletons

save(merged_phylo_counts, file=paste0(outputdir,"merged_phylo_counts_noSingletons.Rda"))

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

# save the data
save(merged_phylo_counts, file=paste0(outputdir,"merged_phylo_counts.Rda"))

## Summarising the data -------------

# Now that we've done all of the filtering, we can plot the final library sizes.

# barplot of library sizes
pdf(paste0(outputdir,"librarySizes.pdf"))
ggplot(meta(merged_phylo_counts), aes(SampleName, SeqDepth_noUnkSk)) + geom_bar(stat = "identity", aes(fill = SamplePop)) +
scale_fill_manual(values = c(IndonesiaCol,MaliCol,UKCol)) + rotate_x_text()
dev.off()

# We can see that the library sizes are highly uneven, with the Indonesian data having the lowest sampling depth (with the exception of a few samples) and the UK dataset having the highest library depth.
# The final step us is to summarise the data and see how many reads we lost at each filtering step.

FilteringSummary = sample_data(merged_phylo_counts)[,c("SamplePop","SeqDepth_Prefilter","SeqDepth_noSingletons","SeqDepth_noViridiplantae","SeqDepth_noChordata","SeqDepth_noMetazoa", "SeqDepth_noUnkSk")]
# save filtering summary
write.table(FilteringSummary, file=paste0(outputdir,"FilteringSummary.txt"))

# melt df and plot
melted_FilteringSummary = melt(FilteringSummary)
pdf(paste0(outputdir,"FilteringSummary.pdf"))
ggplot(melted_FilteringSummary, aes(x=variable, y=log(value), fill=SamplePop)) +
  geom_violin(alpha=0.8) + theme_bw() + ylab("Spearman pairwise correlation") +
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90)) + scale_fill_manual(values = c(IndonesiaCol,MaliCol,UKCol)) +
  geom_boxplot(color="black",width=0.2, alpha = 0.7) + facet_wrap(~ SamplePop, scales = "free")
dev.off()

# We can see that in the Indonesian and Malian dataset, most of the reads are removed when removing Chordates, however for the UK dataset, most reads are removed when removing Metazoa. For the UK dataset, this is due to a high number of reads mapping to molluscs, as I discuss in my script 'RemovingViridiplantaeAndMetazoa.R'.

# Finally, fill in taxa that CCMetagen couldn't identify with taxize
index=which(tax_table(merged_phylo_counts)[,"Phylum"]=="unk_p")
specieslist <- c(unname(tax_table(merged_phylo_counts)[index,"Family"]))
a=tax_name(query = c(specieslist), get = c("phylum"), db = "ncbi")
a[is.na(a)]="unk_p"
tax_table(merged_phylo_counts)[index,"Phylum"]=a[,"phylum"]

################
# Data summary #
################

# First, UK population
UK=subset_samples(merged_phylo_counts, SamplePop == "UK")
any(taxa_sums(UK) == 0)
# remove any 0s from the data
UK <- prune_taxa(taxa_sums(UK) > 0, UK)
sort(rowSums(otu_table(UK))/sum(rowSums(otu_table(UK)))*100, decreasing = T)

# now Malian population at the Phylum level
data_summary_UK_Family <- UK %>%
  tax_glom("Family")

taxa_names(data_summary_UK_Family)=make.unique(paste(tax_table(data_summary_UK_Family)[,c("Superkingdom")], tax_table(data_summary_UK_Family)[,c("Phylum")], tax_table(data_summary_UK_Family)[,c("Class")], tax_table(data_summary_UK_Family)[,c("Order")], tax_table(data_summary_UK_Family)[,c("Family")], sep="_"))
rev(sort(rowSums(otu_table(data_summary_UK_Family))/sum(rowSums(otu_table(data_summary_UK_Family)))))*100

# now Malian population at the Phylum level
data_summary_UK_Phylum <- UK %>%
  tax_glom("Phylum")

taxa_names(data_summary_UK_Phylum)=make.unique(paste(tax_table(data_summary_UK_Phylum)[,c("Superkingdom")], tax_table(data_summary_UK_Phylum)[,c("Phylum")], sep="_"))
rev(sort(rowSums(otu_table(data_summary_UK_Phylum))/sum(rowSums(otu_table(data_summary_UK_Phylum)))))*100


Mali=subset_samples(merged_phylo_counts, SamplePop == "Mali")
any(taxa_sums(UK) == 0)
# remove any 0s from the data
Mali <- prune_taxa(taxa_sums(Mali) > 0, Mali)

# now Malian population at the Phylum level
data_summary_Mali_Phylum <- Mali %>%
  tax_glom("Phylum")

taxa_names(data_summary_Mali_Phylum )=make.unique(paste(tax_table(data_summary_Mali_Phylum )[,c("Superkingdom")], tax_table(data_summary_Mali_Phylum )[,c("Phylum")], sep="_"))
rev(sort(rowSums(otu_table(data_summary_Mali_Phylum ))/sum(rowSums(otu_table(data_summary_Mali_Phylum )))))*100

# # now Indonesian population at the Phylum level
# data_summary_allPops <- merged_phylo_counts %>%
#   tax_glom("Phylum")

# taxa_names(data_summary_allPops)=make.unique(paste(tax_table(data_summary_allPops)[,c("Superkingdom")], tax_table(data_summary_allPops)[,c("Phylum")], sep="_"))

# Indo=subset_samples(data_summary_allPops, SamplePop == "Indonesia")
# any(taxa_sums(Indo) == 0)
# # remove any 0s from the data
# Indo <- prune_taxa(taxa_sums(Indo) > 0, Indo)
# sort(rowSums(otu_table(Indo))/sum(rowSums(otu_table(Indo)))*100, decreasing = T)


###################### 
# Data normalisation #
######################

# The library sizes between samples and groups is highly variable, and therefore comparing the data to each other will result in biased results. 

# There are two ways of handling this: 
# 1. Performing a transformation of the data or
# 2. rarefying the data.

## Centered log-ration transformation -----------------


### Sample grouping ------------

# Now that we've transformed our data, we can make a PCA plot to see how each sample clusters. The current obect we have is a CLR-class object. You can plot this type of data object easily with mixOmics, however I prefer the visualisation that phyloseq offers (you can't alter the PCA plots that much in mixOmics). So, we'll turn the clr object back into a phyloseq object and make an ordination plot of the data.
# !Note: When Euclidean distances are used in PCoA plots, it is [equivalent to a PCA plot](http://ordination.okstate.edu/overview.htm). 

pop_comparison <- merged_phylo_counts %>%
  tax_glom("Phylum")

pop_comparison

# otu_table()   OTU Table:         [ 50 taxa and 183 samples ]
# sample_data() Sample Data:       [ 183 samples by 9 sample variables ]
# tax_table()   Taxonomy Table:    [ 50 taxa by 8 taxonomic ranks ]

# change taxa names to reflect Phylum level
taxa_names(pop_comparison) <- paste(tax_table(pop_comparison)[,"Superkingdom"], tax_table(pop_comparison)[,"Kingdom"], tax_table(pop_comparison)[,"Phylum"], sep="_")

zComposition_estimate <- otu_table(pop_comparison) + 1
zComposition_clr <- microbiome::transform(zComposition_estimate, "clr")

# add in zCompositions information to new phyloseq object
merged_phylo_counts_zComposition <- pop_comparison
taxa_zComposition <- otu_table(zComposition_clr, taxa_are_rows = TRUE)
otu_table(merged_phylo_counts_zComposition) <- taxa_zComposition

# save the clr-transformed object
save(merged_phylo_counts_zComposition, file=paste0(outputdir,"AllPopsSE_physeq_clr.Rda"))

# Make an ordination plot using euclidean distances
zComposition.ord <- ordinate(merged_phylo_counts_zComposition, method = "PCoA", distance = "euclidean")

# plot by Island
pdf(paste0(outputdir,"PCA_dim1_Island.pdf"), width=7)
plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="SamplePop", axes = 1:2) + scale_colour_manual(values=c(IndonesiaCol,MaliCol,UKCol)) + geom_point(aes(), alpha=0.6, size=6) + theme(plot.title = element_text(hjust = 0, size = 16), axis.text=element_text(size=16), axis.title=element_text(size=17), legend.text = element_text(size=14), legend.title = element_text(size=16), legend.position = "bottom")
dev.off()

pdf(paste0(outputdir,"PCA_dim2to5_Island.pdf"), width=7)
plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="SamplePop", axes = 2:3) + scale_colour_manual(values=c(IndonesiaCol,MaliCol,UKCol)) + geom_point(aes(), alpha=0.6, size=6) + theme(plot.title = element_text(hjust = 0, size = 16), axis.text=element_text(size=16), axis.title=element_text(size=17), legend.text = element_text(size=14), legend.title = element_text(size=16), legend.position = "bottom")
plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="SamplePop", axes = 3:4, label="SampleName") + scale_colour_manual(values=c(IndonesiaCol,MaliCol,UKCol)) + geom_point(aes(), alpha=0.6, size=4) + theme(plot.title = element_text(hjust = 0, size = 12))
plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="SamplePop", axes = 4:5, label="SampleName") + scale_colour_manual(values=c(IndonesiaCol,MaliCol,UKCol)) + geom_point(aes(), alpha=0.6, size=4) + theme(plot.title = element_text(hjust = 0, size = 12))
plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="SamplePop", axes = 5:6, label="SampleName") + scale_colour_manual(values=c(IndonesiaCol,MaliCol,UKCol)) + geom_point(aes(), alpha=0.6, size=4) + theme(plot.title = element_text(hjust = 0, size = 12))
dev.off()

# We can see that the first principal component is separating the three studies apart, while PC2 separates the Indonesian and Malian study from the UK. 
# We can also see that PC3 separates the samples by Plasmodium load. 

logged_phyla_counts = log10(colSums(otu_table(merged_phylo_counts)[grep("Apicomplexa",tax_table(merged_phylo_counts)[,"Phylum"])])+1)
sample_data(merged_phylo_counts_zComposition)[["Apicomplexa"]] = logged_phyla_counts
out.wuf.log <- ordinate(merged_phylo_counts_zComposition, method = "PCoA", distance = "euclidean")

# plot by Apicomplexa
pdf(paste0(outputdir,"PCA_dim1to5_Apicomplexa.pdf"), width=8)
plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="Apicomplexa", axes = 1:2, label="SampleName") + geom_point(aes(), alpha=0.6, size=4) + theme(plot.title = element_text(hjust = 0, size = 12)) + scale_colour_gradientn(colours = colorRampPalette(c("#FFCCBB","#800026"))(10), limits=c(1, 6))
plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="Apicomplexa", axes = 2:3, label="SampleName") + geom_point(aes(), alpha=0.6, size=4) + theme(plot.title = element_text(hjust = 0, size = 12)) + scale_colour_gradientn(colours = colorRampPalette(c("#FFCCBB","#800026"))(10), limits=c(1, 6))
plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="Apicomplexa", axes = 3:4, label="SampleName") + geom_point(aes(), alpha=0.6, size=4) + theme(plot.title = element_text(hjust = 0, size = 12)) + scale_colour_gradientn(colours = colorRampPalette(c("#FFCCBB","#800026"))(10), limits=c(1, 6))
plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="Apicomplexa", axes = 4:5, label="SampleName") + geom_point(aes(), alpha=0.6, size=4) + theme(plot.title = element_text(hjust = 0, size = 12)) + scale_colour_gradientn(colours = colorRampPalette(c("#FFCCBB","#800026"))(10), limits=c(1, 6))
plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="Apicomplexa", axes = 5:6, label="SampleName") + geom_point(aes(), alpha=0.6, size=4) + theme(plot.title = element_text(hjust = 0, size = 12)) + scale_colour_gradientn(colours = colorRampPalette(c("#FFCCBB","#800026"))(10), limits=c(1, 6))
dev.off()

# We can also see that PC4 separates samples by virus load.

logged_phyla_counts = log10(colSums(otu_table(merged_phylo_counts)[grep("Kitrinoviricota",tax_table(merged_phylo_counts)[,"Phylum"])])+1)
sample_data(merged_phylo_counts_zComposition)[["Viruses"]] = logged_phyla_counts
out.wuf.log <- ordinate(merged_phylo_counts_zComposition, method = "PCoA", distance = "euclidean")

# plot by virus load
pdf(paste0(outputdir,"PCA_dim1to5_Viruses.pdf"), width=8)
plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="Viruses", axes = 1:2, label="SampleName") + geom_point(aes(), alpha=0.6, size=4) + theme(plot.title = element_text(hjust = 0, size = 12)) + scale_colour_gradientn(colours = colorRampPalette(c("#78c679","#006837"))(20), limits=c(1, 4))
plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="Viruses", axes = 2:3, label="SampleName") + geom_point(aes(), alpha=0.6, size=4) + theme(plot.title = element_text(hjust = 0, size = 12)) + scale_colour_gradientn(colours = colorRampPalette(c("#78c679","#006837"))(20), limits=c(1, 4))
plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="Viruses", axes = 3:4, label="SampleName") + geom_point(aes(), alpha=0.6, size=4) + theme(plot.title = element_text(hjust = 0, size = 12)) + scale_colour_gradientn(colours = colorRampPalette(c("#78c679","#006837"))(20), limits=c(1, 4))
plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="Viruses", axes = 4:5, label="SampleName") + geom_point(aes(), alpha=0.6, size=4) + theme(plot.title = element_text(hjust = 0, size = 12)) + scale_colour_gradientn(colours = colorRampPalette(c("#78c679","#006837"))(20), limits=c(1, 4))
plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="Viruses", axes = 5:6, label="SampleName") + geom_point(aes(), alpha=0.6, size=4) + theme(plot.title = element_text(hjust = 0, size = 12)) + scale_colour_gradientn(colours = colorRampPalette(c("#78c679","#006837"))(20), limits=c(1, 4))
dev.off()

# Finally, let's see how sample grouping looks like by another clustering method - hierarchical clustering by Euclidean distance. Again, this is done on the CLR-transformed data to correct for sequencing depth.

ps_otu <- data.frame(phyloseq::otu_table(merged_phylo_counts_zComposition))
ps_otu <- t(ps_otu)
bc_dist <- vegan::vegdist(ps_otu, method = "euclidean")
ward <- as.dendrogram(hclust(bc_dist, method = "ward.D2"))
#Provide color codes
meta <- data.frame(phyloseq::sample_data(merged_phylo_counts_zComposition))
colorCode <- c(Indonesia = IndonesiaCol, `Mali` = MaliCol, `UK` = UKCol)
labels_colors(ward) <- colorCode[meta$SamplePop][order.dendrogram(ward)]
#Plot
pdf(paste0(outputdir,"hierarchicalSlutering_IslandLabels.pdf"), width=15)
plot(ward)
dev.off()

# From hierarchical clustering, we can see that samples group by study, which confirms what we saw in the PCA analysis.

##############################
# Relative frequency of taxa #
##############################


# add a new column containing family names and superkingdom
tax_table(merged_phylo_counts)[,"Superkingdom"] = paste(tax_table(merged_phylo_counts)[,"Superkingdom"], tax_table(merged_phylo_counts)[,"Family"], sep="_")
tax_table(merged_phylo_counts)[,"Superkingdom"] <- gsub("Bacteria_$", "Bacteria_unclassified", tax_table(merged_phylo_counts)[,"Superkingdom"])
tax_table(merged_phylo_counts)[,"Superkingdom"] <- gsub("Eukaryota_$", "Eukaryota_unclassified", tax_table(merged_phylo_counts)[,"Superkingdom"])
tax_table(merged_phylo_counts)[,"Superkingdom"] <- gsub("Viruses_$", "Viruses_unclassified", tax_table(merged_phylo_counts)[,"Superkingdom"])

# As pointed out, we have a lot of taxa at the family level, and it would be hard to look over everything at once. Instead, we can focus on the most prevalent taxa and 
# highlight everything else in another colour. Here, I chose to highlight the top 20 taxa, since that's still representative while not being too visually exhausting.

aggregated_phyloCounts <- aggregate_top_taxa(merged_phylo_counts, "Superkingdom", top = 20)
# transform to relative counts
relative_phyloCounts <- microbiome::transform(aggregated_phyloCounts, "compositional")
# Remove weird extra family names added at the end of Superkingdom names
tax_table(relative_phyloCounts)[,"Superkingdom"] <- paste(sapply(strsplit(taxa_names(relative_phyloCounts), "[_.]"), `[`, 1), sapply(strsplit(taxa_names(relative_phyloCounts), "[_.]"), `[`, 2), sep="_")
# Change "Other_NA" to just "Other"
tax_table(relative_phyloCounts)[,"Superkingdom"][grep("Other", taxa_names(relative_phyloCounts))] = "Other"

# Plot
p=plot_bar(relative_phyloCounts, fill = "Superkingdom")

# set colour palette
families=levels(p$data$Superkingdom)
# get number of families in each kingdom
table(sapply(strsplit(families, "[_.]"), `[`, 1))

PaletteArchaea = colorRampPalette(c("#DDCC77"))(1)
PaletteBacteria = colorRampPalette(c("#023858","#74a9cf"))(15)
PaletteEukaryote = colorRampPalette(c("#b20000","#ff4c4c"))(2)
PaletteOther = colorRampPalette(c("black"))(1)
PaletteVirus = colorRampPalette(c("#78c679","#006837"))(2)

Merged_Palette <- c(PaletteArchaea,PaletteBacteria,PaletteEukaryote,PaletteOther,PaletteVirus)

pdf(paste0(outputdir,"relativeTaxa_Compositional.pdf"), width=15)
phyloseq::plot_bar(relative_phyloCounts, fill = "Superkingdom") +
  geom_bar(aes(fill = Superkingdom), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") + theme_bw(base_size = 15) +
  facet_wrap(~ SamplePop, scales = "free") + scale_fill_manual(values=Merged_Palette) + theme(axis.text.x = element_text(angle = 90))
  theme(panel.background = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
dev.off()

# We can see that all populations have a high bacterial load, however in the Malian and Indonesian dataset, there'a also a 
# high abundance of Plasmodium and flavivirus. We can also see that the Malian dataset has a high proportion of reads mapping to archaea. 

# Let's plot this using a circular barplot
data=as.matrix(as.data.frame(otu_table(relative_phyloCounts)))
data=t(data)
data=as.data.frame(data)
data$group=as.character(as.matrix(as.data.frame(phyloseq::sample_data(merged_phylo_counts)[,"SamplePop"])))
data$individual=rownames(data)
data$individual=as.factor(data$individual)
data$group=as.factor(data$group)
 
# Transform data in a tidy format (long format)
data = data %>% gather(key = "observation", value="value", -c(22,23)) 

# Set a number of 'empty bar' to add at the end of each group
empty_bar=3
nObsType=nlevels(as.factor(data$observation))
to_add = data.frame( matrix(NA, empty_bar*nlevels(data$group)*nObsType, ncol(data)) )
colnames(to_add) = colnames(data)
to_add$group=rep(levels(data$group), each=empty_bar*nObsType )
data=rbind(data, to_add)
data=data %>% arrange(group, individual)
data$id=rep( seq(1, nrow(data)/nObsType) , each=nObsType)
 
# prepare a data frame for base lines
base_data=data %>% 
  group_by(group) %>% 
  dplyr::summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data = base_data
grid_data$end = grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start = grid_data$start - 1
 
# Make the plot
p = ggplot(data) +      
  
  # Add the stacked bar
  geom_bar(aes(x=as.factor(id), y=value, fill=observation), stat="identity", alpha=0.8) +
  # Add a valu=100/75/50/25 lines. I do it at the beginning to make sure barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.25, xend = start, yend = 0.25), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.50, xend = start, yend = 0.50), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.75, xend = start, yend = 0.75), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 1, xend = start, yend = 1), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  #annotate("text", x = rep(max(data$id),5), y = c(0, 0.25, 0.50, 0.75, 1), label = c("0", "0.25", "0.50", "0.75", "1") , color="black", size=2 , angle=0, fontface="bold", hjust=1) +

  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=base_data, aes(x = title, y = 1.2, label=group), hjust=c(0,0,0), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE) +
  scale_fill_manual(values=Merged_Palette)

pdf(paste0(outputdir,"relativeTaxa_Compositional_CircularBarplot_AllPops.pdf"), width=15)
p
dev.off()


###################################
# Differential abundance testing #
###################################

# We're interested in testing whether the species composition between populations is significantly different. Visually, we saw that the populations look different, but we can't say that for sure. One way to test this is through differential abundance testing.
# One of the best packages I've found so far to test this is called [Aldex2](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0067019) (ANOVA-Like Differential Expression). Aldex2 works by generatating a posterior probability (128 by default) for an observed instance of a taxon (adding a small prior to deal with zeros), then performing a centered log-ratio transformation on the data as a normalisation step (this deals with uneven library sizes). To identify differentially expressed genes, Aldex2 then performs a significance test using a Wilcoxon rank test, and finally, the probability of the taxon being differentially abundant is adjusted with FDR correction (by Benjamini–Hochberg). 
# Aldex2 corrects for uneven library sizes, so rarefying is not necessary. The only input we need is the data with singletons removed.
# Since Aldex2 works by comparing two datasets, we need to subset the dataset. I'll first compare the Indonesian dataset to the UK dataset, then the Indonesian dataset to the Malian dataset.

# Phylum level ----------------------

# Differential abundance testing
IndoVsUK=subset_samples(pop_comparison, SamplePop != "Mali")
any(taxa_sums(IndoVsUK) == 0)
# remove any 0s from the data
IndoVsUK <- prune_taxa(taxa_sums(IndoVsUK) > 0, IndoVsUK)
taxa_names(IndoVsUK)=make.unique(tax_table(IndoVsUK)[,"Phylum"])
# Run aldex2
aldex2_IndoVsUK <- ALDEx2::aldex(data.frame(phyloseq::otu_table(IndoVsUK)), phyloseq::sample_data(IndoVsUK)$SamplePop, test="t", effect = TRUE)

# Let's now pull out the significant values after getting the ALdex2 ouput. We're interested in the columns 

sig_aldex2_IndoVsUK <- aldex2_IndoVsUK %>%
  rownames_to_column(var = "OTU") %>%
  filter(we.eBH < 0.05) %>%
  arrange(effect, we.eBH) %>%
  dplyr::select(OTU, diff.btw, diff.win, effect, we.ep, we.eBH)

# add in taxonomic information to the significant IndoVsUK object
taxa_info <- data.frame(tax_table(IndoVsUK)[,c("Superkingdom","Kingdom","Phylum")])
taxa_info <- taxa_info %>% rownames_to_column(var = "OTU")
sig_aldex2_IndoVsUK <- left_join(sig_aldex2_IndoVsUK, taxa_info)

# now for the non-significant df
aldex2_IndoVsUK$OTU <- rownames(aldex2_IndoVsUK)
aldex2_IndoVsUK <- left_join(aldex2_IndoVsUK, taxa_info)

# save file
write.table(sig_aldex2_IndoVsUK, file=paste0(outputdir,"sig_aldex2_IndoVsUK_phylum.txt"))

# set significance colours
aldex2_IndoVsUK$threshold <- aldex2_IndoVsUK$we.eBH <= 0.05
aldex2_IndoVsUK$threshold = as.numeric(aldex2_IndoVsUK$threshold) + 1

# add label names
labels = c()
counter = 0
for (name in as.character(aldex2_IndoVsUK$Phylum)){
	counter = counter + 1
	if (name == "unk_p"){
		name = as.character(aldex2_IndoVsUK$Superkingdom[counter])
	}
	labels = c(labels, name)
}
# get superkingdom names
taxa_superkingdom = sapply(strsplit(tax_table(IndoVsUK)[,"Superkingdom"], "[_.]"), `[`, 1)

# plot
pdf(paste0(outputdir,"differentialAbundance_IndoVsUK_phylum.pdf"))
ggplot(aldex2_IndoVsUK) +
  geom_point(aes(x = effect, y = -log10(we.eBH)), color = ifelse(aldex2_IndoVsUK$we.eBH <= 0.05, c("grey","#023858","#800026","grey","#78c679")[as.numeric(as.factor(taxa_superkingdom))],"black"), alpha = 0.65, size=8) +
  #geom_text_repel(aes(x = effect, y = -log10(wi.eBH), label = rownames(aldex2_IndoVsDutch))) +
  geom_text_repel(aes(x = effect, y = -log10(we.eBH), label = ifelse(we.eBH <= 0.05, labels,""))) +
  ggtitle("UK Versus Indonesia") +
  xlab("Effect Size") + 
  ylab("-log10 adjusted p-value") + theme_bw(base_size = 18) +
  theme(legend.position = "bottom",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))
dev.off()

# Malian comparison -------------------

IndoVsMali=subset_samples(pop_comparison, SamplePop != "UK")
any(taxa_sums(IndoVsMali) == 0)
IndoVsMali <- prune_taxa(taxa_sums(IndoVsMali) > 0, IndoVsMali)
taxa_names(IndoVsMali)=make.unique(tax_table(IndoVsMali)[,"Phylum"])
aldex2_IndoVsMali <- ALDEx2::aldex(data.frame(phyloseq::otu_table(IndoVsMali)), phyloseq::sample_data(IndoVsMali)$SamplePop, test="t", effect = TRUE)
sig_aldex2_IndoVsMali <- aldex2_IndoVsMali %>%
  rownames_to_column(var = "OTU") %>%
  filter(we.eBH < 0.05) %>%
  arrange(effect, we.eBH) %>%
  dplyr::select(OTU, diff.btw, diff.win, effect, we.ep, we.eBH)

# add in taxonomic information to the significant IndoVsMali object
taxa_info <- data.frame(tax_table(IndoVsMali)[,c("Superkingdom","Kingdom","Phylum")])
taxa_info <- taxa_info %>% rownames_to_column(var = "OTU")
sig_aldex2_IndoVsMali <- left_join(sig_aldex2_IndoVsMali, taxa_info)
# now for the non-significant df
aldex2_IndoVsMali$OTU <- rownames(aldex2_IndoVsMali)
aldex2_IndoVsMali <- left_join(aldex2_IndoVsMali, taxa_info)

# save file
write.table(sig_aldex2_IndoVsMali, file=paste0(outputdir,"sig_aldex2_IndoVsMali_phylum.txt"))

# set significance colours
aldex2_IndoVsMali$threshold <- aldex2_IndoVsMali$we.eBH <= 0.05
aldex2_IndoVsMali$threshold = as.numeric(aldex2_IndoVsMali$threshold) + 1

# add label names
labels = c()
counter = 0
for (name in as.character(aldex2_IndoVsMali$Phylum)){
	counter = counter + 1
	if (name == "unk_p"){
		name = as.character(aldex2_IndoVsMali$Superkingdom[counter])
	}
	labels = c(labels, name)
}
# get superkingdom names
taxa_superkingdom = sapply(strsplit(tax_table(IndoVsMali)[,"Superkingdom"], "[_.]"), `[`, 1)

# plot
pdf(paste0(outputdir,"differentialAbundance_IndoVsMali_phylum_05threshold.pdf"), width=10)
ggplot(aldex2_IndoVsMali) +
  geom_point(aes(x = effect, y = -log10(we.eBH)), color = ifelse(aldex2_IndoVsMali$we.eBH <= 0.05, c("#DDCC77","#023858","#800026","#78c679")[as.numeric(as.factor(taxa_superkingdom))],"black"), alpha = 0.65, size=8) +
  geom_text_repel(aes(x = effect, y = -log10(we.eBH), label = ifelse(we.eBH <= 0.05, labels,""))) +
  ggtitle("Mali Versus Indonesia") +
  xlab("Effect Size") + 
  ylab("-log10 adjusted p-value") + theme_bw(base_size = 18) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))
dev.off()

pdf(paste0(outputdir,"differentialAbundance_IndoVsMali_phylum_001threshold.pdf"))
ggplot(aldex2_IndoVsMali) +
  geom_point(aes(x = effect, y = -log10(we.eBH)), color = ifelse(aldex2_IndoVsMali$we.eBH <= 0.05, c("#DDCC77","#023858","#800026","#78c679")[as.numeric(as.factor(taxa_superkingdom))],"black"), alpha = 0.65, size=8) +
  geom_text_repel(aes(x = effect, y = -log10(we.eBH), label = ifelse(we.eBH <= 0.001, labels,""))) +
  ggtitle("Mali Versus Indonesia") +
  xlab("Effect Size") + 
  ylab("-log10 adjusted p-value") + theme_bw(base_size = 18) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))
dev.off()

# Korowai only versus Malian population ----------

# first, assign island information in order to separate out the Korowai
Island <- sapply(strsplit(rownames(sample_data(pop_comparison)), "[-.]"), `[`, 1)
sample_data(pop_comparison)$Island = Island
sample_data(pop_comparison)$Island[grepl("SRR1", sample_data(pop_comparison)$Island)] <- "Mali"
sample_data(pop_comparison)$Island[grepl("SRR6", sample_data(pop_comparison)$Island)] <- "UK"

# Differential abundance testing
KorVsMali=subset_samples(pop_comparison, Island != "UK") %>% subset_samples(., Island != "SMB") %>% subset_samples(., Island != "MTW")
any(taxa_sums(KorVsMali) == 0)
# remove any 0s from the data
KorVsMali <- prune_taxa(taxa_sums(KorVsMali) > 0, KorVsMali)
taxa_names(KorVsMali)=make.unique(tax_table(KorVsMali)[,"Phylum"])
# Run aldex2
aldex2_KorVsMali <- ALDEx2::aldex(data.frame(phyloseq::otu_table(KorVsMali)), phyloseq::sample_data(KorVsMali)$SamplePop, test="t", effect = TRUE)

# Let's now pull out the significant values after getting the ALdex2 ouput. We're interested in the columns 

sig_aldex2_KorVsMali <- aldex2_KorVsMali %>%
  rownames_to_column(var = "OTU") %>%
  filter(we.eBH < 0.05) %>%
  arrange(effect, we.eBH) %>%
  dplyr::select(OTU, diff.btw, diff.win, effect, we.ep, we.eBH)

# add in taxonomic information to the significant IndoVsUK object
taxa_info <- data.frame(tax_table(KorVsMali)[,c("Superkingdom","Kingdom","Phylum")])
taxa_info <- taxa_info %>% rownames_to_column(var = "OTU")
sig_aldex2_KorVsMali <- left_join(sig_aldex2_KorVsMali, taxa_info)

# now for the non-significant df
aldex2_KorVsMali$OTU <- rownames(aldex2_KorVsMali)
aldex2_KorVsMali <- left_join(aldex2_KorVsMali, taxa_info)

# save file
write.table(sig_aldex2_KorVsMali, file=paste0(outputdir,"sig_aldex2_KorVsMali_phylum.txt"))

# set significance colours
aldex2_KorVsMali$threshold <- aldex2_KorVsMali$we.eBH <= 0.05
aldex2_KorVsMali$threshold = as.numeric(aldex2_KorVsMali$threshold) + 1

# add label names
labels = c()
counter = 0
for (name in as.character(aldex2_KorVsMali$Phylum)){
	counter = counter + 1
	if (name == "unk_p"){
		name = as.character(aldex2_KorVsMali$Superkingdom[counter])
	}
	labels = c(labels, name)
}
# get superkingdom names
taxa_superkingdom = sapply(strsplit(tax_table(KorVsMali)[,"Superkingdom"], "[_.]"), `[`, 1)

# plot
pdf(paste0(outputdir,"differentialAbundance_KorVsMali_phylum.pdf"))
ggplot(aldex2_KorVsMali ) +
  geom_point(aes(x = effect, y = -log10(we.eBH)), color = ifelse(aldex2_KorVsMali$we.eBH <= 0.05, c("#DDCC77","#023858","#800026","#78c679")[as.numeric(as.factor(taxa_superkingdom))],"black"), alpha = 0.65, size=8) +
  #geom_text_repel(aes(x = effect, y = -log10(wi.eBH), label = rownames(aldex2_IndoVsDutch))) +
  geom_text_repel(aes(x = effect, y = -log10(we.eBH), label = ifelse(we.eBH <= 0.05, labels,""))) +
  ggtitle("Korowai Versus Mali") +
  xlab("Effect Size") + 
  ylab("-log10 adjusted p-value") + theme_bw(base_size = 18) +
  theme(legend.position = "bottom",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) + xlim(c(-1.5,7.5)) + ylim(c(0,17))
dev.off()

# Korowai only versus UK population ----------

# Differential abundance testing
KorVsUK=subset_samples(pop_comparison, Island != "Mali") %>% subset_samples(., Island != "SMB") %>% subset_samples(., Island != "MTW")
any(taxa_sums(KorVsUK) == 0)
# remove any 0s from the data
KorVsUK <- prune_taxa(taxa_sums(KorVsUK) > 0, KorVsUK)
taxa_names(KorVsUK)=make.unique(tax_table(KorVsUK)[,"Phylum"])
# Run aldex2
aldex2_KorVsUK <- ALDEx2::aldex(data.frame(phyloseq::otu_table(KorVsUK)), phyloseq::sample_data(KorVsUK)$SamplePop, test="t", effect = TRUE)

# Let's now pull out the significant values after getting the ALdex2 ouput. We're interested in the columns 

sig_aldex2_KorVsUK <- aldex2_KorVsUK %>%
  rownames_to_column(var = "OTU") %>%
  filter(we.eBH < 0.05) %>%
  arrange(effect, we.eBH) %>%
  dplyr::select(OTU, diff.btw, diff.win, effect, we.ep, we.eBH)

# add in taxonomic information to the significant IndoVsUK object
taxa_info <- data.frame(tax_table(KorVsUK)[,c("Superkingdom","Kingdom","Phylum")])
taxa_info <- taxa_info %>% rownames_to_column(var = "OTU")
sig_aldex2_KorVsUK <- left_join(sig_aldex2_KorVsUK, taxa_info)

# now for the non-significant df
aldex2_KorVsUK$OTU <- rownames(aldex2_KorVsUK)
aldex2_KorVsUK <- left_join(aldex2_KorVsUK, taxa_info)

# save file
write.table(sig_aldex2_KorVsUK, file=paste0(outputdir,"sig_aldex2_KorVsUK_phylum.txt"))

# set significance colours
aldex2_KorVsUK$threshold <- aldex2_KorVsUK$we.eBH <= 0.05
aldex2_KorVsUK$threshold = as.numeric(aldex2_KorVsUK$threshold) + 1

# add label names
labels = c()
counter = 0
for (name in as.character(aldex2_KorVsUK$Phylum)){
	counter = counter + 1
	if (name == "unk_p"){
		name = as.character(aldex2_KorVsUK$Superkingdom[counter])
	}
	labels = c(labels, name)
}
# get superkingdom names
taxa_superkingdom = sapply(strsplit(tax_table(KorVsUK)[,"Superkingdom"], "[_.]"), `[`, 1)

# plot
pdf(paste0(outputdir,"differentialAbundance_KorVsUK_phylum.pdf"))
ggplot(aldex2_KorVsUK) +
  geom_point(aes(x = effect, y = -log10(we.eBH)), color = ifelse(aldex2_KorVsUK$we.eBH <= 0.05, c("#DDCC77","#023858","#800026","#78c679")[as.numeric(as.factor(taxa_superkingdom))],"black"), alpha = 0.65, size=8) +
  #geom_text_repel(aes(x = effect, y = -log10(wi.eBH), label = rownames(aldex2_IndoVsDutch))) +
  geom_text_repel(aes(x = effect, y = -log10(we.eBH), label = ifelse(we.eBH <= 0.05, labels,""))) +
  ggtitle("Korowai Versus UK") +
  xlab("Effect Size") + 
  ylab("-log10 adjusted p-value") + theme_bw(base_size = 18) +
  theme(legend.position = "bottom",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) + xlim(c(-1.5,7.5)) + ylim(c(0,17))
dev.off()

# Now test this at the Family level --------------

pop_comparison_family <- merged_phylo_counts %>%
  tax_glom("Family")

# Differential abundance testing
IndoVsUK=subset_samples(pop_comparison_family, SamplePop != "Mali")
any(taxa_sums(IndoVsUK) == 0)
# remove any 0s from the data
IndoVsUK <- prune_taxa(taxa_sums(IndoVsUK) > 0, IndoVsUK)
taxa_names(IndoVsUK)=make.unique(tax_table(IndoVsUK)[,"Family"])
# Run aldex2
aldex2_IndoVsUK <- ALDEx2::aldex(data.frame(phyloseq::otu_table(IndoVsUK)), phyloseq::sample_data(IndoVsUK)$SamplePop, test="t", effect = TRUE)

sig_aldex2_IndoVsUK <- aldex2_IndoVsUK %>%
  rownames_to_column(var = "OTU") %>%
  filter(we.eBH < 0.05) %>%
  arrange(effect, we.eBH) %>%
  dplyr::select(OTU, diff.btw, diff.win, effect, we.ep, we.eBH)

# add in taxonomic information to the significant IndoVsUK object
taxa_info <- data.frame(tax_table(IndoVsUK)[,c("Superkingdom","Kingdom","Phylum","Class","Order","Family")])
taxa_info <- taxa_info %>% rownames_to_column(var = "OTU")
sig_aldex2_IndoVsUK <- left_join(sig_aldex2_IndoVsUK, taxa_info)

# now for the non-significant df
aldex2_IndoVsUK$OTU <- rownames(aldex2_IndoVsUK)
aldex2_IndoVsUK <- left_join(aldex2_IndoVsUK, taxa_info)

# set significance colours
aldex2_IndoVsUK$threshold <- aldex2_IndoVsUK$we.eBH <= 0.05
aldex2_IndoVsUK$threshold = as.numeric(aldex2_IndoVsUK$threshold) + 1

# add label names
labels = c()
counter = 0
for (name in as.character(aldex2_IndoVsUK$Family)){
	counter = counter + 1
	if (grepl("unk", name)){
		name = as.character(aldex2_IndoVsUK$Superkingdom[counter])
	}
	labels = c(labels, name)
}
# get superkingdom names
taxa_superkingdom = sapply(strsplit(tax_table(IndoVsUK)[,"Superkingdom"], "[_.]"), `[`, 1)

# plot
pdf(paste0(outputdir,"differentialAbundance_IndoVsUK_Family_05threshold.pdf"))
ggplot(aldex2_IndoVsUK) +
  geom_point(aes(x = effect, y = -log10(we.eBH)), color = ifelse(aldex2_IndoVsUK$we.eBH <= 0.05, c("grey","#023858","#800026","grey","#78c679")[as.numeric(as.factor(taxa_superkingdom))],"black"), alpha = 0.65, size=8) +
  #geom_text_repel(aes(x = effect, y = -log10(wi.eBH), label = rownames(aldex2_IndoVsDutch))) +
  geom_text_repel(aes(x = effect, y = -log10(we.eBH), label = ifelse(we.eBH <= 0.05, labels,""))) +
  ggtitle("UK Versus Indonesia") +
  xlab("Effect Size") + 
  ylab("-log10 adjusted p-value") + theme_bw(base_size = 18) +
  theme(legend.position = "bottom",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))
dev.off()

IndoVsMali=subset_samples(pop_comparison_family, SamplePop != "UK")
any(taxa_sums(IndoVsMali) == 0)
IndoVsMali <- prune_taxa(taxa_sums(IndoVsMali) > 0, IndoVsMali)
taxa_names(IndoVsMali)=make.unique(tax_table(IndoVsMali)[,"Family"])
aldex2_IndoVsMali <- ALDEx2::aldex(data.frame(phyloseq::otu_table(IndoVsMali)), phyloseq::sample_data(IndoVsMali)$SamplePop, test="t", effect = TRUE)
sig_aldex2_IndoVsMali <- aldex2_IndoVsMali %>%
  rownames_to_column(var = "OTU") %>%
  filter(we.eBH < 0.05) %>%
  arrange(effect, we.eBH) %>%
  dplyr::select(OTU, diff.btw, diff.win, effect, we.ep, we.eBH)

# add in taxonomic information to the significant IndoVsMali object
taxa_info <- data.frame(tax_table(IndoVsMali)[,c("Superkingdom","Kingdom","Phylum","Class","Order","Family")])
taxa_info <- taxa_info %>% rownames_to_column(var = "OTU")
sig_aldex2_IndoVsMali <- left_join(sig_aldex2_IndoVsMali, taxa_info)

# now for the non-significant df
aldex2_IndoVsMali$OTU <- rownames(aldex2_IndoVsMali)
aldex2_IndoVsMali <- left_join(aldex2_IndoVsMali, taxa_info)

# set significance colours
aldex2_IndoVsMali$threshold <- aldex2_IndoVsMali$we.eBH <= 0.05
aldex2_IndoVsMali$threshold = as.numeric(aldex2_IndoVsMali$threshold) + 1

# add label names
labels = c()
counter = 0
for (name in as.character(aldex2_IndoVsMali$Family)){
	counter = counter + 1
	if (grepl("unk", name)){
		name = as.character(aldex2_IndoVsMali$Superkingdom[counter])
	}
	labels = c(labels, name)
}
# get superkingdom names
taxa_superkingdom = sapply(strsplit(tax_table(IndoVsMali)[,"Superkingdom"], "[_.]"), `[`, 1)

# plot
pdf(paste0(outputdir,"differentialAbundance_IndoVsMali_Family_05threshold.pdf"))
ggplot(aldex2_IndoVsMali) +
  geom_point(aes(x = effect, y = -log10(we.eBH)), color = ifelse(aldex2_IndoVsMali$we.eBH <= 0.05, c("#DDCC77","#023858","#800026","grey","#78c679")[as.numeric(as.factor(taxa_superkingdom))],"black"), alpha = 0.65, size=8) +
  geom_text_repel(aes(x = effect, y = -log10(we.eBH), label = ifelse(we.eBH <= 0.05, labels,""))) +
  ggtitle("Mali Versus Indonesia") +
  xlab("Effect Size") + 
  ylab("-log10 adjusted p-value") + theme_bw(base_size = 18) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))
dev.off()

pdf(paste0(outputdir,"differentialAbundance_IndoVsMali_Family_001threshold.pdf"))
ggplot(aldex2_IndoVsMali) +
  geom_point(aes(x = effect, y = -log10(we.eBH)), color = ifelse(aldex2_IndoVsMali$we.eBH <= 0.05, c("#DDCC77","#023858","#800026","grey","#78c679")[as.numeric(as.factor(taxa_superkingdom))],"black"), alpha = 0.65, size=8) +
  geom_text_repel(aes(x = effect, y = -log10(we.eBH), label = ifelse(we.eBH <= 0.001, labels,""))) +
  ggtitle("Mali Versus Indonesia") +
  xlab("Effect Size") + 
  ylab("-log10 adjusted p-value") + theme_bw(base_size = 18) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))
dev.off()


###################
# Alpha diversity #
###################

# Current tools to estimate alpha diversity either underestimate richness or underestimate uncertainty, however DivNet is a package that adresses these problems. Divnet is a method for estimating within- and between-community diversity in ecosystems where taxa interact via an ecological network. It accounts for differences in sequencing depth and estimates the number of missing species based on the sequence depth and
# number of rare taxa in the data.

# To use DivNet, you need unsubsampled data without removing singletons.

# Get phyloseq object without singletons by merging original phyloseq objects
merged_phylo_counts_withSingletons <- merge_phyloseq(AllREadsSE_Indo_Counts_physeq, Mali_Counts_physeq, European_Counts_physeq)
# Remove Viridiplantae and Metazoa
merged_phylo_counts_withSingletons <- subset_taxa(merged_phylo_counts_withSingletons, (Kingdom!="Viridiplantae"))
merged_phylo_counts_withSingletons <- subset_taxa(merged_phylo_counts_withSingletons, (Kingdom!="Metazoa"))
merged_phylo_counts_withSingletons <- subset_taxa(merged_phylo_counts_withSingletons, (Superkingdom!="unk_sk"))
# remove any empty rows
merged_phylo_counts_withSingletons <- prune_taxa(taxa_sums(merged_phylo_counts_withSingletons) > 0, merged_phylo_counts_withSingletons)

# Now that we have data without singletons, we now need to merge our data at a specified taxonomic level. DivNet is computationally expensive, and therefore a higher level is much, much faster.
# We'll therefore test how our groups look like at the Phylum level. Then, we'll run DivNet without specifying any hypothesis testing.

# comparing diversity at the phylum level
pop_comparison <- merged_phylo_counts_withSingletons %>%
  tax_glom("Phylum")

# If we don't change the sample names here from hyphens to periods, we'll get an error later
sample_names(pop_comparison) <- gsub("\\-", ".", sample_names(pop_comparison))

# Run divnet without specifying any hypothesis testing
dv_pop_comparison <- divnet(pop_comparison, ncores = 4)

# DivNet will output an object with estimates for multiple different alpha (and beta) diversity measures (we'll get to the beta diversity estimates later).
# The Shannon and Simpson index are two popular alpha diversity indices to measure species richness. For Shannon diversity, the importance of rare taxa are downeighted, since they do not play a large role in the community or they could potentially be due to error. For this reason, the Shannon index is one of the most popular alpha diversity metrics.
# To interpret the index, a higher Shannon index means higher diversity, whereas a lowed index number means lower diversity.
# Let's take out the Shannon diversity metric from DivNet and plot it.

# Now let's plot the results of shannon and Simpson diversity
summary_df_shannon <- as.data.frame(dv_pop_comparison$shannon %>%
  summary %>%
  add_column("SampleNames" = pop_comparison %>% otu_table %>% sample_names) %>%
  add_column("SamplePop" =  pop_comparison %>% sample_data %>% .[,"SamplePop"] %>% as.matrix(.) %>% .[,1] %>% unname(.)))

# save shannon diversity summary
write.table(summary_df_shannon, file=paste0(outputdir,"summary_df_shannon_noHypothesisTesting.txt"))

pdf(paste0(outputdir,"shannonDiversity_noHypothesisTesting.pdf"), width=7)
ggplot(summary_df_shannon, aes(y = estimate, x = SamplePop, fill = SamplePop)) + geom_violin(alpha=0.7) + 
  geom_boxplot(width=0.08, outlier.color = NA) + 
  scale_fill_manual(values=c(IndonesiaCol,MaliCol,UKCol)) + ggtitle("Shannon Diversity") +
  ylab("Estimate of Shannon Diversity") + xlab("") + theme_bw(base_size = 16) + theme(legend.position="bottom")
dev.off()


# We can see that the Malian samples, on average, have the highest estimates of Shannon alpha diversity, followed by the Indonesian population, then the UK population.
# We can also plot each individual sample, along with their standard deviation (another cool, and imprortant feature that DivNet calculates and uses in their hypothesis testing). 

pdf(paste0(outputdir,"shannonDiversity_noHypothesisTesting_indivPopulations.pdf"), width=10)
plot(dv_pop_comparison$shannon, pop_comparison, col = "SamplePop") + scale_colour_manual(values=c(IndonesiaCol,MaliCol,UKCol)) + 
	theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
dev.off()

# We can see that the Indonesian samples have the highest amount of SE around their estimates, but this makes sense given that this is the population with the lowest library size.
# Now let's see how the population diversity looks like when we use the Simpson diversity index. The Simpson diversity index is a similarity index where the higher the value, the lower the diversity. It measures the probability that two individuals randomly selected from a sample will belong to the same species. With this index, 0 represents infinite diversity and 1, no diversity.

summary_df_simpson <- as.data.frame(dv_pop_comparison$simpson %>%
  summary %>%
  add_column("SampleNames" = pop_comparison %>% otu_table %>% sample_names) %>%
  add_column("SamplePop" =  pop_comparison %>% sample_data %>% .[,"SamplePop"] %>% as.matrix(.) %>% .[,1] %>% unname(.)))

write.table(summary_df_simpson, file=paste0(outputdir,"summary_df_simpson_noHypothesisTesting.txt"))


pdf(paste0(outputdir,"simpsonDiversity_noHypothesisTesting.pdf"), width=10)
ggplot(summary_df_simpson, aes(y = estimate, x = SamplePop, fill = SamplePop)) + geom_violin(alpha=0.7) + 
  geom_boxplot(width=0.08, outlier.color = NA) +
  scale_fill_manual(values=c(IndonesiaCol,MaliCol,UKCol)) + ggtitle("Simpson's Diversity Index") +
  ylab("Estimate of Simpson Diversity") + xlab("")
dev.off()

# Again, we can see that the same trend holds true - the Malian population has the highest diversity (remember, and index of 0 equates to infinite diversity), while the UK population has the lowest diversity.
# This is how the diversity looks like with SE included for each sample.

pdf(paste0(outputdir,"simpsonDiversity_noHypothesisTesting_indivPopulations.pdf"), width=10)
plot(dv_pop_comparison$simpson, pop_comparison, col = "SamplePop") + scale_colour_manual(values=c(IndonesiaCol,MaliCol,UKCol)) +
	theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
dev.off()

# Since a larger Simpson index value equates to a lower diversity index, many people find this confusing and not
# very intuitive. Therefore, the inverse Simpsone Index, or 1 - Simpson Index, is also commonly used.
# Let's plot that now. 

# Subtract the Simpson estimate from one
summary_df_simpson$estimate = 1-summary_df_simpson$estimate

# Plot
pdf(paste0(outputdir,"InverseSimpsonDiversity_noHypothesisTesting.pdf"), width=7)
ggplot(summary_df_simpson, aes(y = estimate, x = SamplePop, fill = SamplePop)) + geom_violin(alpha=0.7) + 
  geom_boxplot(width=0.08, outlier.color = NA) +
  scale_fill_manual(values=c(IndonesiaCol,MaliCol,UKCol)) + ggtitle("Inverse Simpson's Diversity Index") +
  ylab("Estimate of Inverse Simpson Diversity") + xlab("") + theme_bw(base_size = 16) + theme(legend.position="bottom")
dev.off()

# Now let's test the hypothesis that the diversity is different between islands, which is really what we're asking - are the taxa communities between our groups different? This is holding the assumption that the populations are part of a shared ecosystem that has similarities in, for example, pathogen load. We are now estimating the diversity of island/population being an ecosystem, so we're focusing on the ecosystem, not just the samples. If we wanted to reproduce this result, it is better to focus on the populations that the samples come from, not just the samples themselves.
# Let's test this first using the Shannon diversity index. 

# test the hypothesis that the diversity is differnet between islands
dv_pop_comparison_cov <- pop_comparison %>%
  divnet(X = "SamplePop", ncores = 8)

# Plot the results for each individual
pdf(paste0(outputdir,"ShannonDiversity_perIndividual_withHypothesisTesting.pdf"), width=10)
plot(dv_pop_comparison_cov$shannon, pop_comparison, col = "SamplePop") + scale_colour_manual(values=c(IndonesiaCol,MaliCol,UKCol))
dev.off()

# We can see that, as a population, the Malian samples have the highest diversity, followed by the Indonesian samples, then the UK samples.
# Let's test this hypothesis formally.

# test that these populations are actually different
diversity_test_shannon = testDiversity(dv_pop_comparison_cov, "shannon")

# write to file
write.table(diversity_test_shannon, file=paste0(outputdir,"HypothesisTest_Shannon.txt"))

# The result tells us that the mean Shannon diversity in the Indonesian population at the Phylum-level is 0.91, and it is significantly higher by 0.7 orders, on average, in the Malian population. We can also see that the UK population is significantly lower than the Indonesian population by 0.4 orders, on average.
# Let's do the same thing for Simpson diversity.

pdf(paste0(outputdir,"SimpsonDiversity_perIndividual_withHypothesisTesting.pdf"), width=10)
plot(dv_pop_comparison_cov$simpson, pop_comparison, col = "SamplePop") + scale_colour_manual(values=c(IndonesiaCol,MaliCol,UKCol))
dev.off()

# With the Simpson index, we see that the UK population has the highest Simpson index (i.e., lowest diversity),
# followed by the Indonesian population, then the Malian population.
# Again, let's test this formally.

diversity_test_simpson = testDiversity(dv_pop_comparison_cov, "simpson")
# write to file
write.table(diversity_test_simpson, file=paste0(outputdir,"HypothesisTest_Shannon.txt"))

# The result tells us that the mean Simpson diversity index in the Indonesian population at the Phylum-level is 0.6, and it is significantly lower by 0.35 orders, on average, in the Malian population. We can also see that the UK population is significantly higher than the Indonesian population by 0.09 orders, on average.

##################
# Beta Diversity #
##################

# Beta diversity is a measure of dissimilarity metric between samples to compare differneces in species composition. It's helpful to know not only how taxonomically/pathogenically rich
# each sample is, but also to see differences in samples and populations.
# There are multiple beta diversity measures to use, including Bray-curtis dissimilarity (based on abundances), Jaccard distance (based on presence or absence), Euclidean distance, and Unifrac (based on sequence distances using a phylogenetic tree).
# The Bray–Curtis dissimilarity metric is probably the most popular beta diversity metric and is bounded between 0 and 1, where 0 means the two sites have the same composition (all species are shared), and 1 means the two sites do not share any species.
# Unlike alpha diversity, beta diversity is not as sensitive to singletons, and it has even been suggested that 'using simple proportions' (i.e., relative abundance) on non-rarefied data [is fine](https://github.com/joey711/phyloseq/issues/470). I haven't found enough information on this one way or another (it seems like people are far more opinionated on alpha diversity than beta diversity), but I belive sticking with clr-transformed data is essential. 
# You can't use negative values (which CLR has) for Bray-Curtis dissimilarity, so therefore I'll use Euclidean distances on CLR-transformed data using the ordinate function from Phyloseq. Then, I'll compare this to Bray-Curtis dissimilarity and Euclidean distances calculated in DivNet.

# Using Euclidean distances on the CLR-transformed data, we can see that the Malian data is the least similar from the Indonesian and UK data. We can also see that the UK samples cluster inside the Indonesian samples (I'm not sure how to interpret this).
# Now let's use our results from DivNet, which corrects for sequencing depth, and see how this compares.
# For the first bray-curtis dissimilarity analysis, we'll look at all of our samples individually, then we'll look
# at how bray-curtis dissimilarity looks like between islands (with hypothesis testing).

# First, let's look at Bray-curtis dissimilarity at the individual sample level
bray_est <- simplifyBeta(dv_pop_comparison, pop_comparison, "bray-curtis", "SamplePop")
write.table(bray_est, file=paste0(outputdir,"bray_diversity_estimate.txt"))

# add in group comparisons and plot
bray_est$group=paste(bray_est$Covar1,bray_est$Covar2,sep="_")
bray_est$Int <- interaction(bray_est$Covar1, bray_est$Covar2)
bray_est$Int=factor(bray_est$Int, levels = c("UK.Mali", "Mali.Indonesia", "Mali.Mali", "Indonesia.Indonesia", "UK.Indonesia", "UK.UK"))

pdf(paste0(outputdir,"BrayDiversity.pdf"), width=15)
ggplot(data= bray_est, aes(x= Int, y=beta_est, fill=group)) + geom_violin(alpha=0.7) + geom_boxplot(width=0.1) + xlab("Population Comparisons") + theme_bw(base_size = 18) +
  theme(legend.position="none") + ggtitle("Bray-Curtis Distance Estimate") + 
  ylab("Bray-Curtis Distance") + scale_fill_manual(values = c("#228833","#EE7733",MaliCol,IndonesiaCol,"#EE3377",UKCol))
dev.off()

# pdf(paste0(outputdir,"BrayDiversity.pdf"), width=15)
# ggplot(bray_est, aes(x = interaction(Covar1, Covar2), y = beta_est, fill=group)) +
#   geom_violin(alpha=0.7) + geom_boxplot(width=0.1) + xlab("Population Comparisons") + theme_bw(base_size = 18) +
#   theme(legend.position="none") + ggtitle("Bray-Curtis Distance Estimate") +
#   ylab("Bray-Curtis Distance") + scale_fill_manual(values = c(KorowaiCol,"#86509c",MentawaiCol,"#f97702","#019680",SumbaCol))
# dev.off()

# From DivNet, again we can see that the greatest dissimilarity is between the Malian population, particularly the UK versus Malia. Unsurpiringly, the least dissimilar samples are comparing populations to themselves (UK samples to UK samples and Indonesian samples to Indonesian samples), with the exception of the Malian samples. When comparing beta diversity within the Malian samples, we can see that the mean Bray-Curtis dissimilarity estimate is nearly 0.5, which is even higher than UK samples compared to Indonesian samples. We also see quite a but of spread in distance estimates when comparing populations to the Indonesian or Malian samples.
# Now let's test beta diversity in DivNet using Euclidean distance.

# First, let's look at Bray-curtis dissimilarity at the individual sample level
bray_est_eucl <- simplifyBeta(dv_pop_comparison, pop_comparison, "euclidean", "SamplePop")
write.table(bray_est_eucl, file=paste0(outputdir,"euclidean_distance_estimate.txt"))

# add in group comparisons and plot
bray_est_eucl$group <- paste(bray_est_eucl$Covar1, bray_est_eucl$Covar2,sep="_")

pdf(paste0(outputdir,"EuclideanDistance.pdf"), width=10)
ggplot(bray_est_eucl, aes(x = interaction(Covar1, Covar2), y = beta_est, fill=group)) +
  geom_violin(alpha=0.7) + geom_boxplot(width=0.1) + xlab("Population Comparisons") + 
  theme(legend.position="none") + ggtitle("Euclidean Distance Estimate") +
  ylab("Euclidean Distance")
dev.off()

# We can see very similar trends in using Euclidean distance, compared to the Bray-Curtis dissimilarity metric. 
# Now let's see how this looks like for island-level comparisons.

# Bray-Curtis dissimilarity
pdf(paste0(outputdir,"BrayDiversity_WIthHypothesisTesting.pdf"), width=10)
simplifyBeta(dv_pop_comparison_cov, pop_comparison, "bray-curtis", "SamplePop") %>%
  ggplot(aes(x = interaction(Covar1, Covar2), 
             y = beta_est)) +
  geom_point() +
  geom_linerange(aes(ymin = lower, ymax = upper)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") + ylab("Estimates of Bray-Curtis distance")
dev.off()

# Euclidean distance
pdf(paste0(outputdir,"EuclideanDistance_WIthHypothesisTesting.pdf"), width=10)
simplifyBeta(dv_pop_comparison_cov, pop_comparison, "euclidean", "SamplePop") %>% filter(beta_var != 0) %>% 
  ggplot(aes(x = interaction(Covar1, Covar2), 
             y = beta_est)) +
  geom_point() +
  geom_linerange(aes(ymin = lower, ymax = upper)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") + ylab("Estimates of Bray-Curtis distance")
dev.off()

# This confirms that indeeed, the Malian population is the most dissimilar from the other two populations.