# Indonesian Sample SE Analysis by KSB
# 2020-10-1

# This part of the study will analyse differences unmapped reads in our Indonesian samples. 
# The aim of this analysis is therefore to see what pathogen can be identified in whole blood in our Indonesian samples and if it differs between islands. We will test this by looking at sample clustering, relative abundance of taxa, differential abundance testing, and diversity estimates. 

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
library(factoextra)
library(taxize)

# set up directories
refdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/EpiStudy/"
filteringDir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/ControlSampleComparison/"
outputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/SEIndonesianSamples/"
batchInfodir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/dataPreprocessing/"

# set ggplot colour theme to white
theme_set(theme_bw())

# Set up colour schemes
KorowaiCol="#F21A00"
MentawaiCol="#3B9AB2"
SumbaCol="#EBCC2A"

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
save(removeReplicates, file=paste0(outputdir,"removeReplicates.Rda"))
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
save(AllREadsSE_Indo_Counts_physeq, file=paste0(outputdir,"AllREadsSE_Indo_Counts_physeq_unfiltered.Rda"))

# We now have a phyloseq object of 117 samples and 2,763 taxa.

###################
# Data processing #
###################

## Removing singletons from the data ---------

# The easiest way to get rid of some error in your data is to throw out any count information below some threshold. Oddly, in microbiomics, there's no set thresholding for this. In the end, it's really a compromise between accuracy and keeping rare taxa. What *is* decided is that filtering out at least singletons is standard, since these are regarded as sources of error or contamination. Some resources I found that were helpful on this can be found [here](http://drive5.com/usearch/manual/singletons.html) and [here](https://forum.qiime2.org/t/do-you-guys-still-remove-singletons-or-doubletons-these-days/7138/2).  
# We can see that there's a high number of low counts within the data.

# histogram of data
pdf(paste0(outputdir,"seqDepthHistpgram.pdf"))
ggplot(meta(AllREadsSE_Indo_Counts_physeq)) + geom_histogram(aes(x = SeqDepth_Prefilter), alpha= 0.6, bins=100)
dev.off()

# For the Indonesian dataset, the sequencing depth is variable and quite low in some samples. Therefore, pushing this threshold up too high will eliminate rare taxa, especially given that we didn't have a high library size to begin with.
# Although our starting library size is small, let's explore the data a bit by looking at the effect of removing singletons.

# Separate species' abundances and taxonomy columns
rarecurve_counts <- otu_table(AllREadsSE_Indo_Counts_physeq)
col <- c(rep(KorowaiCol,sum(sample_data(AllREadsSE_Indo_Counts_physeq)[,"SamplePop"]=="MPI")),rep(MentawaiCol,sum(sample_data(AllREadsSE_Indo_Counts_physeq)[,"SamplePop"]=="MTW")),rep(SumbaCol,sum(sample_data(AllREadsSE_Indo_Counts_physeq)[,"SamplePop"]=="SMB")))
# Try with different filtering thresholds:
pdf(paste0(outputdir,"rarefaction.pdf"))
for (i in c(1,5)){
 	rarecurve_counts[rarecurve_counts<=i]<-0
	rarecurve(t(otu_table(rarecurve_counts, taxa_are_rows = TRUE)), step=200, col=col,label=F, xlab="Counts",ylab="Number Species",main=paste("Removing Reads",i,"And Below",sep=" "),xlim=c(0,50000))
}
dev.off()

# From the rarefaction curves, we can see that most samples qre starting to asymptote at or under 10,000 reads (however, it also looks like we need a higher sequencing depth, as there isn't a definitive asymptote, meaning we didn't capture everything). 
# We can also see that when removing reads 5 and below (shown on the right hand side), the curve asymptotes much sooner. This is because you need more reads to detect rare species, and therefore removing reads 5 and below eliminates some of these rare species.
# As mentioned before, because our starting read depth is small, we will stick with removing singletons. We will also add the sequencing depth information to the phyloseq object to keep track of oir library size after filtering.

# assign unfiltered phyloseq object to a new object for downstream use
AllREadsSE_Indo_Counts_physeq_withSingletons <- AllREadsSE_Indo_Counts_physeq
# Filter out singletons
otu_table(AllREadsSE_Indo_Counts_physeq)[otu_table(AllREadsSE_Indo_Counts_physeq)<=1]<-0
AllREadsSE_Indo_Counts_physeq <- prune_taxa(taxa_sums(AllREadsSE_Indo_Counts_physeq) > 0, AllREadsSE_Indo_Counts_physeq)
# add sequencing depth information after filtering out singletons
SeqDepth_noSingletons = colSums(otu_table(AllREadsSE_Indo_Counts_physeq))
sample_data(AllREadsSE_Indo_Counts_physeq)$SeqDepth_noSingletons = SeqDepth_noSingletons

save(AllREadsSE_Indo_Counts_physeq, file=paste0(outputdir,"AllREadsSE_Indo_Counts_physeq_noSingletons.Rda"))

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

# save the data
save(AllREadsSE_Indo_Counts_physeq, file=paste0(outputdir,"AllREadsSE_Indo_Counts_physeq_filtered.Rda"))

# Also save the compositional data
relative_phyloCounts_compositional <- microbiome::transform(AllREadsSE_Indo_Counts_physeq, "compositional")
save(relative_phyloCounts_compositional, file=paste0(outputdir,"relative_phyloCounts_compositional_filtered.Rda"))
relative_phyloCounts_clr <- microbiome::transform(AllREadsSE_Indo_Counts_physeq, "clr")
save(relative_phyloCounts_clr, file=paste0(outputdir,"relative_phyloCounts_clr_filtered.Rda"))

## Summarising the data ---------------

# Now that we've done all of the filtering, we can plot the final library sizes.

# barplot of library sizes
pdf(paste0(outputdir,"librarySizes.pdf"))
ggplot(meta(AllREadsSE_Indo_Counts_physeq), aes(SampleName, SeqDepth_noUnkSk)) + geom_bar(stat = "identity", aes(fill = SamplePop)) +
scale_fill_manual(values = c(KorowaiCol,MentawaiCol,SumbaCol)) + rotate_x_text()
dev.off()

# We can see that the library sizes are highly uneven, with some samples having a  very high library depth (many in the Korowai), and many having a very low library depth. 
# The final step us is to summarise the data and see how many reads we lost at each filtering step. We'll show this one the log scale since there is a large variance in depth between samples. 

FilteringSummary = sample_data(AllREadsSE_Indo_Counts_physeq)[,c("SamplePop","SeqDepth_Prefilter","SeqDepth_noSingletons","SeqDepth_noViridiplantae","SeqDepth_noChordata","SeqDepth_noMetazoa", "SeqDepth_noUnkSk")]
# save filtering summary
write.table(FilteringSummary, file=paste0(outputdir,"FilteringSummary.txt"))
# melt df and plot
melted_FilteringSummary = melt(FilteringSummary)

pdf(paste0(outputdir,"FilteringSummary.pdf"))
ggplot(melted_FilteringSummary, aes(x=variable, y=log(value))) +
  geom_violin(alpha=0.6, fill="#AA3377") + theme_bw(base_size = 18) + ylab("Spearman pairwise correlation") +
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90)) +
  geom_boxplot(color="black",width=0.2, alpha = 0.7)
dev.off()

# We can see that most of the reads are removed when removing Chordates.

# Finally, fill in missing Phlya names using the package 'taxize'
index=which(tax_table(AllREadsSE_Indo_Counts_physeq)[,"Phylum"]=="unk_p")
specieslist <- c(unname(tax_table(AllREadsSE_Indo_Counts_physeq)[index,"Family"]))
a=tax_name(query = c(specieslist), get = c("phylum"), db = "ncbi")
a[is.na(a)]="unk_p"
tax_table(AllREadsSE_Indo_Counts_physeq)[index,"Phylum"]=a[,"phylum"]

################
# Data summary #
################

AllREadsSE_Indo_Counts_physeq
# otu_table()   OTU Table:         [ 1390 taxa and 117 samples ]
# sample_data() Sample Data:       [ 117 samples by 10 sample variables ]
# tax_table()   Taxonomy Table:    [ 1390 taxa by 8 taxonomic ranks ]

data_summary_allPops_Family <- AllREadsSE_Indo_Counts_physeq %>%
  tax_glom("Family")

taxa_names(data_summary_allPops_Family)=make.unique(paste(tax_table(data_summary_allPops_Family)[,c("Superkingdom")], tax_table(data_summary_allPops_Family)[,c("Phylum")], tax_table(data_summary_allPops_Family)[,c("Class")], tax_table(data_summary_allPops_Family)[,c("Order")], tax_table(data_summary_allPops_Family)[,c("Family")], sep="_"))

round(rev(sort(rowSums(otu_table(data_summary_allPops_Family))/sum(rowSums(otu_table(data_summary_allPops_Family)))))*100, 3)

data_summary_allPops_Family
# otu_table()   OTU Table:         [ 271 taxa and 117 samples ]
# sample_data() Sample Data:       [ 117 samples by 10 sample variables ]
# tax_table()   Taxonomy Table:    [ 271 taxa by 8 taxonomic ranks ]

data_summary_allPops <- AllREadsSE_Indo_Counts_physeq %>%
  tax_glom("Phylum")

taxa_names(data_summary_allPops)=make.unique(paste(tax_table(data_summary_allPops)[,c("Superkingdom")], tax_table(data_summary_allPops)[,c("Phylum")], sep="_"))

round(rev(sort(rowSums(otu_table(data_summary_allPops))/sum(rowSums(otu_table(data_summary_allPops)))))*100, 3)

######################
# Data normalisation #
######################

# The library sizes between samples and groups is highly variable, and therefore comparing the data to each other will result in biased results. 

# There are two ways of handling this: 
# 1. Performing a transformation of the data or
# 2. rarefying the data.

## Centered log-ration transformation --------

pop_comparison <- AllREadsSE_Indo_Counts_physeq %>%
  tax_glom("Phylum")

# change taxa names to reflect Phylum level
taxa_names(pop_comparison) <- paste(tax_table(pop_comparison)[,"Superkingdom"], tax_table(pop_comparison)[,"Kingdom"], tax_table(pop_comparison)[,"Phylum"], sep="_")

zComposition_estimate <- otu_table(pop_comparison) + 1
zComposition_clr <- microbiome::transform(zComposition_estimate, "clr")

# add in zCompositions information to new phyloseq object
merged_phylo_counts_zComposition <- pop_comparison
taxa_zComposition <- otu_table(zComposition_clr, taxa_are_rows = TRUE)
otu_table(merged_phylo_counts_zComposition) <- taxa_zComposition

# save the AllREadsSE_Indo_Counts_physeq_clr object
save(merged_phylo_counts_zComposition, file=paste0(outputdir,"AllREadsSE_Indo_Counts_physeq_clr.Rda"))


### Sample grouping ----------

# Now that we've transformed our data, we can make a PCA plot to see how each sample clusters. 
# !Note: When Euclidean distances are used in PCoA plots, it is [equivalent to a PCA plot](http://ordination.okstate.edu/overview.htm). 

# Make an ordination plot using euclidean distances
zComposition.ord <- ordinate(merged_phylo_counts_zComposition, method = "PCoA", distance = "euclidean")

# plot by Island
pdf(paste0(outputdir,"PCA_Dimension1to5_CLRTransform_ByIsland.pdf"), width = 8)
plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="SamplePop", axes = 1:2) + geom_point(aes(), alpha=0.6, size=6) + theme(plot.title = element_text(hjust = 0, size = 16), axis.text=element_text(size=16), axis.title=element_text(size=17), legend.text = element_text(size=14), legend.title = element_text(size=16)) + scale_colour_manual(values=c(KorowaiCol,MentawaiCol,SumbaCol))
plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="SamplePop", axes = 2:3, label="SampleName") + scale_colour_manual(values=c(KorowaiCol,MentawaiCol,SumbaCol)) + geom_point(aes(), alpha=0.6, size=4) + theme(plot.title = element_text(hjust = 0, size = 12))
plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="SamplePop", axes = 3:4, label="SampleName") + scale_colour_manual(values=c(KorowaiCol,MentawaiCol,SumbaCol)) + geom_point(aes(), alpha=0.6, size=4) + theme(plot.title = element_text(hjust = 0, size = 12))
plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="SamplePop", axes = 4:5, label="SampleName") + scale_colour_manual(values=c(KorowaiCol,MentawaiCol,SumbaCol)) + geom_point(aes(), alpha=0.6, size=4) + theme(plot.title = element_text(hjust = 0, size = 12))
dev.off()

# plot by plasmodium load
logged_phyla_counts = log10(colSums(otu_table(AllREadsSE_Indo_Counts_physeq)[grep("Apicomplexa",tax_table(AllREadsSE_Indo_Counts_physeq)[,"Phylum"])])+1)
sample_data(merged_phylo_counts_zComposition)[["Apicomplexa"]] = logged_phyla_counts
zComposition.ord <- ordinate(merged_phylo_counts_zComposition, method = "PCoA", distance = "euclidean")
pdf(paste0(outputdir,"PCA_Dimension1to5_CLRTransform_ByPlasmodiumLoad.pdf"), width = 8)
plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="Apicomplexa", axes = 1:2, shape = "SamplePop") + geom_point(aes(), alpha=0.6, size=6) + theme(plot.title = element_text(hjust = 0, size = 16), axis.text=element_text(size=16), axis.title=element_text(size=17), legend.text = element_text(size=14), legend.title = element_text(size=16)) + scale_colour_gradientn(colours = colorRampPalette(c("#FFCCBB","#800026"))(10), limits=c(1, 6))
plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="Apicomplexa", axes = 2:3, label="SampleName") + geom_point(aes(), alpha=0.6, size=4) + theme(plot.title = element_text(hjust = 0, size = 12)) + scale_colour_gradientn(colours = colorRampPalette(c("#FFCCBB","#800026"))(10), limits=c(1, 6)) + theme_bw(base_size = 14)
plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="Apicomplexa", axes = 3:4, label="SampleName") + geom_point(aes(), alpha=0.6, size=4) + theme(plot.title = element_text(hjust = 0, size = 12)) + scale_colour_gradientn(colours = colorRampPalette(c("#FFCCBB","#800026"))(10), limits=c(1, 6))
plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="Apicomplexa", axes = 4:5, label="SampleName") + geom_point(aes(), alpha=0.6, size=4) + theme(plot.title = element_text(hjust = 0, size = 12)) + scale_colour_gradientn(colours = colorRampPalette(c("#FFCCBB","#800026"))(10), limits=c(1, 6))
dev.off()

# plot by Flavivirus load
logged_phyla_counts_flavivirus = log10(colSums(otu_table(AllREadsSE_Indo_Counts_physeq)[grep("Kitrinoviricota",tax_table(AllREadsSE_Indo_Counts_physeq)[,"Phylum"])])+1)
sample_data(merged_phylo_counts_zComposition)[["Flavivirus"]] = logged_phyla_counts_flavivirus
zComposition.ord <- ordinate(merged_phylo_counts_zComposition, method = "PCoA", distance = "euclidean")
pdf(paste0(outputdir,"PCA_Dimension1to5_CLRTransform_ByFlavivirusLoad.pdf"), width = 8)
#plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="Flavivirus", axes = 1:2, shape = "SamplePop") + geom_point(aes(), alpha=0.6, size=6) + theme(plot.title = element_text(hjust = 0, size = 12)) + scale_colour_gradientn(colours = colorRampPalette(c("#78c679","#006837"))(20), limits=c(1, 4))
plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="Flavivirus", axes = 1:2, shape = "SamplePop") + geom_point(aes(), alpha=0.6, size=6) + theme(plot.title = element_text(hjust = 0, size = 16), axis.text=element_text(size=16), axis.title=element_text(size=17), legend.text = element_text(size=14), legend.title = element_text(size=16)) + scale_colour_gradientn(colours = colorRampPalette(c("#78c679","#006837"))(20), limits=c(1, 4))
plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="Flavivirus", axes = 2:3, label="SampleName") + geom_point(aes(), alpha=0.6, size=4) + theme(plot.title = element_text(hjust = 0, size = 12)) + scale_colour_gradientn(colours = colorRampPalette(c("#78c679","#006837"))(20), limits=c(1, 4))
plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="Flavivirus", axes = 3:4, label="SampleName") + geom_point(aes(), alpha=0.6, size=4) + theme(plot.title = element_text(hjust = 0, size = 12)) + scale_colour_gradientn(colours = colorRampPalette(c("#78c679","#006837"))(20), limits=c(1, 4))
plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="Flavivirus", axes = 4:5, label="SampleName") + geom_point(aes(), alpha=0.6, size=4) + theme(plot.title = element_text(hjust = 0, size = 12)) + scale_colour_gradientn(colours = colorRampPalette(c("#78c679","#006837"))(20), limits=c(1, 4))
dev.off()

# plot by Actinobacteria load
logged_phyla_counts_actinobacter = log10(colSums(otu_table(AllREadsSE_Indo_Counts_physeq)[grep("Actinobacteria",tax_table(AllREadsSE_Indo_Counts_physeq)[,"Phylum"])])+1)
sample_data(merged_phylo_counts_zComposition)[["Actinobacteria"]] = logged_phyla_counts_actinobacter 
zComposition.ord <- ordinate(merged_phylo_counts_zComposition, method = "PCoA", distance = "euclidean")
pdf(paste0(outputdir,"PCA_Dimension1to5_CLRTransform_ByFlavivirusLoad.pdf"), width = 8)
#plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="Flavivirus", axes = 1:2, shape = "SamplePop") + geom_point(aes(), alpha=0.6, size=6) + theme(plot.title = element_text(hjust = 0, size = 12)) + scale_colour_gradientn(colours = colorRampPalette(c("#78c679","#006837"))(20), limits=c(1, 4))
plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="Actinobacteria", axes = 3:4, shape = "SamplePop") + geom_point(aes(), alpha=0.6, size=6) + theme(plot.title = element_text(hjust = 0, size = 16), axis.text=element_text(size=16), axis.title=element_text(size=17), legend.text = element_text(size=14), legend.title = element_text(size=16)) + scale_colour_gradientn(colours = colorRampPalette(c("lightblue","blue"))(5), limits=c(0.5, 3))
plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="Flavivirus", axes = 2:3, label="SampleName") + geom_point(aes(), alpha=0.6, size=4) + theme(plot.title = element_text(hjust = 0, size = 12)) + scale_colour_gradientn(colours = colorRampPalette(c("#78c679","#006837"))(20), limits=c(1, 4))
plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="Flavivirus", axes = 3:4, label="SampleName") + geom_point(aes(), alpha=0.6, size=4) + theme(plot.title = element_text(hjust = 0, size = 12)) + scale_colour_gradientn(colours = colorRampPalette(c("#78c679","#006837"))(20), limits=c(1, 4))
plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="Flavivirus", axes = 4:5, label="SampleName") + geom_point(aes(), alpha=0.6, size=4) + theme(plot.title = element_text(hjust = 0, size = 12)) + scale_colour_gradientn(colours = colorRampPalette(c("#78c679","#006837"))(20), limits=c(1, 4))
dev.off()

# plot by batch
pdf(paste0(outputdir,"PCA_Dimension1to5_CLRTransform_ByBatch.pdf"))
plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="batch", axes = 1:2, label="SampleName") + scale_colour_manual(values=c("#cb54d6","#3d45c4","black")) + geom_point(aes(), alpha=0.6, size=4) + theme(plot.title = element_text(hjust = 0, size = 12))
plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="batch", axes = 2:3, label="SampleName") + scale_colour_manual(values=c("#cb54d6","#3d45c4","black")) + geom_point(aes(), alpha=0.6, size=4) + theme(plot.title = element_text(hjust = 0, size = 12))
plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="batch", axes = 3:4, label="SampleName") + scale_colour_manual(values=c("#cb54d6","#3d45c4","black")) + geom_point(aes(), alpha=0.6, size=4) + theme(plot.title = element_text(hjust = 0, size = 12))
plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="batch", axes = 4:5, label="SampleName") + scale_colour_manual(values=c("#cb54d6","#3d45c4","black")) + geom_point(aes(), alpha=0.6, size=4) + theme(plot.title = element_text(hjust = 0, size = 12))
dev.off()

# Perform ANOVA and get matrix of significance levels for all taxa
all_taxa <- sapply(strsplit(rownames(tax_table(merged_phylo_counts_zComposition)), "[_.]"), `[`, 4)
all_taxa <- all_taxa %>% gsub("unk", NA, .) %>% gsub("\\bp\\b", NA, .)
all_taxa <- all_taxa[-which(is.na(all_taxa))]
zComposition.ord <- ordinate(merged_phylo_counts_zComposition, method = "PCoA", distance = "euclidean")
taxa_matrix <- matrix(ncol=5, nrow=length(all_taxa))
rownames(taxa_matrix) <- all_taxa
for (i in 1:5){
  for (taxa in all_taxa){
   logged_phyla_counts <- log10(colSums(otu_table(AllREadsSE_Indo_Counts_physeq)[grep(taxa,tax_table(AllREadsSE_Indo_Counts_physeq)[,"Phylum"])])+1)
   taxa_matrix[taxa,i] <- anova(lm(zComposition.ord$vectors[,i] ~ logged_phyla_counts))[1,5]
 }
}
write.table(taxa_matrix, file=paste0(outputdir,"ANOVA_taxa_matrix.txt"), sep="\t")

# Finally, let's see how sample grouping looks like by another clustering method - hierarchical clustering by Euclidean distance. Again, this is done on the CLR-transformed data to correct for sequencing depth.
ps_otu <- data.frame(phyloseq::otu_table(merged_phylo_counts_zComposition))
ps_otu <- t(ps_otu)
bc_dist <- vegan::vegdist(ps_otu, method = "euclidean")
ward <- as.dendrogram(hclust(bc_dist, method = "ward.D2"))
#Provide color codes
meta <- data.frame(phyloseq::sample_data(merged_phylo_counts_zComposition))
colorCode <- c(MPI = KorowaiCol, `MTW` = MentawaiCol, `SMB` = SumbaCol)
labels_colors(ward) <- colorCode[meta$SamplePop][order.dendrogram(ward)]
#Plot
pdf(paste0(outputdir,"hierarchicalClustering_ByIsland.pdf"))
plot(ward, main="Island")
dev.off()

# Label colours by Plasomodium load
initial = .bincode(meta$Apicomplexa, breaks=seq(min(meta$Apicomplexa, na.rm=T), max(meta$Apicomplexa, na.rm=T), len = 80),include.lowest = TRUE)
plasmoCol <- colorRampPalette(c("darkblue", "red", "firebrick"))(80)[initial]
labels_colors(ward)=plasmoCol[order.dendrogram(ward)]
pdf(paste0(outputdir,"hierarchicalClustering_ByPlasmodiumLoad.pdf"))
plot(ward, main="Plasmodium")
dev.off()

# Label colours by Flavivirus load
ps_otu <- data.frame(phyloseq::otu_table(merged_phylo_counts_zComposition))
colnames(ps_otu) <-  gsub("MPI", "KOR", colnames(ps_otu))
ps_otu <- t(ps_otu)
bc_dist <- vegan::vegdist(ps_otu, method = "euclidean")
ward <- as.dendrogram(hclust(bc_dist, method = "ward.D2"))
#Provide color codes
rownames(meta) <- gsub("MPI", "KOR", rownames(meta))
initial = .bincode(meta$Flavivirus, breaks=seq(min(meta$Flavivirus, na.rm=T), max(meta$Flavivirus, na.rm=T), len = 80),include.lowest = TRUE)
flavoCol <- colorRampPalette(c("darkblue", "#78c679","#006837"))(79)[initial]
labels_colors(ward)=flavoCol[order.dendrogram(ward)]
pdf(paste0(outputdir,"hierarchicalClustering_ByFlavivirusLoad.pdf"))
plot(ward, main="Flavivirus")
dev.off()

pdf(paste0(outputdir,"hierarchicalClustering_ByAllPathogens.pdf"), width = 15)
fviz_dend(ward, k = 3,                 # Cut in four groups
          cex = 0.9,                 # label size
          k_colors = c("#78c679", "#b20000","#616161"), 
          color_labels_by_k = TRUE,  # color labels by groups
          ggtheme = theme_classic(base_size=24),     # Change theme
          main = ""
          )
dev.off()

# In the PCA and hierarchical clustering, there is no obvious clustering by island, however some samples do cluster by plasmodium load.

##############################
# Relative frequency of taxa #
##############################

# One of the questions we're most interested in when investigating these samples is: what is in the data? One of the ways to do this is by visualising the data itself. A common way to do this is by looking at a stacked barplot, one for each sample, composed of the relative frequency of taxa in that sample.
# Why do we look at relative frequency? As we've seen before, our library sizes are uneven, and therefore we want to see what the proportion of each taxa is in each sample. However, compositional data does not account for the fact that as one species goes up, it will force another species to go down (i.e., it is bounded).
# I'll try to solve this problem in a few ways. The first way is by plotting the relative taxa, then plotting the CLR-transformed data (which is unbounded), and finally, seeing how this compares to rarefied data.
# To visualise relative frequency, I want to show my taxa at the family level, however since there are over 200 unique taxa at the family level, this would be difficult to visualise with so many colours. Instead, we'll make a new taxa variable which combines Superkingdom information with Family-level information. This way, we can highlight colours by superkingdom (which isn't so visually overwhelming), and still preserve Family-level information.

## Unsubsampled compositional data -----------

# add a new column containing family names and superkingdom
tax_table(AllREadsSE_Indo_Counts_physeq)[,"Superkingdom"] = paste(tax_table(AllREadsSE_Indo_Counts_physeq)[,"Superkingdom"], tax_table(AllREadsSE_Indo_Counts_physeq)[,"Family"], sep="_")
tax_table(AllREadsSE_Indo_Counts_physeq)[,"Superkingdom"] <- gsub("Bacteria_$", "Bacteria_unclassified", tax_table(AllREadsSE_Indo_Counts_physeq)[,"Superkingdom"])
tax_table(AllREadsSE_Indo_Counts_physeq)[,"Superkingdom"] <- gsub("Eukaryota_$", "Eukaryota_unclassified", tax_table(AllREadsSE_Indo_Counts_physeq)[,"Superkingdom"])
tax_table(AllREadsSE_Indo_Counts_physeq)[,"Superkingdom"] <- gsub("Viruses_$", "Viruses_unclassified", tax_table(AllREadsSE_Indo_Counts_physeq)[,"Superkingdom"])

# As pointed out, we have a lot of taxa at the family level, and it would be hard to look over everything at once. Instead, we can focus on the most prevalent taxa and highlight everything else in another colour.
# Here, I chose to highlight the top 20 taxa, since that's still representative while not being too visually exhausting.

aggregated_phyloCounts <- aggregate_top_taxa(AllREadsSE_Indo_Counts_physeq, "Superkingdom", top = 20)
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

PaletteBacteria = colorRampPalette(c("#023858","#74a9cf"))(14)
PaletteEukaryote = colorRampPalette(c("#4c0000","#b20000","#ff4c4c"))(4)
PaletteOther = colorRampPalette(c("black"))(1)
PaletteVirus = colorRampPalette(c("#78c679","#006837"))(2)

Merged_Palette <- c(PaletteBacteria,PaletteEukaryote,PaletteOther,PaletteVirus)

pdf(paste0(outputdir,"relativeTaxa_Compositional.pdf"), width=15)
phyloseq::plot_bar(relative_phyloCounts, fill = "Superkingdom") +
  geom_bar(aes(fill = Superkingdom), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +  theme_bw(base_size = 15) +
  facet_wrap(~ SamplePop, scales = "free") + scale_fill_manual(values=Merged_Palette) +
  theme(panel.background = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
dev.off()

# Now plot with a circular barplot
data=as.matrix(as.data.frame(otu_table(relative_phyloCounts)))
data=t(data)
data=as.data.frame(data)
data$group=sapply(strsplit(rownames(data), "[-.]"), `[`, 1)
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

pdf(paste0(outputdir,"relativeTaxa_Compositional_CircularBarplot.pdf"), width=15)
p
dev.off()

##################################
# Differential abundance testing #
##################################

# We're interested in testing whether the species composition between islands is significantly different. One way to test this is through differential abundance testing.
# One of the best packages I've found so far to test this is called [Aldex2](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0067019) (ANOVA-Like Differential Expression). Aldex2 works by generatating a posterior probability (128 by default) for an observed instance of a taxon (adding a small prior to deal with zeros), then performing a centered log-ratio transformation on the data as a normalisation step (this deals with uneven library sizes). To identify differentially expressed genes, Aldex2 then performs a significance test using a Wilcoxon rank test, and finally, the probability of the taxon being differentially abundant is adjusted with FDR correction (by Benjamini–Hochberg). 
# Aldex2 corrects for uneven library sizes, so rarefying is not necessary. The only input we need is the data with singletons removed.
# Since Aldex2 works by comparing two datasets, we need to subset the datasets into groups of two. 

# Differential abundance testing
MPIVsMTW=subset_samples(pop_comparison, SamplePop != "SMB")
any(taxa_sums(MPIVsMTW) == 0)
# remove any 0s from the data
MPIVsMTW <- prune_taxa(taxa_sums(MPIVsMTW) > 0, MPIVsMTW)
taxa_names(MPIVsMTW)=make.unique(tax_table(MPIVsMTW)[,"Phylum"])
# Run aldex2
aldex2_MPIVsMTW <- ALDEx2::aldex(data.frame(phyloseq::otu_table(MPIVsMTW)), phyloseq::sample_data(MPIVsMTW)$SamplePop, test="t", effect = TRUE)

# we.eBH - Expected Benjamini-Hochberg corrected P value of Welch’s t test

sig_aldex2_MPIVsMTW <- aldex2_MPIVsMTW %>%
  rownames_to_column(var = "OTU") %>%
  filter(we.eBH < 0.05) %>%
  arrange(effect, we.eBH) %>%
  dplyr::select(OTU, diff.btw, diff.win, effect, we.ep, we.eBH)

labels = sapply(strsplit(rownames(aldex2_MPIVsMTW), "[..]"), `[`, 1) %>% gsub("unk_p","Uncharacterised",.)
taxa_superkingdom = sapply(strsplit(tax_table(MPIVsMTW)[,"Superkingdom"], "[_.]"), `[`, 1) %>% gsub("unk","Uncharacterised",.)

pdf(paste0(outputdir,"differentialAbundance_MPIVsMTW.pdf"))
ggplot(aldex2_MPIVsMTW) +
  geom_point(aes(x = effect, y = -log10(we.eBH)), color = ifelse(aldex2_MPIVsMTW$we.eBH <= 0.05, c("#023858","#800026","grey","#78c679")[as.numeric(as.factor(taxa_superkingdom))],"black"), alpha = 0.65, size=8) +
  #geom_text_repel(aes(x = effect, y = -log10(wi.eBH), label = rownames(aldex2_IndoVsDutch))) +
  geom_text_repel(aes(x = effect, y = -log10(we.eBH), label = ifelse(we.eBH <= 0.05, labels,""))) +
  ggtitle("Mentawai Versus Korowai") + xlim(c(-1,0.5)) + ylim(c(0,2.75)) +
  xlab("Effect Size") +  
  ylab("-log10 adjusted p-value") + theme_bw(base_size = 18) +
  theme(legend.position = "bottom",
        plot.title = element_text(size = rel(1.5), hjust = 0.75),
        axis.title = element_text(size = rel(1.25)))
dev.off()

# Now let's run it on the Korowai versus people from the island of Sumba.

MPIVsSMB=subset_samples(pop_comparison, SamplePop != "MTW")
any(taxa_sums(MPIVsSMB) == 0)
MPIVsSMB <- prune_taxa(taxa_sums(MPIVsSMB) > 0, MPIVsSMB)
taxa_names(MPIVsSMB)=make.unique(tax_table(MPIVsSMB)[,"Phylum"])
aldex2_MPIVsSMB <- ALDEx2::aldex(data.frame(phyloseq::otu_table(MPIVsSMB)), phyloseq::sample_data(MPIVsSMB)$SamplePop, test="t", effect = TRUE)
sig_aldex2_MPIVsSMB <- aldex2_MPIVsSMB %>%
  rownames_to_column(var = "OTU") %>%
  filter(we.eBH < 0.05) %>%
  arrange(effect, we.eBH) %>%
  dplyr::select(OTU, diff.btw, diff.win, effect, we.ep, we.eBH)

labels = sapply(strsplit(rownames(aldex2_MPIVsSMB), "[..]"), `[`, 1) %>% gsub("unk_p","Uncharacterised",.)
taxa_superkingdom = sapply(strsplit(tax_table(MPIVsSMB)[,"Superkingdom"], "[_.]"), `[`, 1) %>% gsub("unk","Uncharacterised",.)

pdf(paste0(outputdir,"differentialAbundance_MPIVsSMB.pdf"))
ggplot(aldex2_MPIVsSMB) +
  geom_point(aes(x = effect, y = -log10(we.eBH)), color = ifelse(aldex2_MPIVsSMB$we.eBH <= 0.05, c("grey","#023858","#800026","grey","#78c679")[as.numeric(as.factor(taxa_superkingdom))],"black"), alpha = 0.65, size=8) +
  #geom_text_repel(aes(x = effect, y = -log10(wi.eBH), label = rownames(aldex2_IndoVsDutch))) +
  geom_text_repel(aes(x = effect, y = -log10(we.eBH), label = ifelse(we.eBH <= 0.05, labels,""))) +
  ggtitle("Sumba Versus Korowai") + xlim(c(-1,0.5)) + ylim(c(0,2.75)) +
  xlab("Effect Size") + 
  ylab("-log10 adjusted p-value") + theme_bw(base_size = 18) +
  theme(legend.position = "bottom",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))
dev.off()

# Finally, test this on Sumba versus Mentawai

MTWVsSMB=subset_samples(pop_comparison, SamplePop != "MPI")
any(taxa_sums(MTWVsSMB) == 0)
MTWVsSMB <- prune_taxa(taxa_sums(MTWVsSMB) > 0, MTWVsSMB)
taxa_names(MTWVsSMB)=make.unique(tax_table(MTWVsSMB)[,"Phylum"])
aldex2_MTWVsSMB <- ALDEx2::aldex(data.frame(phyloseq::otu_table(MTWVsSMB)), phyloseq::sample_data(MTWVsSMB)$SamplePop, test="t", effect = TRUE)
sig_aldex2_MTWVsSMB <- aldex2_MTWVsSMB %>%
  rownames_to_column(var = "OTU") %>%
  filter(we.eBH < 0.05) %>%
  arrange(effect, we.eBH) %>%
  dplyr::select(OTU, diff.btw, diff.win, effect, we.ep, we.eBH)

labels = sapply(strsplit(rownames(aldex2_MTWVsSMB), "[..]"), `[`, 1) %>% gsub("unk_p","Uncharacterised",.)
taxa_superkingdom = sapply(strsplit(tax_table(MTWVsSMB)[,"Superkingdom"], "[_.]"), `[`, 1) %>% gsub("unk","Uncharacterised",.)

pdf(paste0(outputdir,"differentialAbundance_MTWVsSMB.pdf"))
ggplot(aldex2_MTWVsSMB) +
  geom_point(aes(x = effect, y = -log10(we.eBH)), color = ifelse(aldex2_MTWVsSMB$we.eBH <= 0.05, c("grey","#023858","#800026","grey","#78c679")[as.numeric(as.factor(taxa_superkingdom))],"black"), alpha = 0.65, size=8) +
  #geom_text_repel(aes(x = effect, y = -log10(wi.eBH), label = rownames(aldex2_IndoVsDutch))) +
  geom_text_repel(aes(x = effect, y = -log10(we.eBH), label = ifelse(we.eBH <= 0.05, labels,""))) +
  ggtitle("Sumba Versus Mentawai") + xlim(c(-1,0.5)) + ylim(c(0,2.75)) +
  xlab("Effect Size") + 
  ylab("-log10 adjusted p-value") + theme_bw(base_size = 18) +
  theme(legend.position = "bottom",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))
dev.off()

# Try this at the family level

pop_comparison_family <- AllREadsSE_Indo_Counts_physeq %>%
  tax_glom("Family")

# Differential abundance testing
MPIVsMTW=subset_samples(pop_comparison_family, SamplePop != "SMB")
any(taxa_sums(MPIVsMTW) == 0)
# remove any 0s from the data
MPIVsMTW <- prune_taxa(taxa_sums(MPIVsMTW) > 0, MPIVsMTW)
taxa_names(MPIVsMTW)=make.unique(tax_table(MPIVsMTW)[,"Family"])
# Run aldex2
aldex2_MPIVsMTW <- ALDEx2::aldex(data.frame(phyloseq::otu_table(MPIVsMTW)), phyloseq::sample_data(MPIVsMTW)$SamplePop, test="t", effect = TRUE)

sig_aldex2_MPIVsMTW <- aldex2_MPIVsMTW %>%
  rownames_to_column(var = "OTU") %>%
  filter(we.eBH < 0.05) %>%
  arrange(effect, we.eBH) %>%
  dplyr::select(OTU, diff.btw, diff.win, effect, we.ep, we.eBH)


# Now let's run it on the Korowai versus people from the island of Sumba.

MPIVsSMB=subset_samples(pop_comparison_family, SamplePop != "MTW")
any(taxa_sums(MPIVsSMB) == 0)
MPIVsSMB <- prune_taxa(taxa_sums(MPIVsSMB) > 0, MPIVsSMB)
taxa_names(MPIVsSMB)=make.unique(tax_table(MPIVsSMB)[,"Family"])
aldex2_MPIVsSMB <- ALDEx2::aldex(data.frame(phyloseq::otu_table(MPIVsSMB)), phyloseq::sample_data(MPIVsSMB)$SamplePop, test="t", effect = TRUE)
sig_aldex2_MPIVsSMB <- aldex2_MPIVsSMB %>%
  rownames_to_column(var = "OTU") %>%
  filter(we.eBH < 0.05) %>%
  arrange(effect, we.eBH) %>%
  dplyr::select(OTU, diff.btw, diff.win, effect, we.ep, we.eBH)

# Finally, let's test if there are differentially abundant taxa in Mentawai versus people from the island of Sumba.

MTWVsSMB=subset_samples(pop_comparison_family, SamplePop != "MPI")
any(taxa_sums(MTWVsSMB) == 0)
MTWVsSMB <- prune_taxa(taxa_sums(MTWVsSMB) > 0, MTWVsSMB)
taxa_names(MTWVsSMB)=make.unique(tax_table(MTWVsSMB)[,"Family"])
aldex2_MTWVsSMB <- ALDEx2::aldex(data.frame(phyloseq::otu_table(MTWVsSMB)), phyloseq::sample_data(MTWVsSMB)$SamplePop, test="t", effect = TRUE)
sig_aldex2_MTWVsSMB <- aldex2_MTWVsSMB %>%
  rownames_to_column(var = "OTU") %>%
  filter(we.eBH < 0.05) %>%
  arrange(effect, we.eBH) %>%
  dplyr::select(OTU, diff.btw, diff.win, effect, we.ep, we.eBH)

# there are no statistically significant differentially abundant taxa.

###################
# Alpha diversity #
###################

# DivNet is a package that adresses these problems. Divnet is a method for estimating within- and between-community diversity in ecosystems where taxa interact via an ecological network. It accounts for differences in sequencing depth and estimates the number of missing species based on the sequence depth and
# number of rare taxa in the data.
# To use DivNet, you need unsubsampled data without removing singletons.

# Remove Viridiplantae and Metazoa
AllREadsSE_Indo_Counts_physeq_withSingletons <- subset_taxa(AllREadsSE_Indo_Counts_physeq_withSingletons, (Kingdom!="Viridiplantae"))
AllREadsSE_Indo_Counts_physeq_withSingletons <- subset_taxa(AllREadsSE_Indo_Counts_physeq_withSingletons, (Kingdom!="Metazoa"))
AllREadsSE_Indo_Counts_physeq_withSingletons <- subset_taxa(AllREadsSE_Indo_Counts_physeq_withSingletons, (Superkingdom!="unk_sk"))
# remove any empty rows
AllREadsSE_Indo_Counts_physeq_withSingletons <- prune_taxa(taxa_sums(AllREadsSE_Indo_Counts_physeq_withSingletons) > 0, AllREadsSE_Indo_Counts_physeq_withSingletons)

# Now that we have data without singletons, we now need to merge our data at a specified taxonomic level. DivNet is computationally expensive, and therefore a higher level is much, much faster.
# We'll therefore test how our groups look like at the Phylum level. Then, we'll run DivNet without specifying any hypothesis testing.

# comparing diversity at the phylum level
pop_comparison <- AllREadsSE_Indo_Counts_physeq_withSingletons %>%
  tax_glom("Phylum")

# If we don't change the sample names here from hyphens to periods, we'll get an error later
sample_names(pop_comparison) <- gsub("\\-", ".", sample_names(pop_comparison))

# Run divnet without specifying any hypothesis testing
dv_pop_comparison <- divnet(pop_comparison, ncores = 4)

# DivNet will output an object with estimates for multiple different alpha (and beta) diversity measures (we'll get to the beta diversity estimates later).
# The Shannon and Simpson index are two popular alpha diversity indices to measure species richness. For Shannon diversity, the importance of rare taxa are downweighted, since they do not play a large role in the community or they could potentially be due to error. For this reason, the Shannon index is one of the most popular alpha diversity metrics.
# To interpret the index, a higher Shannon index means higher diversity, whereas a lowed index number means lower diversity.
# Let's take out the Shannon diversity metric from DivNet and plot it.

# Now let's plot the results of shannon and Simpson diversity
summary_df_shannon <- as.data.frame(dv_pop_comparison$shannon %>%
  summary %>%
  add_column("SampleNames" = pop_comparison %>% otu_table %>% sample_names) %>%
  add_column("SamplePop" =  pop_comparison %>% sample_data %>% .[,"SamplePop"] %>% as.matrix(.) %>% .[,1] %>% unname(.)))

pdf(paste0(outputdir,"shannonDiversity_DivNet_noHypothesisTesting.pdf"))
ggplot(summary_df_shannon, aes(y = estimate, x = SamplePop, fill = SamplePop)) + geom_violin(alpha=0.7) + 
  geom_jitter(height = 0, width = .2) + geom_boxplot(width=0.08, outlier.color = NA) +
  scale_fill_manual(values=c(KorowaiCol,MentawaiCol,SumbaCol)) + ggtitle("Shannon Diversity") +
  ylab("Estimate of Shannon Diversity")
dev.off()

# We can see that, on average, the Sumba and Mentawai populations have a higher Shannon diversity while the Korowai are slightly lower. 
# We can also plot each individual sample, along with their standard deviation (another cool, and imprortant feature that DivNet calculates and uses in their hypothesis testing). 

pdf(paste0(outputdir,"shannonDiversity_DivNet_noHypothesisTesting_bySample.pdf"))
plot(dv_pop_comparison$shannon, pop_comparison, col = "SamplePop") + scale_colour_manual(values=c(KorowaiCol,MentawaiCol,SumbaCol))
dev.off()

# Now let's see how the population diversity looks like when we use the Simpson diversity index. The Simpson diversity index is a similarity index where the higher the value, the lower the diversity. It measures the probability that two individuals randomly selected from a sample will belong to the same species. With this index, 0 represents infinite diversity and 1, no diversity.

summary_df_simpson <- as.data.frame(dv_pop_comparison$simpson %>%
  summary %>%
  add_column("SampleNames" = pop_comparison %>% otu_table %>% sample_names) %>%
  add_column("SamplePop" =  pop_comparison %>% sample_data %>% .[,"SamplePop"] %>% as.matrix(.) %>% .[,1] %>% unname(.)))

pdf(paste0(outputdir,"simpsonDiversity_DivNet_noHypothesisTesting.pdf"))
ggplot(summary_df_simpson, aes(y = estimate, x = SamplePop, fill = SamplePop)) + geom_violin(alpha=0.7) + 
  geom_jitter(height = 0, width = .2) + geom_boxplot(width=0.08, outlier.color = NA) +
  scale_fill_manual(values=c(KorowaiCol,MentawaiCol,SumbaCol)) + ggtitle("Simpson's Diversity Index") +
  ylab("Estimate of Simpson Diversity")
dev.off()

# Again, we can see that the same trend holds true - Mentawai and Sumba have the highest diversity (remember, and index of 0 equates to infinite diversity), while the Korowai population has the lowest diversity.
# This is how the diversity looks like with SE included for each sample.

pdf(paste0(outputdir,"simpsonDiversity_DivNet_noHypothesisTesting_bySample.pdf"))
plot(dv_pop_comparison$simpson, pop_comparison, col = "SamplePop") + scale_colour_manual(values=c(KorowaiCol,MentawaiCol,SumbaCol))
dev.off()

# Since a larger Simpson index value equates to a lower diversity index, many people find this confusing and not
# very intuitive. Therefore, the inverse Simpsone Index, or 1 - Simpson Index, is also commonly used.
# Let's plot that now. 

# Subtract the Simpson estimate from one
summary_df_simpson$estimate = 1-summary_df_simpson$estimate
# Plot

pdf(paste0(outputdir,"InverseSimpsonDiversity_DivNet_noHypothesisTesting.pdf"))
ggplot(summary_df_simpson, aes(y = estimate, x = SamplePop, fill = SamplePop)) + geom_violin(alpha=0.7) + 
  geom_jitter(height = 0, width = .2) + geom_boxplot(width=0.08, outlier.color = NA) +
  scale_fill_manual(values=c(KorowaiCol,MentawaiCol,SumbaCol)) + ggtitle("Simpson's Diversity Index") +
  ylab("Estimate of Simpson Diversity")
dev.off()

# Now let's test the hypothesis that the diversity is different between islands. We are now estimating the diversity of island/population being an ecosystem, so we're focusing on the ecosystem, not just the samples. 
# Let's test this first using the Shannon diversity index. 

# test the hypothesis that the diversity is differnet between islands
dv_pop_comparison_cov <- pop_comparison %>%
  divnet(X = "SamplePop", ncores = 8)

pdf(paste0(outputdir,"shannonDiversity_DivNet_withHypothesisTesting.pdf"))
# Plot the results for each individual
plot(dv_pop_comparison_cov$shannon, pop_comparison, col = "SamplePop") + scale_colour_manual(values=c(KorowaiCol,MentawaiCol,SumbaCol))
dev.off()

# Interestingly, we now see that, as a population, the Korowai samples have the highest Shannon diversity, followed by the Mentawai and Sumba samples. However, you can also see that there are large error bars around the Korowai samples.
# Let's test the hypothesis that the islands are different formally.

# test that these populations are actually different
diversityTest_shannon=testDiversity(dv_pop_comparison_cov, "shannon")
write.table(diversityTest_shannon, file=paste0(outputdir,"diversityTest_Shannon.txt"))

# The result tells us that the mean Shannon diversity in the Korowai at the Phylum-level is 1.21, and it is significantly lower by 0.15 orders, on average, in the Mentawai population. We can also see that the Sumba population is significantly lower than the Korowai population by 0.15 orders, on average.
# Let's do the same thing for Simpson diversity.

pdf(paste0(outputdir,"simpsonDiversity_DivNet_withHypothesisTesting.pdf"))
plot(dv_pop_comparison_cov$simpson, pop_comparison, col = "SamplePop") + scale_colour_manual(values=c(KorowaiCol,MentawaiCol,SumbaCol))
dev.off()

# With the Simpson index, we see that the Mentawai and Sumba populations have the highest Simpson index (i.e., lowest diversity).
# Again, let's test this formally.

diversityTest_simpson=testDiversity(dv_pop_comparison_cov, "simpson")
write.table(diversityTest_simpson, file=paste0(outputdir,"diversityTest_Simpson.txt"))

# The result tells us that the mean Simpson diversity index in the korowai population at the Phylum-level is 0.5, and it is significantly higher by 0.05 orders, on average, in the Mentawai population. We can also see that the Sumba population is significantly higher than the Korowai population by 0.05 orders, on average.

##################
# Beta Diversity #
##################

# let's explore the results from DivNet to see how the Bray-Curtis dissimilarity looks like.

# First, let's look at Bray-curtis dissimilarity at the individual sample level
bray_est <- simplifyBeta(dv_pop_comparison, pop_comparison, "bray-curtis", "SamplePop")

# add in group comparisons and plot
bray_est$group=paste(bray_est$Covar1,bray_est$Covar2,sep="_")
pdf(paste0(outputdir,"BrayCurtis_DivNet_noHypothesisTesting.pdf"), width = 15)
ggplot(bray_est, aes(x = interaction(Covar1, Covar2), y = beta_est, fill=group)) +
  geom_violin(alpha=0.7) + geom_boxplot(width=0.1) + xlab("Population Comparisons") + theme_bw(base_size = 18) +
  theme(legend.position="none") + ggtitle("Bray-Curtis Distance Estimate") +
  ylab("Bray-Curtis Distance") + scale_fill_manual(values = c(KorowaiCol,"#86509c",MentawaiCol,"#f97702","#019680",SumbaCol))
dev.off()

# First, let's look at Bray-curtis dissimilarity at the individual sample level
bray_est_eucl <- simplifyBeta(dv_pop_comparison, pop_comparison, "euclidean", "SamplePop")

# add in group comparisons and plot
bray_est_eucl$group <- paste(bray_est_eucl$Covar1, bray_est_eucl$Covar2,sep="_")

pdf(paste0(outputdir,"Euclidean_DivNet_noHypothesisTesting.pdf"))
ggplot(bray_est_eucl, aes(x = interaction(Covar1, Covar2), y = beta_est, fill=group)) +
  geom_violin(alpha=0.7) + geom_boxplot(width=0.1) + xlab("Population Comparisons") + 
  theme(legend.position="none") + ggtitle("Euclidean Distance Estimate") +
  ylab("Euclidean Distance")
dev.off()

# Now let's see how this looks like for island-level comparisons.

# Bray-Curtis dissimilarity
pdf(paste0(outputdir,"BrayCurtis_DivNet_withHypothesisTesting.pdf"))
simplifyBeta(dv_pop_comparison_cov, pop_comparison, "bray-curtis", "SamplePop") %>%
  ggplot(aes(x = interaction(Covar1, Covar2), 
             y = beta_est)) +
  geom_point() +
  geom_linerange(aes(ymin = lower, ymax = upper)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") + ylab("Estimates of Bray-Curtis distance")
dev.off()

# Euclidean distance
pdf(paste0(outputdir,"Euclidean_DivNet_withHypothesisTesting.pdf"))
simplifyBeta(dv_pop_comparison_cov, pop_comparison, "euclidean", "SamplePop") %>% filter(beta_var != 0) %>% 
  ggplot(aes(x = interaction(Covar1, Covar2), 
             y = beta_est)) +
  geom_point() +
  geom_linerange(aes(ymin = lower, ymax = upper)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") + ylab("Estimates of Bray-Curtis distance")
dev.off()
