require(ggplot2)
require(RColorBrewer)
library(dplyr)
library(plyr)
library(reshape2)
library(ggpubr)
library(tidyverse)             
library(phyloseq)                   

# set up directories
outputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/ControlSampleComparison/"
refdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/EpiStudy/"

###################
# Indonesian data #
###################

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

AllREadsSE_Indo_Counts_physeq

####################
## Malian Samples ##
####################

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

###################
## Dutch samples ##
###################

control_100K_Counts <- read.csv(paste0(refdir,"Controls_NoFiltering_Counts.csv"),check.names=FALSE)

# Separate species' abundances and taxonomy columns
taxa_raw <- as.matrix(control_100K_Counts[,c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species")])
abund_raw <- as.matrix(control_100K_Counts[,-which(colnames(control_100K_Counts) %in% c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species"))])
# convert to Phyloseq object
tax = tax_table(taxa_raw)
taxa = otu_table(abund_raw, taxa_are_rows = TRUE)
control_100K_Counts_physeq = phyloseq(taxa, tax)

# get information on phyloseq object
control_100K_Counts_physeq

# text file made from script 'GetDiseaseStatus_ControlSamples.sh'
diseaseStatus=read.table(paste0(refdir,"ControlSamples_DiseaseStatus.txt"))
colnames(diseaseStatus)=c("Samples","diseaseStatus")
controlSamples=as.character(diseaseStatus[grep("CSF", diseaseStatus$diseaseStatus),"Samples"])

# prune out samples we don't want
control_100K_Counts_physeq=prune_samples(controlSamples,control_100K_Counts_physeq)

# add in sample information, i.e., the sample names and population they're from
samplenames <- colnames(otu_table(control_100K_Counts_physeq))
pop <- rep("Netherlands",ncol(otu_table(control_100K_Counts_physeq)))

# make this into a df and add to the Phloseq object
samples_df=data.frame(SampleName=colnames(otu_table(control_100K_Counts_physeq)), SamplePop=pop)
samples = sample_data(samples_df)
rownames(samples)=samples$SampleName
sample_data(control_100K_Counts_physeq) <- samples

# remove taxa with only 0's in the phyloseq object
any(taxa_sums(control_100K_Counts_physeq) == 0)
control_100K_Counts_physeq=prune_taxa(taxa_sums(control_100K_Counts_physeq) > 0, control_100K_Counts_physeq)

# get phyloseq summary information
control_100K_Counts_physeq

#################
## USA samples ##
#################

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

# Summarise Data Before Filtering

nOTUs_prefilter=c(nrow(tax_table(AllREadsSE_Indo_Counts_physeq)), nrow(tax_table(Mali_Counts_physeq)), nrow(tax_table(control_100K_Counts_physeq)), nrow(tax_table(Americans_Counts_physeq)))
nReads_prefilter=c(sum(colSums(otu_table(AllREadsSE_Indo_Counts_physeq))), sum(colSums(otu_table(Mali_Counts_physeq))), sum(colSums(otu_table(control_100K_Counts_physeq))), sum(colSums(otu_table(Americans_Counts_physeq))))

#########################
# Filter out Singletons #
#########################

AllREadsSE_Indo_Counts_physeq_noSingletons = AllREadsSE_Indo_Counts_physeq
Mali_Counts_physeq_noSingletons = Mali_Counts_physeq
control_100K_Counts_physeq_noSingletons = control_100K_Counts_physeq
Americans_Counts_physeq_noSingletons = Americans_Counts_physeq

otu_table(AllREadsSE_Indo_Counts_physeq_noSingletons)[otu_table(AllREadsSE_Indo_Counts_physeq_noSingletons)<=1]<-0
AllREadsSE_Indo_Counts_physeq_noSingletons <- prune_taxa(taxa_sums(AllREadsSE_Indo_Counts_physeq_noSingletons) > 0, AllREadsSE_Indo_Counts_physeq_noSingletons)
otu_table(Mali_Counts_physeq_noSingletons)[otu_table(Mali_Counts_physeq_noSingletons)<=1]<-0
Mali_Counts_physeq_noSingletons <- prune_taxa(taxa_sums(Mali_Counts_physeq_noSingletons) > 0, Mali_Counts_physeq_noSingletons)
otu_table(control_100K_Counts_physeq_noSingletons)[otu_table(control_100K_Counts_physeq_noSingletons)<=1]<-0
control_100K_Counts_physeq_noSingletons <- prune_taxa(taxa_sums(control_100K_Counts_physeq_noSingletons) > 0, control_100K_Counts_physeq_noSingletons)
otu_table(Americans_Counts_physeq_noSingletons)[otu_table(Americans_Counts_physeq_noSingletons)<=1]<-0
Americans_Counts_physeq_noSingletons <- prune_taxa(taxa_sums(Americans_Counts_physeq_noSingletons) > 0, Americans_Counts_physeq_noSingletons)

# Summarise data

nOTUs_noSingletons=c(nrow(tax_table(AllREadsSE_Indo_Counts_physeq_noSingletons)), nrow(tax_table(Mali_Counts_physeq_noSingletons)), nrow(tax_table(control_100K_Counts_physeq_noSingletons)), nrow(tax_table(Americans_Counts_physeq_noSingletons)))
nReads_noSingletons=c(sum(colSums(otu_table(AllREadsSE_Indo_Counts_physeq_noSingletons))), sum(colSums(otu_table(Mali_Counts_physeq_noSingletons))), sum(colSums(otu_table(control_100K_Counts_physeq_noSingletons))), sum(colSums(otu_table(Americans_Counts_physeq_noSingletons))))

############################
# Filter out Viridiplantae #
############################

AllREadsSE_Indo_Counts_physeq_noPlants <- subset_taxa(AllREadsSE_Indo_Counts_physeq_noSingletons, (Kingdom!="Viridiplantae"))
AllREadsSE_Indo_Counts_physeq_noPlants <- prune_taxa(taxa_sums(AllREadsSE_Indo_Counts_physeq_noPlants) > 0, AllREadsSE_Indo_Counts_physeq_noPlants)
Mali_Counts_physeq_noPlants <- subset_taxa(Mali_Counts_physeq_noSingletons, (Kingdom!="Viridiplantae"))
Mali_Counts_physeq_noPlants <- prune_taxa(taxa_sums(Mali_Counts_physeq_noPlants) > 0, Mali_Counts_physeq_noPlants)
control_100K_Counts_physeq_noPlants <- subset_taxa(control_100K_Counts_physeq_noSingletons, (Kingdom!="Viridiplantae"))
control_100K_Counts_physeq_noPlants <- prune_taxa(taxa_sums(control_100K_Counts_physeq_noPlants) > 0, control_100K_Counts_physeq_noPlants)
Americans_Counts_physeq_noPlants <- subset_taxa(Americans_Counts_physeq_noSingletons, (Kingdom!="Viridiplantae"))
Americans_Counts_physeq_noPlants <- prune_taxa(taxa_sums(Americans_Counts_physeq_noPlants) > 0, Americans_Counts_physeq_noPlants)

# Summarise data

nOTUs_noPlants=c(nrow(tax_table(AllREadsSE_Indo_Counts_physeq_noPlants)), nrow(tax_table(Mali_Counts_physeq_noPlants)), nrow(tax_table(control_100K_Counts_physeq_noPlants)), nrow(tax_table(Americans_Counts_physeq_noPlants)))
nReads_noPlants=c(sum(colSums(otu_table(AllREadsSE_Indo_Counts_physeq_noPlants))), sum(colSums(otu_table(Mali_Counts_physeq_noPlants))), sum(colSums(otu_table(control_100K_Counts_physeq_noPlants))), sum(colSums(otu_table(Americans_Counts_physeq_noPlants))))

#######################
# Filter out Metazoa #
#######################

AllREadsSE_Indo_Counts_physeq_noAnimals <- subset_taxa(AllREadsSE_Indo_Counts_physeq_noPlants, (Kingdom!="Metazoa"))
AllREadsSE_Indo_Counts_physeq_noAnimals <- prune_taxa(taxa_sums(AllREadsSE_Indo_Counts_physeq_noAnimals) > 0, AllREadsSE_Indo_Counts_physeq_noAnimals)
Mali_Counts_physeq_noAnimals <- subset_taxa(Mali_Counts_physeq_noPlants, (Kingdom!="Metazoa"))
Mali_Counts_physeq_noAnimals <- prune_taxa(taxa_sums(Mali_Counts_physeq_noAnimals) > 0, Mali_Counts_physeq_noAnimals)
control_100K_Counts_physeq_noAnimals <- subset_taxa(control_100K_Counts_physeq_noPlants, (Kingdom!="Metazoa"))
control_100K_Counts_physeq_noAnimals <- prune_taxa(taxa_sums(control_100K_Counts_physeq_noAnimals) > 0, control_100K_Counts_physeq_noAnimals)
Americans_Counts_physeq_noAnimals <- subset_taxa(Americans_Counts_physeq_noPlants, (Kingdom!="Metazoa"))
Americans_Counts_physeq_noAnimals <- prune_taxa(taxa_sums(Americans_Counts_physeq_noAnimals) > 0, Americans_Counts_physeq_noAnimals)

# Summarise data

nOTUs_noMetazoa=c(nrow(tax_table(AllREadsSE_Indo_Counts_physeq_noAnimals)), nrow(tax_table(Mali_Counts_physeq_noAnimals)), nrow(tax_table(control_100K_Counts_physeq_noAnimals)), nrow(tax_table(Americans_Counts_physeq_noAnimals)))
nReads_noMetazoa=c(sum(colSums(otu_table(AllREadsSE_Indo_Counts_physeq_noAnimals))), sum(colSums(otu_table(Mali_Counts_physeq_noAnimals))), sum(colSums(otu_table(control_100K_Counts_physeq_noAnimals))), sum(colSums(otu_table(Americans_Counts_physeq_noAnimals))))

##################
# Visualisation #
##################

# Make table and visualise
filtering_summary=data.frame(OTUs_Prefiltering=nOTUs_prefilter,Reads_Prefiltering=nReads_prefilter,OTUs_noSingletons=nOTUs_noSingletons, Reads_noSingletons=nReads_noSingletons, OTUs_noPlants=nOTUs_noPlants, Reads_noPlants=nReads_noPlants, OTUs_noMetazoa=nOTUs_noMetazoa, Reads_noMetazoa=nReads_noMetazoa)
rownames(filtering_summary) = c("Indonesia","Mali", "Netherlands", "USA")

# save file
write.table(filtering_summary,file=paste0(outputdir,"Filtering_Summary.txt"),sep="\t")

# Plot

# make this into a figure
filtering_summary=t(filtering_summary)
# Just interested in Reads
filtering_summary=filtering_summary[grep("Reads",rownames(filtering_summary)),]
filtering_summary=as.data.frame(filtering_summary)
stage=sapply(strsplit(rownames(filtering_summary), "[_.]"), `[`, 2)
type=sapply(strsplit(rownames(filtering_summary), "[_.]"), `[`, 1)
filtering_summary$Stage=stage
filtering_summary$Type=type
melted_filtering_summary=melt(filtering_summary)
melted_filtering_summary$Stage <- factor(melted_filtering_summary$Stage, levels = c("Prefiltering", "noSingletons", "noPlants", "noMetazoa"))

pdf(paste0(outputdir,"filtering_summary.pdf"), width=10)
ggplot(data=melted_filtering_summary, aes(x=Stage, y=log10(value), fill=variable)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired") + facet_wrap(~Type, scales = "free") + theme_bw()
dev.off()

# Add in new column, Superkimgdom + family, to all phyloseq objects and save data
AllREadsSE_Indo_Counts_physeq <- tax_glom(AllREadsSE_Indo_Counts_physeq_noAnimals, "Family")
any(taxa_sums(AllREadsSE_Indo_Counts_physeq) == 0)
taxa_names(AllREadsSE_Indo_Counts_physeq) <- paste(tax_table(AllREadsSE_Indo_Counts_physeq)[,"Superkingdom"], tax_table(AllREadsSE_Indo_Counts_physeq)[,"Kingdom"], tax_table(AllREadsSE_Indo_Counts_physeq)[,"Phylum"], tax_table(AllREadsSE_Indo_Counts_physeq)[,"Class"], tax_table(AllREadsSE_Indo_Counts_physeq)[,"Order"], tax_table(AllREadsSE_Indo_Counts_physeq)[,"Family"], sep="_")
Mali_Counts_physeq <- tax_glom(Mali_Counts_physeq_noAnimals, "Family")
any(taxa_sums(Mali_Counts_physeq) == 0)
taxa_names(Mali_Counts_physeq) <- paste(tax_table(Mali_Counts_physeq)[,"Superkingdom"], tax_table(Mali_Counts_physeq)[,"Kingdom"], tax_table(Mali_Counts_physeq)[,"Phylum"], tax_table(Mali_Counts_physeq)[,"Class"], tax_table(Mali_Counts_physeq)[,"Order"], tax_table(Mali_Counts_physeq)[,"Family"], sep="_")
control_100K_Counts_physeq <- tax_glom(control_100K_Counts_physeq_noAnimals, "Family")
any(taxa_sums(control_100K_Counts_physeq) == 0)
taxa_names(control_100K_Counts_physeq) <- paste(tax_table(control_100K_Counts_physeq)[,"Superkingdom"], tax_table(control_100K_Counts_physeq)[,"Kingdom"], tax_table(control_100K_Counts_physeq)[,"Phylum"], tax_table(control_100K_Counts_physeq)[,"Class"], tax_table(control_100K_Counts_physeq)[,"Order"], tax_table(control_100K_Counts_physeq)[,"Family"], sep="_")
Americans_Counts_physeq <- tax_glom(Americans_Counts_physeq_noAnimals, "Family")
any(taxa_sums(Americans_Counts_physeq) == 0)
taxa_names(Americans_Counts_physeq) <- paste(tax_table(Americans_Counts_physeq)[,"Superkingdom"], tax_table(Americans_Counts_physeq)[,"Kingdom"], tax_table(Americans_Counts_physeq)[,"Phylum"], tax_table(Americans_Counts_physeq)[,"Class"], tax_table(Americans_Counts_physeq)[,"Order"], tax_table(Americans_Counts_physeq)[,"Family"], sep="_")

merged_phylo_counts=merge_phyloseq(AllREadsSE_Indo_Counts_physeq, Mali_Counts_physeq, control_100K_Counts_physeq, Americans_Counts_physeq)
a=subset_samples(merged_phylo_counts, SamplePop == "Netherlands")
any(taxa_sums(merged_phylo_counts) == 0)
# TRUE
merged_phylo_counts <- prune_taxa(taxa_sums(merged_phylo_counts) > 0, merged_phylo_counts)
any(taxa_sums(pruned_data) == 0)







tax_table(Mali_Counts_physeq)[,"Superkingdom"] = paste(tax_table(Mali_Counts_physeq)[,"Superkingdom"], tax_table(Mali_Counts_physeq)[,"Family"], sep="_")
tax_table(Mali_Counts_physeq)[,"Superkingdom"] <- gsub("Bacteria_$", "Bacteria_unclassified", tax_table(Mali_Counts_physeq)[,"Superkingdom"])
tax_table(Mali_Counts_physeq)[,"Superkingdom"] <- gsub("Eukaryota_$", "Eukaryota_unclassified", tax_table(Mali_Counts_physeq)[,"Superkingdom"])
tax_table(Mali_Counts_physeq)[,"Superkingdom"] <- gsub("Viruses_$", "Viruses_unclassified", tax_table(Mali_Counts_physeq)[,"Superkingdom"])

# For some reason, top is actually top + 1, so here it would be 20
aggregated_phyloCounts <- aggregate_top_taxa(Mali_Counts_physeq, "Superkingdom", top = 20)
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

PaletteArchaea = colorRampPalette(c("grey"))(1)
PaletteBacteria = colorRampPalette(c("#023858","#74a9cf"))(15)
PaletteEukaryote = colorRampPalette(c("#fd8d3c","#800026"))(2)
PaletteOther = colorRampPalette(c("black"))(1)
#PaletteUnk = colorRampPalette(c("grey"))(1)
PaletteVirus = colorRampPalette(c("#78c679","#006837"))(2)

Merged_Palette <- c(PaletteArchaea,PaletteBacteria,PaletteEukaryote,PaletteOther,PaletteVirus)

phyloseq::plot_bar(relative_phyloCounts, fill = "Superkingdom") +
  geom_bar(aes(fill = Superkingdom), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ SamplePop, scales = "free") + scale_fill_manual(values=Merged_Palette) +
  theme(panel.background = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())




