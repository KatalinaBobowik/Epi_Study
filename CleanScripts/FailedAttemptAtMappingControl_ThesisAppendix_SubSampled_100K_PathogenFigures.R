# created by KSB, 20.01.19
# Create a PhyloSeq object for the sub-sampled, 100K Indonesian data

# load packages
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
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/EpiStudy/"
outputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/Controls/"

# load in data --------------------------
AllREadsSE_Indo_Counts <- read.csv(paste0(inputdir,"Controls_NoFiltering_Counts.csv"),check.names=FALSE)

# Separate species' abundances and taxonomy columns
taxa_raw <- as.matrix(AllREadsSE_Indo_Counts[,c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species")])
abund_raw <- as.matrix(AllREadsSE_Indo_Counts[,-which(colnames(AllREadsSE_Indo_Counts) %in% c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species"))])

# convert to Phyloseq object
tax = tax_table(taxa_raw)
taxa = otu_table(abund_raw, taxa_are_rows = TRUE)
AllREadsSE_Indo_Counts_physeq = phyloseq(taxa, tax)

# make this into a df and add to the Phloseq object
samples_df=data.frame(SampleName=colnames(otu_table(AllREadsSE_Indo_Counts_physeq)))
samples = sample_data(samples_df)
rownames(samples)=samples$SampleName
sample_data(AllREadsSE_Indo_Counts_physeq) <- samples

# get information on phyloseq objects
AllREadsSE_Indo_Counts_physeq

# Data processing ------------------

merged_phylo_counts = AllREadsSE_Indo_Counts_physeq

# Filter out singletons
otu_table(merged_phylo_counts)[otu_table(merged_phylo_counts)<=1]<-0
merged_phylo_counts <- prune_taxa(taxa_sums(merged_phylo_counts) > 0, merged_phylo_counts)

# Filter out Viridiplantae 
merged_phylo_counts <- subset_taxa(merged_phylo_counts, (Kingdom!="Viridiplantae"))
merged_phylo_counts <- prune_taxa(taxa_sums(merged_phylo_counts) > 0, merged_phylo_counts)

# Filter out Chordata
merged_phylo_counts <- subset_taxa(merged_phylo_counts, (Phylum!="Chordata"))
merged_phylo_counts <- prune_taxa(taxa_sums(merged_phylo_counts) > 0, merged_phylo_counts)

# Filter out Metazoa
merged_phylo_counts <- subset_taxa(merged_phylo_counts, (Kingdom!="Metazoa"))
merged_phylo_counts <- prune_taxa(taxa_sums(merged_phylo_counts) > 0, merged_phylo_counts)

# Filter out taxa which are unassigned at the Kingdom level
merged_phylo_counts <- subset_taxa(merged_phylo_counts, (Superkingdom!="unk_sk"))
merged_phylo_counts <- prune_taxa(taxa_sums(merged_phylo_counts) > 0, merged_phylo_counts)

# add a new column containing family names and superkingdom
tax_table(merged_phylo_counts)[,"Superkingdom"] = paste(tax_table(merged_phylo_counts)[,"Superkingdom"], tax_table(merged_phylo_counts)[,"Family"], sep="_")
tax_table(merged_phylo_counts)[,"Superkingdom"] <- gsub("Bacteria_$", "Bacteria_unclassified", tax_table(merged_phylo_counts)[,"Superkingdom"])
tax_table(merged_phylo_counts)[,"Superkingdom"] <- gsub("Eukaryota_$", "Eukaryota_unclassified", tax_table(merged_phylo_counts)[,"Superkingdom"])
tax_table(merged_phylo_counts)[,"Superkingdom"] <- gsub("Viruses_$", "Viruses_unclassified", tax_table(merged_phylo_counts)[,"Superkingdom"])

# As pointed out, we have a lot of taxa at the family level, and it would be hard to look over everything at once. Instead, we can focus on the most prevalent taxa and 
# highlight everything else in another colour. Here, I chose to highlight the top 20 taxa, since that's still representative while not being too visually exhausting.

aggregated_phyloCounts <- aggregate_top_taxa(merged_phylo_counts, "Superkingdom", top = 20)
# # transform to relative counts
# relative_phyloCounts <- microbiome::transform(aggregated_phyloCounts, "compositional")
# Remove weird extra family names added at the end of Superkingdom names
tax_table(aggregated_phyloCounts)[,"Superkingdom"] <- paste(sapply(strsplit(taxa_names(aggregated_phyloCounts), "[_.]"), `[`, 1), sapply(strsplit(taxa_names(aggregated_phyloCounts), "[_.]"), `[`, 2), sep="_")
# Change "Other_NA" to just "Other"
tax_table(aggregated_phyloCounts)[,"Superkingdom"][grep("Other", taxa_names(aggregated_phyloCounts))] = "Other"

# Plot
p=plot_bar(aggregated_phyloCounts, fill = "Superkingdom")

# set colour palette
families=levels(p$data$Superkingdom)
# get number of families in each kingdom
table(sapply(strsplit(families, "[_.]"), `[`, 1))

PaletteBacteria = colorRampPalette(c("#023858","#74a9cf"))(15)
PaletteEukaryote = colorRampPalette(c("#b20000","#ff4c4c"))(5)
PaletteOther = colorRampPalette(c("black"))(1)

Merged_Palette <- c(PaletteBacteria,PaletteEukaryote,PaletteOther)

pdf(paste0(outputdir,"CountsofTaxa_Barplot_Controls.pdf"), width=15)
phyloseq::plot_bar(aggregated_phyloCounts, fill = "Superkingdom") +
  geom_bar(aes(fill = Superkingdom), stat = "identity", position = "stack") +
  labs(x = "", y = "Count Abundance\n") + theme_bw(base_size = 15) +
  scale_fill_manual(values=Merged_Palette) + 
  theme(panel.background = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
dev.off()

# mean total for entire library
mean(colSums(otu_table(merged_phylo_counts)))
# [1] 13.48438


