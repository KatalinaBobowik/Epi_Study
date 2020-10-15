# created by KSB, 20.01.19
# Create a PhyloSeq object for the sub-sampled, 100K Indonesian data

# load packages
require(ggplot2)
require(RColorBrewer)
library(dplyr)
library(plyr)
library(phyloseq)
library(reshape2)
library(ggpubr)
library(vegan)
library(ggfortify)

# set up directories
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/100KSamples/"
outputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/100KSamples/"
controldir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/Controls/"

#########################
##### Indo RPM Data #####
#########################

# Indonesian Count
load(paste0(inputdir,"subSampled_100K_Indo_Counts_physeq.Rda"))
# Indonesian RPM
load(paste0(inputdir,"subSampled_100K_Indo_RPM_physeq.Rda"))
# control counts
load(paste0(inputdir,"control_100K_Counts_physeq.Rda"))
# control RPM
load(paste0(inputdir,"control_100K_RPM_physeq.Rda"))

# Take out singletons for controls
otu_table(control_100K_RPM_physeq)[otu_table(control_100K_Counts_physeq)<=1]<-0
# remove 0s
control_100K_RPM_physeq <- prune_taxa(taxa_sums(control_100K_RPM_physeq) > 0, control_100K_RPM_physeq)

# Mapping reads by RPM ---------------------

control_100K_RPM_physeq=subset_taxa(control_100K_RPM_physeq, (Kingdom!="Viridiplantae"))
control_100K_RPM_physeq=subset_taxa(control_100K_RPM_physeq, (Phylum!="Chordata"))
control_100K_RPM_physeq <- prune_taxa(taxa_sums(control_100K_RPM_physeq) > 0, control_100K_RPM_physeq)

# add a new column containing family names and superkingdom
tax_table(control_100K_RPM_physeq)[,"Superkingdom"] = paste(tax_table(control_100K_RPM_physeq)[,"Superkingdom"], tax_table(control_100K_RPM_physeq)[,"Family"], sep="_")
tax_table(control_100K_RPM_physeq)[,"Superkingdom"] <- gsub("Bacteria_$", "Bacteria_unclassified", tax_table(control_100K_RPM_physeq)[,"Superkingdom"])
tax_table(control_100K_RPM_physeq)[,"Superkingdom"] <- gsub("Eukaryota_$", "Eukaryota_unclassified", tax_table(control_100K_RPM_physeq)[,"Superkingdom"])
tax_table(control_100K_RPM_physeq)[,"Superkingdom"] <- gsub("Viruses_$", "Viruses_unclassified", tax_table(control_100K_RPM_physeq)[,"Superkingdom"])

# Plot
TopNOTUs <- names(sort(taxa_sums(control_100K_RPM_physeq), TRUE)[1:20])
TopFamilies <- prune_taxa(TopNOTUs, control_100K_RPM_physeq)
p=plot_bar(TopFamilies, fill = "Superkingdom")

# set colour palette
families=levels(p$data$Superkingdom)
# get number of families in each kingdom
table(sapply(strsplit(families, "[_.]"), `[`, 1))
 # Bacteria Eukaryota  Viruses
 #    9         2      1      

PaletteBacteria = colorRampPalette(c("#023858","#74a9cf"))(4)
PaletteEukaryote = colorRampPalette(c("#fd8d3c","#800026"))(8)
PaletteVirus = colorRampPalette(c("#78c679","#006837"))(1)

Merged_Palette <- c(PaletteBacteria,PaletteEukaryote,PaletteVirus)

fig <- p + scale_fill_manual(values=Merged_Palette) +
  geom_bar(aes(fill=Superkingdom), stat="identity", position="stack") +
  guides(fill=guide_legend(ncol=2)) + theme_bw()

pdf(paste0(outputdir,"Barplot_NoPlanst_NoChordates.pdf"), width=10)
fig
dev.off()

##########
# INdo ###
##########

# Take out singletons for controls
otu_table(subSampled_100K_Indo_RPM_physeq)[otu_table(subSampled_100K_Indo_Counts_physeq)<=1]<-0
# remove 0s
subSampled_100K_Indo_RPM_physeq <- prune_taxa(taxa_sums(subSampled_100K_Indo_RPM_physeq) > 0, subSampled_100K_Indo_RPM_physeq)

# Mapping reads by RPM ---------------------

subSampled_100K_Indo_RPM_physeq=subset_taxa(subSampled_100K_Indo_RPM_physeq, (Kingdom!="Viridiplantae"))
subSampled_100K_Indo_RPM_physeq=subset_taxa(subSampled_100K_Indo_RPM_physeq, (Phylum!="Chordata"))
subSampled_100K_Indo_RPM_physeq <- prune_taxa(taxa_sums(subSampled_100K_Indo_RPM_physeq) > 0, subSampled_100K_Indo_RPM_physeq)

# add a new column containing family names and superkingdom
tax_table(subSampled_100K_Indo_RPM_physeq)[,"Superkingdom"] = paste(tax_table(subSampled_100K_Indo_RPM_physeq)[,"Superkingdom"], tax_table(subSampled_100K_Indo_RPM_physeq)[,"Family"], sep="_")
tax_table(subSampled_100K_Indo_RPM_physeq)[,"Superkingdom"] <- gsub("Bacteria_$", "Bacteria_unclassified", tax_table(subSampled_100K_Indo_RPM_physeq)[,"Superkingdom"])
tax_table(subSampled_100K_Indo_RPM_physeq)[,"Superkingdom"] <- gsub("Eukaryota_$", "Eukaryota_unclassified", tax_table(subSampled_100K_Indo_RPM_physeq)[,"Superkingdom"])
tax_table(subSampled_100K_Indo_RPM_physeq)[,"Superkingdom"] <- gsub("Viruses_$", "Viruses_unclassified", tax_table(subSampled_100K_Indo_RPM_physeq)[,"Superkingdom"])

# Plot
TopNOTUs <- names(sort(taxa_sums(subSampled_100K_Indo_RPM_physeq), TRUE)[1:20])
TopFamilies <- prune_taxa(TopNOTUs, subSampled_100K_Indo_RPM_physeq)
p=plot_bar(TopFamilies, fill = "Superkingdom")

# set colour palette
families=levels(p$data$Superkingdom)
# get number of families in each kingdom
table(sapply(strsplit(families, "[_.]"), `[`, 1))
 # Bacteria Eukaryota  Viruses
 #    9         2      1      

PaletteBacteria = colorRampPalette(c("#023858","#74a9cf"))(7)
PaletteEukaryote = colorRampPalette(c("#fd8d3c","#800026"))(4)
PaletteUnk = colorRampPalette(c("burlywood"))(1)
PaletteVirus = colorRampPalette(c("#78c679","#006837"))(2)

Merged_Palette <- c(PaletteBacteria,PaletteEukaryote,PaletteUnk,PaletteVirus)

fig <- p + scale_fill_manual(values=Merged_Palette) +
  geom_bar(aes(fill=Superkingdom), stat="identity", position="stack") +
  guides(fill=guide_legend(ncol=2)) + theme_bw()

pdf(paste0(outputdir,"Barplot_NoPlanst_NoChordates_indo.pdf"), width=10)
fig
dev.off()

# without Plasmodium
plasmodium_rows=grep("Plasmodiidae",tax_table(subSampled_100K_Indo_RPM_physeq)[,"Superkingdom"])
plasmodium_otu=rownames(tax_table(subSampled_100K_Indo_RPM_physeq)[plasmodium_rows,])
allTaxa = taxa_names(subSampled_100K_Indo_RPM_physeq)
allTaxa <- allTaxa[!(allTaxa %in% plasmodium_otu)]
subSampled_100K_Indo_RPM_physeq_noPlasmo=prune_taxa(allTaxa, subSampled_100K_Indo_RPM_physeq)

# Plot
TopNOTUs <- names(sort(taxa_sums(subSampled_100K_Indo_RPM_physeq_noPlasmo), TRUE)[1:20])
TopFamilies <- prune_taxa(TopNOTUs, subSampled_100K_Indo_RPM_physeq_noPlasmo)
p=plot_bar(TopFamilies, fill = "Superkingdom")

# set colour palette
families=levels(p$data$Superkingdom)
# get number of families in each kingdom
table(sapply(strsplit(families, "[_.]"), `[`, 1))
 # Bacteria Eukaryota  Viruses
 #    9         2      1      

PaletteBacteria = colorRampPalette(c("#023858","#74a9cf"))(9)
PaletteEukaryote = colorRampPalette(c("#fd8d3c","#800026"))(3)
PaletteUnk = colorRampPalette(c("burlywood"))(1)
PaletteVirus = colorRampPalette(c("#78c679","#006837"))(2)

Merged_Palette <- c(PaletteBacteria,PaletteEukaryote,PaletteUnk,PaletteVirus)

fig <- p + scale_fill_manual(values=Merged_Palette) +
  geom_bar(aes(fill=Superkingdom), stat="identity", position="stack") +
  guides(fill=guide_legend(ncol=2)) + theme_bw()

pdf(paste0(outputdir,"Barplot_NoPlants_NoChordates_noPlasmo_indo.pdf"), width=10)
fig
dev.off()


