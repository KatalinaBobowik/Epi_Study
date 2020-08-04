# created by KSB, 20.01.19
# Script created to inspect where reads are mapping to. 


# load packages
require(ggplot2)
require(RColorBrewer)
library(dplyr)
library(plyr)
library(phyloseq)
library(reshape2)

# set ggplot colour theme to white
theme_set(theme_bw())

# set up directories
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/EpiStudy/"
outputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/Controls"

# Import data and convert to a phyloseq object
# raw_CCMetagen_data <- read.csv(paste0(inputdir,"IndoUnmapped_species_table_Helminths.csv"),check.names=FALSE)
raw_CCMetagen_data <- read.csv(paste0(inputdir,"ControlUnmapped_species_table_RPM.csv"),check.names=FALSE)
raw_CCMetagen_data <- read.csv("/Users/katalinabobowik/Desktop/ControlUnmapped_species_table_noRepeats_noRNA_RPM.csv",check.names=FALSE)

# remove plants from the analysis
raw_CCMetagen_data=raw_CCMetagen_data[-which(raw_CCMetagen_data$Kingdom=="Viridiplantae"),]
raw_CCMetagen_data=raw_CCMetagen_data[-which(raw_CCMetagen_data$Phylum %in% c("Bacillariophyta","Cnidaria","Echinodermata","Mollusca")),]

# add a new column containing family names and superkingdom
raw_CCMetagen_data$SuperKFamily <- paste(raw_CCMetagen_data$Superkingdom, raw_CCMetagen_data$Family, sep="_")
# add an 'unclassified' string at the end of these to avoid confusion
raw_CCMetagen_data <- data.frame(lapply(raw_CCMetagen_data, function(x) {sub("Bacteria_$", "Bacteria_unclassified", x)}))
raw_CCMetagen_data <- data.frame(lapply(raw_CCMetagen_data, function(x) {sub("Eukaryota_$", "Eukaryota_unclassified", x)}))
raw_CCMetagen_data <- data.frame(lapply(raw_CCMetagen_data, function(x) {sub("Viruses_$", "Viruses_unclassified", x)}))

# delete taxonomic ranks for which there are multiple taxa merged into 'Eukaryota_unclassified' or 'Bacteria_unclassified' - (delete Phylum, Class, Order and the previous Family)
CCMetagen_data <-raw_CCMetagen_data[,-which(names(raw_CCMetagen_data) %in% c("Phylum","Class","Order","Family"))]
# convert the abundances to numeric
CCMetagen_data[-which(colnames(CCMetagen_data) %in% c("Superkingdom","Kingdom","Genus","Species","SuperKFamily"))] = mutate_all(CCMetagen_data[,-which(colnames(CCMetagen_data) %in% c("Superkingdom","Kingdom","Genus","Species","SuperKFamily"))], function(x) as.numeric(as.character(x)))
# aggregate the unclassified taxa
CCMetagen_data <- aggregate(. ~ Superkingdom+Kingdom+SuperKFamily,CCMetagen_data, sum)
# rename the new family column
colnames(CCMetagen_data)[which(names(CCMetagen_data)=="SuperKFamily")] <- "Family"
# find out which names in family are duplicated (if duplicated, this will give us troubles later on)
as.character(unique(CCMetagen_data[which(duplicated(CCMetagen_data$Family)),]$Family))
# "Eukaryota_unclassified"

# get first instance of unclassified eukaryotes and merge with other unclassified eukaryotes 
CCMetagen_data[grep("Eukaryota_unk_f|Eukaryota_unclassified",CCMetagen_data$Family)[1],4:ncol(CCMetagen_data)]=colSums(CCMetagen_data[grep("Eukaryota_unk_f|Eukaryota_unclassified",CCMetagen_data$Family),4:ncol(CCMetagen_data)])
# delete the rest
CCMetagen_data=CCMetagen_data[-grep("Eukaryota_unk_f|Eukaryota_unclassified",CCMetagen_data$Family)[2:length(grep("Eukaryota_unk_f|Eukaryota_unclassified",CCMetagen_data$Family))],]
CCMetagen_data=droplevels(CCMetagen_data)

# Convert to PhyloSeq object ------------------------------------------

taxa_raw <- as.matrix(CCMetagen_data[,c("Superkingdom","Kingdom","Family")])
rownames(taxa_raw) <- taxa_raw[,"Family"]
abund_raw <- as.matrix(CCMetagen_data[,-which(colnames(CCMetagen_data) %in% c("Superkingdom","Kingdom","Family","Genus","Species"))])
rownames(abund_raw) <- CCMetagen_data[,which(names(CCMetagen_data)=="Family")]

# convert to Phyloseq object
tax = tax_table(taxa_raw)
taxa = otu_table(abund_raw, taxa_are_rows = TRUE)

CCMeta_physeq = phyloseq(taxa, tax)
plot_bar(CCMeta_physeq, fill = "Superkingdom")

TopNOTUs <- names(sort(taxa_sums(CCMeta_physeq), TRUE)[1:16])
TopFamilies <- prune_taxa(TopNOTUs, CCMeta_physeq)
plot_bar(TopFamilies, fill = "Family")

p = plot_bar(TopFamilies, fill="Family")
# set colour palette
families=levels(p$data$Family)
# get number of families in each kingdom
table(sapply(strsplit(families, "[_.]"), `[`, 1))
 # Bacteria Eukaryota   Viruses 
 #    1        23          3 
p$data$Family <- factor(p$data$Family) 

PaletteBacteria = colorRampPalette(c("#023858","#74a9cf"))(1)
PaletteEukaryote = colorRampPalette(c("#fd8d3c","#800026"))(13)
PaletteVirus = colorRampPalette(c("#78c679","#006837"))(2)

Merged_Palette <- c(PaletteBacteria,PaletteEukaryote,PaletteVirus)

fig <- p + scale_fill_manual(values=Merged_Palette) +
  geom_bar(aes(fill=Family), stat="identity", position="stack") +
  guides(fill=guide_legend(ncol=2))

pdf(paste0(outputdir,"Allpathogens.pdf"), width=15)
fig
dev.off()

# save OTU table
scaledOTUs = sapply(1:123, function(x) otu_table(TopFamilies)[,x]/sum(otu_table(TopFamilies)[,x]))
colnames(scaledOTUs) = colnames(otu_table(TopFamilies))
rownames(scaledOTUs) = rownames(otu_table(TopFamilies))
write.table(scaledOTUs, file=paste0(outputdir,"scaledOTUs.txt"))

# get plot of frequency of each pathogen
melted_taxa=melt(taxa)
colnames(melted_taxa)=c("Pathogen","samples","RPM")
#melted_taxa$samples=gsub("\\..*","",melted_taxa$samples)
pdf(paste0(outputdir,"pathogensByIsland_plasmodium_Boxplot.pdf"))
ggboxplot(melted_taxa, x = "Pathogen", y = "RPM", fill="Island", add=c("boxplot"),add.params = c(list(fill = "white"), list(width=0.05)))
dev.off()

