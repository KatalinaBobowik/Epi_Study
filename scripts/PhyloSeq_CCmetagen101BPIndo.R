# created by KSB, 20.01.19
# Script created to inspect where reads are mapping to. 

# load packages
require(ggplot2)
require(RColorBrewer)
library(dplyr)
library(plyr)
library(phyloseq)
library(reshape2)
library(ggpubr)

# set ggplot colour theme to white
theme_set(theme_bw())

# set up directories
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/EpiStudy/"
outputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/"

# Import data and convert to a phyloseq object
raw_CCMetagen_data <- read.csv(paste0(inputdir,"101BPIndo_Unmapped_species_table_noRepeats_noRNA_RPM.csv"),check.names=FALSE)

# Data preprocessing ------------------------------------

# remove plants from the analysis
raw_CCMetagen_data=raw_CCMetagen_data[-which(raw_CCMetagen_data$Kingdom %in% c("Viridiplantae","Fungi")),]
raw_CCMetagen_data=raw_CCMetagen_data[-which(raw_CCMetagen_data$Phylum %in% c("Bacillariophyta","Cnidaria","Echinodermata","Mollusca")),]
raw_CCMetagen_data=raw_CCMetagen_data[-which(raw_CCMetagen_data$Order %in% "Albuginales"),]

# remove all reference samples that had a low number of accession numbers
# lowAccessiontaxa=read.table(paste0(outputdir,"lowAccessiontaxa.txt"))
# lowAccessiontaxa=as.vector(lowAccessiontaxa[,1])
# raw_CCMetagen_data=raw_CCMetagen_data[-which(raw_CCMetagen_data$Phylum %in% lowAccessiontaxa),]

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
CCMetagen_data[grep("Eukaryota_unk_f|Eukaryota_unclassified",CCMetagen_data$Family)[1],4:128]=colSums(CCMetagen_data[grep("Eukaryota_unk_f|Eukaryota_unclassified",CCMetagen_data$Family),4:128])
# delete the rest
CCMetagen_data=CCMetagen_data[-grep("Eukaryota_unk_f|Eukaryota_unclassified",CCMetagen_data$Family)[2:length(grep("Eukaryota_unk_f|Eukaryota_unclassified",CCMetagen_data$Family))],]
CCMetagen_data=droplevels(CCMetagen_data)


# Convert to PhyloSeq object ------------------------------------------

# Separate species' abundances and taxonomy columns
taxa_raw <- as.matrix(CCMetagen_data[,c("Superkingdom","Kingdom","Family")])
rownames(taxa_raw) <- taxa_raw[,"Family"]
abund_raw <- as.matrix(CCMetagen_data[,-which(colnames(CCMetagen_data) %in% c("Superkingdom","Kingdom","Family","Genus","Species"))])
rownames(abund_raw) <- CCMetagen_data[,which(names(CCMetagen_data)=="Family")]

# convert to Phyloseq object
tax = tax_table(taxa_raw)
taxa = otu_table(abund_raw, taxa_are_rows = TRUE)

# remove unclassified taxa. These have a low number of reads mapping to them and add no extra information
taxa=taxa[-which(rownames(taxa) %in% c("unk_sk_unk_f")),]

# exploratory stuff on species per island

df=melt(taxa)
df$island=NA
df$island[grep("MPI",df$Var2)]="MPI"
df$island[grep("SMB",df$Var2)]="SMB"
df$island[grep("MTW",df$Var2)]="MTW"

pdf(paste0(outputdir,"pathogensByIsland_noRepeatsNolowAccession.pdf"))
ggplot(df, aes(x=island, y=value, fill=Var1)) + geom_bar(width = 1, stat = "identity") + labs(y="RPM", x = "Island")
dev.off()

# The high number of reads mapping to Plasmodium in Sumba is really just due to one sample
taxa_noLowReads=taxa[which(rownames(taxa) %in% c("Eukaryota_Plasmodiidae")),]
melted_taxa=melt(taxa_noLowReads)
colnames(melted_taxa)=c("Pathogen","Island","RPM")
melted_taxa$Island=gsub("\\..*","",melted_taxa$Island)
melted_taxa$Pathogen=gsub("Bacteria_","",melted_taxa$Pathogen) %>% gsub("Eukaryota_","",.) %>% gsub("Viruses_","",.)
pdf(paste0(outputdir,"pathogensByIsland_plasmodium_Boxplot.pdf"))
ggboxplot(melted_taxa, x = "Island", y = "RPM", fill="Island", add=c("boxplot"),add.params = c(list(fill = "white"), list(width=0.05)))
dev.off()

# back to Phyloseq pipeline -----------------------------

CCMeta_physeq = phyloseq(taxa, tax)
# TopNOTUs <- names(sort(taxa_sums(CCMeta_physeq), TRUE))
# TopFamilies <- prune_taxa(TopNOTUs, CCMeta_physeq)
p=plot_bar(CCMeta_physeq, fill = "Family")

# set colour palette
families=levels(p$data$Family)
# get number of families in each kingdom
table(sapply(strsplit(families, "[_.]"), `[`, 1))
# Bacteria Eukaryota   Viruses 
#      18         9         4 
p$data$Family <- factor(p$data$Family) 

PaletteBacteria = colorRampPalette(c("#023858","#74a9cf"))(18)
PaletteEukaryote = colorRampPalette(c("#fd8d3c","#800026"))(9)
PaletteVirus = colorRampPalette(c("#78c679","#006837"))(4)

Merged_Palette <- c(PaletteBacteria,PaletteEukaryote,PaletteVirus)

fig <- p + scale_fill_manual(values=Merged_Palette) +
  geom_bar(aes(fill=Family), stat="identity", position="stack") +
  guides(fill=guide_legend(ncol=2))

pdf(paste0(outputdir,"Allpathogens.pdf"), width=15)
fig
dev.off()

# save OTU table
scaledOTUs = sapply(1:123, function(x) otu_table(CCMeta_physeq)[,x]/sum(otu_table(CCMeta_physeq)[,x]))
colnames(scaledOTUs) = colnames(otu_table(CCMeta_physeq))
rownames(scaledOTUs) = rownames(otu_table(CCMeta_physeq))
write.table(scaledOTUs, file=paste0(outputdir,"scaledOTUs.txt"))

# plot trends for each taxa
melted_taxa=melt(taxa)
colnames(melted_taxa)=c("Pathogen","Island","RPM")
melted_taxa$Island=gsub("\\..*","",melted_taxa$Island)
melted_taxa$Pathogen=gsub("Bacteria_","",melted_taxa$Pathogen) %>% gsub("Eukaryota_","",.) %>% gsub("Viruses_","",.)
pdf(paste0(outputdir,"pathogensByIsland_AllPathogens_Boxplot.pdf"))
ggboxplot(melted_taxa, x = "Island", y = "RPM", fill="Island", add=c("boxplot"),add.params = c(list(fill = "white"), list(width=0.05))) + facet_wrap(~Pathogen, scales = "free")
dev.off()

# without Plasmodium
taxa=taxa[-which(rownames(taxa) %in% "Eukaryota_Plasmodiidae"),]

CCMeta_physeq = phyloseq(taxa, tax)
p=plot_bar(CCMeta_physeq, fill = "Family")

PaletteBacteria = colorRampPalette(c("#023858","#74a9cf"))(18)
PaletteEukaryote = colorRampPalette(c("#fd8d3c","#800026"))(8)
PaletteVirus = colorRampPalette(c("#78c679","#006837"))(4)


Merged_Palette <- c(PaletteBacteria,PaletteEukaryote,PaletteVirus)

fig <- p + scale_fill_manual(values=Merged_Palette) +
  geom_bar(aes(fill=Family), stat="identity", position="stack") +
  guides(fill=guide_legend(ncol=2))

pdf(paste0(outputdir,"Allpathogens_NoPlasmo.pdf"), width=15)
fig
dev.off()
