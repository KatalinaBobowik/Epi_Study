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
library(vegan)

# set ggplot colour theme to white
theme_set(theme_bw())

# set up directories
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/EpiStudy/"
outputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/Controls/"

# QC -----------------

# make rarefaction curves to look at the overall diversity of the control samples

raw_CCMetagen_data <- read.csv(paste0(inputdir,"Controls_NoFiltering_Counts.csv"),check.names=FALSE)

# Data preprocessing ------------------------------------

# add a new column containing family names and superkingdom
raw_CCMetagen_data$SuperKFamily <- paste(raw_CCMetagen_data$Superkingdom, raw_CCMetagen_data$Species, sep="_")

# Separate species' abundances and taxonomy columns
abund_raw_noSingletons <- as.matrix(raw_CCMetagen_data[,c(1:192)])

# get rid of anything with 5 counts or less
abund_raw_noSingletons[abund_raw_noSingletons==1]<-0
abund_raw_noSingletons[abund_raw_noSingletons==2]<-0
abund_raw_noSingletons[abund_raw_noSingletons==3]<-0
abund_raw_noSingletons[abund_raw_noSingletons==4]<-0
abund_raw_noSingletons[abund_raw_noSingletons==5]<-0

pdf(paste0(outputdir,"AllSickaandHealthyControls_Rarefaction.pdf"), width=15)
rarecurve(t(otu_table(abund_raw_noSingletons, taxa_are_rows = TRUE)), step=20, label=T, xlab="Counts",ylab="Number Species")
dev.off()

# rarefaction curve for just control sample
diseaseStatus=read.table(paste0(inputdir,"ControlSamples_DiseaseStatus.txt"))
colnames(diseaseStatus)=c("Samples","diseaseStatus")
controlSamples=as.character(diseaseStatus[grep("CSF", diseaseStatus$diseaseStatus),"Samples"])

pdf(paste0(outputdir,"ControlSamplesOnly_Rarefaction.pdf"), width=15)
rarecurve(t(otu_table(abund_raw_noSingletons[,controlSamples], taxa_are_rows = TRUE)), step=20, label=T, xlab="Counts",ylab="Number Species")
dev.off()

# Phyloseq pipeline ------------------------------

# Import data and convert to a phyloseq object
raw_CCMetagen_data <- read.csv(paste0(inputdir,"controlDataset_merged.csv"),check.names=FALSE)
 
# remove plants from the analysis
raw_CCMetagen_data=raw_CCMetagen_data[-which(raw_CCMetagen_data$Kingdom=="Viridiplantae"),]
raw_CCMetagen_data=raw_CCMetagen_data[-which(raw_CCMetagen_data$Phylum %in% c("Bacillariophyta","Cnidaria","Echinodermata","Mollusca")),]

# add a new column containing family names and superkingdom
raw_CCMetagen_data$SuperKFamily <- paste(raw_CCMetagen_data$Superkingdom, raw_CCMetagen_data$Family, sep="_")
# add an 'unclassified' string at the end of these to avoid confusion
raw_CCMetagen_data <- data.frame(lapply(raw_CCMetagen_data, function(x) {sub("Bacteria_$", "Bacteria_unclassified", x)}))
raw_CCMetagen_data <- data.frame(lapply(raw_CCMetagen_data, function(x) {sub("Eukaryota_$", "Eukaryota_unclassified", x)}))
raw_CCMetagen_data <- data.frame(lapply(raw_CCMetagen_data, function(x) {sub("Viruses_$", "Viruses_unclassified", x)}))

# one taxa, unk_sk_unk_f, comes from the study M10823.1 and is a Transposon Tn28 beta-lactmase gene
# which are found in bacteria. I'll change the name from "unk_sk_unk_f" to Bacteria_unclassified 
raw_CCMetagen_data$SuperKFamily[which(raw_CCMetagen_data$SuperKFamily %in% c("unk_sk_unk_f"))]="Bacteria_unclassified"
raw_CCMetagen_data$SuperKFamily=droplevels(raw_CCMetagen_data$SuperKFamily)

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
# "Bacteria_unclassified"  "Eukaryota_unclassified" "Eukaryota_unk_f"

# get first instance of unclassified eukaryotes and merge with other unclassified eukaryotes 
CCMetagen_data[grep("Eukaryota_unk_f|Eukaryota_unclassified",CCMetagen_data$Family)[1],4:ncol(CCMetagen_data)]=colSums(CCMetagen_data[grep("Eukaryota_unk_f|Eukaryota_unclassified",CCMetagen_data$Family),4:ncol(CCMetagen_data)])
# delete the rest
CCMetagen_data=CCMetagen_data[-grep("Eukaryota_unk_f|Eukaryota_unclassified",CCMetagen_data$Family)[2:length(grep("Eukaryota_unk_f|Eukaryota_unclassified",CCMetagen_data$Family))],]

# do this for bacteria
# get first instance of unclassified bacteria and merge with other unclassified bacteria
CCMetagen_data[grep("Bacteria_unclassified",CCMetagen_data$Family)[1],4:ncol(CCMetagen_data)]=colSums(CCMetagen_data[grep("Bacteria_unclassified",CCMetagen_data$Family),4:ncol(CCMetagen_data)])
# delete the rest
CCMetagen_data=CCMetagen_data[-grep("Bacteria_unclassified",CCMetagen_data$Family)[2:length(grep("Bacteria_unclassified",CCMetagen_data$Family))],]
CCMetagen_data=droplevels(CCMetagen_data)

# only get control samples
CCMetagen_data_Controls=CCMetagen_data[,c("Superkingdom","Kingdom","Family",controlSamples,"Genus","Species")]
# save this df
write.table(CCMetagen_data_Controls, file=paste0(outputdir,"ControlDF.txt"))

# Convert to PhyloSeq object ------------------------------------------

taxa_raw <- as.matrix(CCMetagen_data_Controls[,c("Superkingdom","Kingdom","Family")])
rownames(taxa_raw) <- taxa_raw[,"Family"]
abund_raw <- as.matrix(CCMetagen_data_Controls[,-which(colnames(CCMetagen_data_Controls) %in% c("Superkingdom","Kingdom","Family","Genus","Species"))])
rownames(abund_raw) <- CCMetagen_data_Controls[,which(names(CCMetagen_data_Controls)=="Family")]

# convert to Phyloseq object
tax = tax_table(taxa_raw)
taxa = otu_table(abund_raw, taxa_are_rows = TRUE)

# save taxa table
write.table(taxa, file=paste0(outputdir,"Control_Taxa_DF.txt"))

CCMeta_physeq_Controls = phyloseq(taxa, tax)
save(CCMeta_physeq_Controls, file = paste0(outputdir, "CCMeta_physeq_Controls.Rda"))

TopNOTUs <- names(sort(taxa_sums(CCMeta_physeq_Controls), TRUE)[1:20])
TopFamilies <- prune_taxa(TopNOTUs, CCMeta_physeq_Controls)
plot_bar(TopFamilies, fill = "Family")

p = plot_bar(TopFamilies, fill="Family")
# set colour palette
families=levels(p$data$Family)
# get number of families in each kingdom
table(sapply(strsplit(families, "[_.]"), `[`, 1))
 # Bacteria Eukaryota   
 #    14        6       
p$data$Family <- factor(p$data$Family) 

PaletteBacteria = colorRampPalette(c("#023858","#74a9cf"))(14)
PaletteEukaryote = colorRampPalette(c("#fd8d3c","#800026"))(6)
#PaletteVirus = colorRampPalette(c("#78c679","#006837"))(2)

Merged_Palette <- c(PaletteBacteria,PaletteEukaryote)

fig <- p + scale_fill_manual(values=Merged_Palette) +
  geom_bar(aes(fill=Family), stat="identity", position="stack") +
  guides(fill=guide_legend(ncol=2))

pdf(paste0(outputdir,"Allpathogens_Controls.pdf"), width=15)
fig
dev.off()

# get plot of frequency of each pathogen
melted_taxa=melt(taxa)
colnames(melted_taxa)=c("Pathogen","samples","RPM")
#melted_taxa$samples=gsub("\\..*","",melted_taxa$samples)
pdf(paste0(outputdir,"pathogensByIsland_plasmodium_Boxplot.pdf"))
ggboxplot(melted_taxa, x = "Pathogen", y = "RPM", fill="Island", add=c("boxplot"),add.params = c(list(fill = "white"), list(width=0.05)))
dev.off()

# Repeat this but for all taxa included ----------------

raw_CCMetagen_data <- read.csv(paste0(inputdir,"controlDataset_merged.csv"),check.names=FALSE)
 
# remove plants from the analysis
raw_CCMetagen_data=raw_CCMetagen_data[-which(raw_CCMetagen_data$Kingdom=="Viridiplantae"),]
raw_CCMetagen_data=raw_CCMetagen_data[-which(raw_CCMetagen_data$Phylum %in% c("Bacillariophyta","Cnidaria","Echinodermata","Mollusca")),]

# only get control samples
diseaseStatus=read.table(paste0(inputdir,"ControlSamples_DiseaseStatus.txt"))
colnames(diseaseStatus)=c("Samples","diseaseStatus")
controlSamples=as.character(diseaseStatus[grep("CSF", diseaseStatus$diseaseStatus),"Samples"])

CCMetagen_data_Controls=raw_CCMetagen_data[,c(controlSamples,"Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species")]
# save this df
write.table(CCMetagen_data_Controls, file=paste0(outputdir,"ControlDF.txt"))

# Convert to PhyloSeq object ------------------------------------------

taxa_raw <- as.matrix(CCMetagen_data_Controls[,c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species")])
#rownames(taxa_raw) <- taxa_raw[,"Family"]
abund_raw <- as.matrix(CCMetagen_data_Controls[,-which(colnames(CCMetagen_data_Controls) %in% c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species"))])
#rownames(abund_raw) <- CCMetagen_data_Controls[,which(names(CCMetagen_data_Controls)=="Family")]

# convert to Phyloseq object
tax = tax_table(taxa_raw)
taxa = otu_table(abund_raw, taxa_are_rows = TRUE)

# save taxa table
write.table(taxa, file=paste0(outputdir,"Control_Taxa_DF.txt"))

CCMeta_physeq_Controls = phyloseq(taxa, tax)
save(CCMeta_physeq_Controls, file = paste0(outputdir, "CCMeta_physeq_Controls_AllTaxaRanks.Rda"))

TopNOTUs <- names(sort(taxa_sums(CCMeta_physeq_Controls), TRUE)[1:20])
TopFamilies <- prune_taxa(TopNOTUs, CCMeta_physeq_Controls)
plot_bar(TopFamilies, fill = "Family")

p = plot_bar(TopFamilies, fill="Family")
# set colour palette
families=levels(p$data$Family)
# get number of families in each kingdom
table(sapply(strsplit(families, "[_.]"), `[`, 1))
 # Bacteria Eukaryota   
 #    14        6       
p$data$Family <- factor(p$data$Family) 

PaletteBacteria = colorRampPalette(c("#023858","#74a9cf"))(14)
PaletteEukaryote = colorRampPalette(c("#fd8d3c","#800026"))(6)
#PaletteVirus = colorRampPalette(c("#78c679","#006837"))(2)

Merged_Palette <- c(PaletteBacteria,PaletteEukaryote)

fig <- p + scale_fill_manual(values=Merged_Palette) +
  geom_bar(aes(fill=Family), stat="identity", position="stack") +
  guides(fill=guide_legend(ncol=2))

pdf(paste0(outputdir,"Allpathogens_Controls.pdf"), width=15)
fig
dev.off()

# get number of empty rows in each rank
taxa=c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
for (level in taxa){
	length=length(which(CCMetagen_data_Controls[,level]==""))
	print(c(level,length))
}

tax_table(CCMeta_physeq_Controls)
for (level in taxa){
	length=length(which(tax_table(CCMeta_physeq_Controls)[,level]==""))
	perc_missing=length/nrow(tax_table(CCMeta_physeq_Controls))
	print(c(level,length, perc_missing))
}


