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
lowAccessiontaxa=read.table(paste0(outputdir,"lowAccessiontaxa.txt"))
lowAccessiontaxa=as.vector(lowAccessiontaxa[,1])
raw_CCMetagen_data=raw_CCMetagen_data[-which(raw_CCMetagen_data$Phylum %in% lowAccessiontaxa),]

# add a new column containing family names and superkingdom
raw_CCMetagen_data$SuperKPhylum <- paste(raw_CCMetagen_data$Superkingdom, raw_CCMetagen_data$Phylum, sep="_")
# add an 'unclassified' string at the end of these to avoid confusion
raw_CCMetagen_data <- data.frame(lapply(raw_CCMetagen_data, function(x) {sub("Bacteria_$", "Bacteria_unclassified", x)}))
raw_CCMetagen_data <- data.frame(lapply(raw_CCMetagen_data, function(x) {sub("Eukaryota_$", "Eukaryota_unclassified", x)}))
raw_CCMetagen_data <- data.frame(lapply(raw_CCMetagen_data, function(x) {sub("Viruses_$", "Viruses_unclassified", x)}))

# delete taxonomic ranks for which there are multiple taxa merged into 'Eukaryota_unclassified' or 'Bacteria_unclassified' - (delete Phylum, Class, Order and the previous Family)
CCMetagen_data <-raw_CCMetagen_data[,-which(names(raw_CCMetagen_data) %in% c("Phylum","Class","Order","Family"))]
# convert the abundances to numeric
CCMetagen_data[-which(colnames(CCMetagen_data) %in% c("Superkingdom","Kingdom","Genus","Species","SuperKPhylum"))] = mutate_all(CCMetagen_data[,-which(colnames(CCMetagen_data) %in% c("Superkingdom","Kingdom","Genus","Species","SuperKPhylum"))], function(x) as.numeric(as.character(x)))
# aggregate the unclassified taxa
CCMetagen_data <- aggregate(. ~ Superkingdom+Kingdom+SuperKPhylum,CCMetagen_data, sum)
# rename the new family column
colnames(CCMetagen_data)[which(names(CCMetagen_data)=="SuperKPhylum")] <- "Phylum"
# find out which names in family are duplicated (if duplicated, this will give us troubles later on)
as.character(unique(CCMetagen_data[which(duplicated(CCMetagen_data$Phylum)),]$Phylum))
# character(0)

# Convert to PhyloSeq object ------------------------------------------

# Separate species' abundances and taxonomy columns
taxa_raw <- as.matrix(CCMetagen_data[,c("Superkingdom","Kingdom","Phylum")])
rownames(taxa_raw) <- taxa_raw[,"Phylum"]
abund_raw <- as.matrix(CCMetagen_data[,-which(colnames(CCMetagen_data) %in% c("Superkingdom","Kingdom","Phylum","Genus","Species"))])
rownames(abund_raw) <- CCMetagen_data[,which(names(CCMetagen_data)=="Phylum")]

# convert to Phyloseq object
tax = tax_table(taxa_raw)
taxa = otu_table(abund_raw, taxa_are_rows = TRUE)

# exploratory stuff on species per island

df=melt(taxa)
df$island=NA
df$island[grep("MPI",df$Var2)]="MPI"
df$island[grep("SMB",df$Var2)]="SMB"
df$island[grep("MTW",df$Var2)]="MTW"

pdf(paste0(outputdir,"pathogensByIsland_Phylum.pdf"))
ggplot(df, aes(x=island, y=value, fill=Var1)) + geom_bar(width = 1, stat = "identity") + labs(y="RPM", x = "Island")
dev.off()

melted_taxa=melt(taxa)
colnames(melted_taxa)=c("Pathogen","Island","RPM")
melted_taxa$Island=gsub("\\..*","",melted_taxa$Island)
melted_taxa$Pathogen=gsub("Bacteria_","",melted_taxa$Pathogen) %>% gsub("Eukaryota_","",.)
pdf(paste0(outputdir,"pathogensByIsland_noRepeatsNolowAccession_Boxplot_Phylum.pdf"))
ggboxplot(melted_taxa, x = "Island", y = "RPM", fill="Island", add=c("boxplot"),add.params = c(list(fill = "white"), list(width=0.05))) + facet_wrap(~Pathogen, scales = "free")
dev.off()

# back to Phyloseq pipeline -----------------------------

CCMeta_physeq = phyloseq(taxa, tax)
p = plot_bar(CCMeta_physeq, fill = "Phylum")

fig <- p + geom_bar(aes(fill=Phylum), stat="identity", position="stack") +
	guides(fill=guide_legend(ncol=2))

pdf(paste0(outputdir,"Allpathogens_Phylum.pdf"), width=15)
fig
dev.off()

# plot trends for each taxa
melted_taxa=melt(taxa)
colnames(melted_taxa)=c("Pathogen","Island","RPM")
melted_taxa$Island=gsub("\\..*","",melted_taxa$Island)
melted_taxa$Pathogen=gsub("Bacteria_","",melted_taxa$Pathogen) %>% gsub("Eukaryota_","",.) %>% gsub("Viruses_","",.)
pdf(paste0(outputdir,"pathogensByIsland_AllPathogens_Phylum_Boxplot.pdf"))
ggboxplot(melted_taxa, x = "Island", y = "RPM", fill="Island", add=c("boxplot"),add.params = c(list(fill = "white"), list(width=0.05))) + facet_wrap(~Pathogen, scales = "free")
dev.off()

# without Plasmodium
taxa=taxa[-which(rownames(taxa) %in% "Eukaryota_Apicomplexa"),]

CCMeta_physeq = phyloseq(taxa, tax)
p=plot_bar(CCMeta_physeq, fill = "Phylum")

fig <- p + geom_bar(aes(fill=Phylum), stat="identity", position="stack") +
  guides(fill=guide_legend(ncol=2))

pdf(paste0(outputdir,"Allpathogens_NoPlasmo_Phylum.pdf"), width=15)
fig
dev.off()
