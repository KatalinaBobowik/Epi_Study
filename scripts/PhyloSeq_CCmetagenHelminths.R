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
outputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/"

# Import data and convert to a phyloseq object
# raw_CCMetagen_data <- read.csv(paste0(inputdir,"IndoUnmapped_species_table_Helminths.csv"),check.names=FALSE)
raw_CCMetagen_data <- read.csv(paste0(inputdir,"IndoUnmapped_species_table_noRepeats_RPM.csv"),check.names=FALSE)

# remove plants form the analysis
raw_CCMetagen_data=raw_CCMetagen_data[-which(raw_CCMetagen_data$Kingdom=="Viridiplantae"),]
# separate bacterial from eukaryotic unclassified taxa
raw_CCMetagen_data$SuperKFamily <- paste(raw_CCMetagen_data$Superkingdom, raw_CCMetagen_data$Family, sep="_")
# add an 'unclassified' string at the end of these to avoid confusion
raw_CCMetagen_data <- data.frame(lapply(raw_CCMetagen_data, function(x) {sub("Bacteria_$", "Bacteria_unclassified", x)}))
raw_CCMetagen_data <- data.frame(lapply(raw_CCMetagen_data, function(x) {sub("Eukaryota_$", "Eukaryota_unclassified", x)}))
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
# get first instance of unclassified eukaryotes and merge with other unclassified eukaryotes 
CCMetagen_data[grep("Eukaryota_unk_f|Eukaryota_unclassified",CCMetagen_data$Family)[1],4:128]=colSums(CCMetagen_data[grep("Eukaryota_unk_f|Eukaryota_unclassified",CCMetagen_data$Family),4:128])
# delete the rest
CCMetagen_data=CCMetagen_data[-grep("Eukaryota_unk_f|Eukaryota_unclassified",CCMetagen_data$Family)[2:length(grep("Eukaryota_unk_f|Eukaryota_unclassified",CCMetagen_data$Family))],]
CCMetagen_data=droplevels(CCMetagen_data)

taxa_raw <- as.matrix(CCMetagen_data[,c("Superkingdom","Kingdom","Family")])
rownames(taxa_raw) <- taxa_raw[,"Family"]
abund_raw <- as.matrix(CCMetagen_data[,-which(colnames(CCMetagen_data) %in% c("Superkingdom","Kingdom","Family","Genus","Species"))])
rownames(abund_raw) <- CCMetagen_data[,which(names(CCMetagen_data)=="Family")]

# convert to Phyloseq object
tax = tax_table(taxa_raw)
taxa = otu_table(abund_raw, taxa_are_rows = TRUE)

# exploratory stuff on species per island
df=melt(taxa)
df$island=NA
df$island[grep("MPI",df$Var2)]="MPI"
df$island[grep("SMB",df$Var2)]="SMB"
df$island[grep("MTW",df$Var2)]="MTW"

pdf(paste0(outputdir,"pathogensByIsland_Helminths.pdf"), width=20)
ggplot(df, aes(x=island, y=value, fill=Var1)) + geom_bar(width = 1, stat = "identity") + labs(y="RPM", x = "Island")
dev.off()

# the strong Plasmo signal is driven by one sample, SMB-PTB-028. Let's take him out
taxa_noSMB.PTB028=taxa[,-which(colnames(taxa) %in% "SMB.PTB028_noRNA")]
df=melt(taxa_noSMB.PTB028)
df$island=NA
df$island[grep("MPI",df$Var2)]="MPI"
df$island[grep("SMB",df$Var2)]="SMB"
df$island[grep("MTW",df$Var2)]="MTW"

# replot
pdf(paste0(outputdir,"pathogensByIsland_NoSMBPTB028_Helminths.pdf"), width=20)
ggplot(df, aes(x=island, y=value, fill=Var1)) + geom_bar(width = 1, stat = "identity")
dev.off()

df=df[order(df$value,decreasing=T),]
topGenes_Plasmo <- by(df, df["island"], head, n=20)
top=by(df, df["island"], head, n=20)
top=rbind(top[["SMB"]],top[["MPI"]],top[["MTW"]])
pdf(paste0(outputdir,"pathogensByIsland_NoSMBPTB028_Top20Helminths.pdf"))
ggplot(top, aes(x=island, y=value, fill=Var1)) + geom_bar(width = 1, stat = "identity")
dev.off()

# back to CCmetagen pipeline
CCMeta_physeq = phyloseq(taxa, tax)
plot_bar(CCMeta_physeq, fill = "Superkingdom")

TopNOTUs <- names(sort(taxa_sums(CCMeta_physeq), TRUE)[1:16])
TopFamilies <- prune_taxa(TopNOTUs, CCMeta_physeq)
plot_bar(TopFamilies, fill = "Family")

p = plot_bar(TopFamilies, fill="Family")
p$data$Sample <- factor(p$data$Sample)

familiesMPI = levels(p$data$Family)
p$data$Family <- factor(p$data$Family) 

PaletteBacteria = colorRampPalette(c("#ebe302", "#eb7a02"))(13)
PaletteEukaryote = c("#1c9c02")
PaletteVirus = c("#7b029c","#9c0273")

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

# without Plasmodium
taxa=taxa[-which(rownames(taxa) %in% "Eukaryota_Plasmodiidae"),]

CCMeta_physeq = phyloseq(taxa, tax)
plot_bar(CCMeta_physeq, fill = "Superkingdom")

TopNOTUs <- names(sort(taxa_sums(CCMeta_physeq), TRUE)[1:16])
TopFamilies <- prune_taxa(TopNOTUs, CCMeta_physeq)
plot_bar(TopFamilies, fill = "Family")

p = plot_bar(TopFamilies, fill="Family")
p$data$Sample <- factor(p$data$Sample)

familiesMPI = levels(p$data$Family)
p$data$Family <- factor(p$data$Family) 

PaletteBacteria = colorRampPalette(c("#ebe302", "#eb7a02"))(13)
PaletteEukaryote = c("#1c9c02")
PaletteVirus = c("#7b029c","#9c0273")

Merged_Palette <- c(PaletteBacteria,PaletteEukaryote,PaletteVirus)

fig <- p + scale_fill_manual(values=Merged_Palette) +
  geom_bar(aes(fill=Family), stat="identity", position="stack") +
  guides(fill=guide_legend(ncol=2))

pdf(paste0(outputdir,"Allpathogens_NoHelminthsOrPlasmo.pdf"), width=15)
fig
dev.off()
