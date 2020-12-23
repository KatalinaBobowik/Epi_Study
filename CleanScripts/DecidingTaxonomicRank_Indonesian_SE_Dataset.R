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
library(stats)
library(ggrepel)

# set ggplot colour theme to white
theme_set(theme_bw())

# set up directories
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/SEIndonesianSamples/"
outputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/Indo_250K/"

# load in single-ended Indonesian phyloseq counts object with all reads included (not subsampled)
load(paste0(inputdir,"AllREadsSE_Indo_Counts_physeq_filtered.Rda"))


# Data preprocessing ------------------------------------

# look at the rate of missingness and number of reads at each rank
Species=tax_glom(AllREadsSE_Indo_Counts_physeq,"Species")
Genus=tax_glom(AllREadsSE_Indo_Counts_physeq,"Genus")
Family=tax_glom(AllREadsSE_Indo_Counts_physeq,"Family")
Order=tax_glom(AllREadsSE_Indo_Counts_physeq,"Order")
Class=tax_glom(AllREadsSE_Indo_Counts_physeq,"Class")
Phylum=tax_glom(AllREadsSE_Indo_Counts_physeq,"Phylum")
Kingdom=tax_glom(AllREadsSE_Indo_Counts_physeq,"Kingdom")
Superkingdom=tax_glom(AllREadsSE_Indo_Counts_physeq,"Superkingdom")

# get number of reads for each species
Species_reads=data.frame(colSums(otu_table(Species)),"Species") %>% `colnames<-`(c("Reads","Taxonomic_Rank"))
Genus_reads=data.frame(colSums(otu_table(Genus)),"Genus") %>% `colnames<-`(c("Reads","Taxonomic_Rank"))
Family_reads=data.frame(colSums(otu_table(Family)),"Family") %>% `colnames<-`(c("Reads","Taxonomic_Rank"))
Order_reads=data.frame(colSums(otu_table(Order)),"Order") %>% `colnames<-`(c("Reads","Taxonomic_Rank"))
Class_reads=data.frame(colSums(otu_table(Class)),"Class") %>% `colnames<-`(c("Reads","Taxonomic_Rank"))
Phylum_reads=data.frame(colSums(otu_table(Phylum)),"Phylum") %>% `colnames<-`(c("Reads","Taxonomic_Rank"))
Kingdom_reads=data.frame(colSums(otu_table(Kingdom)),"Kingdom") %>% `colnames<-`(c("Reads","Taxonomic_Rank"))
Superkingdom_reads=data.frame(colSums(otu_table(Superkingdom)),"Superkingdom") %>% `colnames<-`(c("Reads","Taxonomic_Rank"))

# make a table with the number of reads at each taxonomic rank for each sample
missingness_summary=data.frame(Species_reads$Reads,Genus_reads$Reads,Family_reads$Reads,Order_reads$Reads,Class_reads$Reads,Phylum_reads$Reads,Kingdom_reads$Reads,Superkingdom_reads$Reads)
colnames(missingness_summary)=c("Species","Genus","Family","Order","Class","Phylum","Kingdom","Superkingdom")
rownames(missingness_summary)=rownames(Species_reads)
# save missingness table
write.table(missingness_summary,file=paste0(outputdir,"Missingness_Summary_PerSample.txt"),sep="\t")

# first plot number of missing reads
missingness_reads=data.frame(Species=c(1:ncol(otu_table(AllREadsSE_Indo_Counts_physeq))), Genus=c(1:ncol(otu_table(AllREadsSE_Indo_Counts_physeq))),Family=c(1:ncol(otu_table(AllREadsSE_Indo_Counts_physeq))),Order=c(1:ncol(otu_table(AllREadsSE_Indo_Counts_physeq))),Class=c(1:ncol(otu_table(AllREadsSE_Indo_Counts_physeq))),Phylum=c(1:ncol(otu_table(AllREadsSE_Indo_Counts_physeq))),Kingdom=c(1:ncol(otu_table(AllREadsSE_Indo_Counts_physeq))),Superkingdom=c(1:ncol(otu_table(AllREadsSE_Indo_Counts_physeq))))
for (column.name in names(missingness_summary)) {
  missingness_reads[column.name] = (missingness_summary$Superkingdom - missingness_summary[column.name])[,1]
}
rownames(missingness_reads)=rownames(Species_reads)

# remove samples which are 0 throughout
missingness_reads=missingness_reads[complete.cases(missingness_reads), ]
missingness_reads=missingness_reads[rowSums(missingness_reads > 0) != 0, ]
missingness_reads$Sample=rownames(missingness_reads)
missingness_reads_melt=melt(missingness_reads)

pdf(paste0(outputdir,"MissingnessSummary_nReads_PerSample.pdf"), width=10)
ggplot(missingness_reads_melt, aes(fill=variable, y=value, x=Sample)) + 
    geom_bar(position="stack", stat="identity") +
    coord_flip() +
    labs(fill = "Taxonomic Rank") +
    ylab("n Reads") + ggtitle("Sample Missingness")
dev.off()

# now make this into a plot
for (column.name in names(missingness_summary)) {
  missingness_summary[column.name] = 1 - (missingness_summary[column.name] / missingness_summary$Superkingdom)
}
missingness_summary = missingness_summary * 100

# remove NA values and samples which are 0 throughout
missingness_summary=missingness_summary[complete.cases(missingness_summary), ]
missingness_summary=missingness_summary[rowSums(missingness_summary > 0) != 0, ]
missingness_summary$Sample=rownames(missingness_summary)
missingness_summary_melt=melt(missingness_summary)

pdf(paste0(outputdir,"MissingnessSummary_PerSample.pdf"), width=10)
ggplot(missingness_summary_melt, aes(fill=variable, y=value, x=Sample)) + 
    geom_bar(position="fill", stat="identity") +
    coord_flip() +
    labs(fill = "Taxonomic Rank") +
    ylab("% Missingness") + ggtitle("Sample Percent Missingness")
dev.off()

# Percent missingness in each taxonomic rank
total_OTUs=nrow(tax_table(AllREadsSE_Indo_Counts_physeq))
# what is the maximum number of OTUs?
total_OTUs
# 1386

# get total OTUs for each species
taxonomicRanks=c("Species","Genus","Family","Order","Class","Phylum","Kingdom","Superkingdom")
missingness_vector=c()
for (rank in taxonomicRanks){
	missingness=sum(tax_table(AllREadsSE_Indo_Counts_physeq)[,rank]=="")/total_OTUs
	missingness_vector=c(missingness_vector, missingness)
}

# merge total OTU information into a df
SpeciesMissingnessDF=data.frame(species=taxonomicRanks, missingness=missingness_vector)
SpeciesMissingnessDF$notMissing=1-SpeciesMissingnessDF$missingness
melt=melt(SpeciesMissingnessDF)
melt$species <- factor(melt$species, levels = taxonomicRanks)

# plot
pdf(paste0(outputdir,"SpeciesPercentMissingness.pdf"), width=10)
ggplot(melt, aes(fill=forcats::fct_rev(variable), y=value, x=species)) + 
    geom_bar(position="fill", stat="identity") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.title = element_blank()) + 
    scale_fill_manual(values=c("grey","blue"), labels = c("Not Missing", "Missing")) +
    ylab("% Missingness") + ggtitle("Species Percent Missingness")
dev.off()

# now get number of reads, add to the df, and plot
nReads_vector=c()
for (rank in taxonomicRanks){
	nReads=sum(get(paste0(rank,"_reads"))$Reads)
	nReads_vector=c(nReads_vector, nReads)
}
SpeciesMissingnessDF$reads=nReads_vector
# remove not missing column
SpeciesMissingnessDF=SpeciesMissingnessDF[,!(names(SpeciesMissingnessDF) %in% "notMissing")]
# save table
write.table(SpeciesMissingnessDF,file=paste0(outputdir,"SpeciesMissingness_ReadsAndPercMissingness.txt"),sep="\t")

# Finally, plot this information
pdf(paste0(outputdir,"SpeciesPercentMissingness.pdf"), width=10)
ggplot(SpeciesMissingnessDF, aes(x=missingness, y=reads, group=species,label=species)) +
  geom_point(aes(color=species,alpha=0.3),size = 6,position = position_jitter(w = 0.005, h = 1)) +
  theme(legend.position = "none") +
  geom_text_repel()
dev.off()