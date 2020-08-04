# created by KSB, 20.01.19
# Script created to inspect where reads are mapping to. 

# load packages
require(ggplot2)
require(RColorBrewer)
library(dplyr)
library(plyr)
library(phyloseq)
library(reshape2)
library(tidyverse)

# set ggplot colour theme to white
theme_set(theme_bw())

# set up directories
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/EpiStudy/"
outputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/"

# read in file and assigtn column names
summaryAccessionTable=read.delim(paste0(inputdir,"summaryAccessionTable_Filtered.txt"),header=F)
colnames(summaryAccessionTable)=c("kingdom","kingdom_instance","phylum","phylum_instance","class","class_instance","order","order_instance","family","family_instance","genus","genus_instance","species","species_instance")

# remove kingdom and kingdom instance, as they are too broad
summaryAccessionTable_noKingdom=subset(summaryAccessionTable, select=-c(kingdom,kingdom_instance))
# add in pphylum column and melt df
summaryAccessionTable_noKingdom$phylum=summaryAccessionTable$phylum
summaryAccessionTable_noKingdom=melt(summaryAccessionTable_noKingdom)

# visualise data
ggplot(data=summaryAccessionTable_noKingdom, aes(x=variable, y=value, group=phylum)) + 
	geom_line(aes(color=phylum))+ 
	geom_point(aes(color=phylum))

# apicomplexa seems to be very high at most taxonomic ranks, so remove it and replot
summaryAccessionTable_noKingdom=summaryAccessionTable_noKingdom %>% filter(str_detect(phylum, "Apicomplexa", negate = TRUE))
ggplot(data=d, aes(x=variable, y=value, group=phylum)) + 
	geom_line(aes(color=phylum)) + 
	geom_point(aes(color=phylum))

# unk_p seems to be very high at phylum, class, and order, so remove it and replot
summaryAccessionTable_noKingdom=summaryAccessionTable_noKingdom %>% filter(str_detect(phylum, "unk_p", negate = TRUE))
ggplot(data=d, aes(x=variable, y=value, group=phylum)) + 
	geom_line(aes(color=phylum)) + 
	geom_point(aes(color=phylum))

# Ascomycota seems ot be very high only at the order level, whereas other taxonomic ranks are still quite low
# This is driven by 'unk_o' in order. We'll leave it for now to keep investigating other signals but Ascomycota should be removed
summaryAccessionTable_noKingdom=summaryAccessionTable_noKingdom %>% filter(str_detect(phylum, "Ascomycota", negate = TRUE))
ggplot(data=d, aes(x=variable, y=value, group=phylum)) + 
	geom_line(aes(color=phylum)) +
	geom_point(aes(color=phylum))

# Proteobacteria is high at phylum, class, order, family, and genus, so remove it
summaryAccessionTable_noKingdom=summaryAccessionTable_noKingdom %>% filter(str_detect(phylum, "Proteobacteria", negate = TRUE))
ggplot(data=d, aes(x=variable, y=value, group=phylum)) + 
	geom_line(aes(color=phylum)) +
	geom_point(aes(color=phylum))

pdf(paste0(outputdir,"~/Desktop/phylum_filtered.pdf"))
for (name in unique(summaryAccessionTable_noKingdom$phylum)){
	phylum_filtered=summaryAccessionTable_noKingdom %>% filter(str_detect(phylum, name, negate = FALSE))
	print(ggplot(data=phylum_filtered, aes(x=variable, y=value, group=phylum)) + 
		geom_line(aes(color=phylum)) +
		geom_point(aes(color=phylum)))
}
dev.off()

# We can confidently filter these out:
unique(summaryAccessionTable_noKingdom$phylum)
#  [1] Actinobacteria  Bacteroidetes   Firmicutes      Basidiomycota  
#  [5] Cnidaria        Echinodermata   Nematoda        Platyhelminthes
#  [9] Priapulida      Bacillariophyta Euglenozoa      Evosea         
# [13] Haptista 

write.table(unique(summaryAccessionTable_noKingdom$phylum), file=paste0(outputdir,"lowAccessiontaxa.txt"))


