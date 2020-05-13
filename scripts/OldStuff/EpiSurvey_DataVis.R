# created by KSB, 20.01.19
# Script created to inspect where reads are mapping to. 


# load packages
require(ggplot2)
require(reshape2)
library(taxize)
library(ape)

# set wd where output will go to
setwd("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epidemiological_Survey")

# define GIT directory
dir="/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/GIT/Epidemiological_Survey"

# plot reads without any filtering --------------------------------

# read in rds file. This file was created on the server from the script 'collateFiles_EpiSurvey.R'
unmappedReads=readRDS(file = paste0(dir,"/data/UnmappedReads_Metadata.rds"))

# merge duplicated organism together and get counts
org_freq=lapply(unmappedReads, function(x) aggregate(Accession_Frequency ~ Organism, data=x, FUN=sum))
# sort by accession frequency order (descending)
org_freq=lapply(org_freq, function(x) x[order(-x[,2]),])

# show the top hits for the first ten most abundant organisms
topTenUnfiltered=lapply(org_freq, function(x) x[1:10,])
topTenUnfiltered=melt(topTenUnfiltered)

# now plot
pdf("topTenUNfilteredGenes.pdf", height=15, width=15)
ggplot(topTenUnfiltered,aes(x=factor(L1), y=value)) + geom_col(aes(fill=factor(topTenUnfiltered$Organism))) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()

# Organise species information ------------------------------

speciesInfo=read.table(paste0(dir,"/data/speciesInfo.txt"),stringsAsFactors=F)

# since 'superkingdom' has the least amount of missing information, we'll use this column to check if there are any NA values
speciesInfo[which(is.na(speciesInfo$superkingdom)),]$Organism

# [1] "synthetic construct"                                        
# [2] "Plasmodium falciparum Vietnam Oak-Knoll (FVO)"              
# [3] "3-Indolyl"                                                  
# [4] "human, urine, chronic renal failure patient, Peptide, 98 aa"
# [5] "GB_virus"                                                   
# [6] "Synthetic plasmid pRK96"                                    
# [7] "ADP-ribose"                                                 
# [8] "alpha V96w"                                                 
# [9] "human, B-CLC cells, Peptide Partial, 182 aa"                
# [10] "Entamoeba histolytica HM-1:IMSS"                            
# [11] "Plasmodium falciparum Tanzania (2000708)"                   
# [12] "hydroxy(Methoxy)methyl"                                     
# [13] "Expression vector pInSRT-hM(V1)"                            
# [14] "Cloning vector pInSRT-GFPV1" 

# let's get rid of the entries where nothing has been found and are not interesting to us
unwantedOrganisms=c(1,3,4,6:9,12,13,14)
speciesInfo[which(is.na(speciesInfo$superkingdom)),][unwantedOrganisms,]$Organism

# [1] "synthetic construct"                                        
# [2] "3-Indolyl"                                                  
# [3] "human, urine, chronic renal failure patient, Peptide, 98 aa"
# [4] "Synthetic plasmid pRK96"                                    
# [5] "ADP-ribose"                                                 
# [6] "alpha V96w"                                                 
# [7] "human, B-CLC cells, Peptide Partial, 182 aa"                
# [8] "hydroxy(Methoxy)methyl"                                     
# [9] "Expression vector pInSRT-hM(V1)"                            
# [10] "Cloning vector pInSRT-GFPV1"

# get rid of unwanted organisms from df
speciesInfo=speciesInfo[-which(is.na(speciesInfo$superkingdom))[unwantedOrganisms],]

# The organisms that are coming up as NA do not have proper names recognised by NCBI. Let's rename them
missingOrgReplacement=c("Plasmodium falciparum","GB virus C","Entamoeba histolytica","Plasmodium falciparum")
speciesInfo[which(is.na(speciesInfo$superkingdom)),]$Organism=missingOrgReplacement

# make a new index and new spcies names after getting rid of entries
index=which(is.na(speciesInfo$superkingdom))
newSpeciesNames=as.character(speciesInfo$Organism[index])

counter=0
for (i in index){
	counter=counter+1
	speciesInfo[index[counter],]=tax_name(query = newSpeciesNames[counter], get = c("superkingdom","phylum", "order", "family", "genus"), db = "ncbi")
}

# check if there are any NA entries in superkingdom
speciesInfo[which(is.na(speciesInfo$superkingdom)),]
#[1] db           Organism     superkingdom phylum       order       
#[6] family       genus       
#<0 rows> (or 0-length row.names)

# since we'll be merging data by genus, we want to make sure there are no 'NA' entries. Let's see how many 'NA' entries we have for genus
nrow(speciesInfo[which(is.na(speciesInfo$genus)),])
# [1] 8

# for all entries with 'NA' values, we can fill in information from the family and phylum columns
# first fill in information from family
speciesInfo[which(is.na(speciesInfo$genus)),"genus"]=speciesInfo[which(is.na(speciesInfo$genus)),"family"]
# now fill in infomration from phylum
speciesInfo[which(is.na(speciesInfo$genus)),"genus"]=speciesInfo[which(is.na(speciesInfo$genus)),"phylum"]
# let's check one last time to see if there are any empty entries in genus
speciesInfo[which(is.na(speciesInfo$genus)),]
#       db            Organism superkingdom phylum order family genus
# 535 ncbi Terrabacteria group     Bacteria   <NA>  <NA>   <NA>  <NA>

# This entry wasn't picked up so we'll fill in the phylum,family, and genus information with 'Bacteria'
speciesInfo[which(is.na(speciesInfo$genus)),c("phylum","family","genus")]="Bacteria"

# add in all species information by merging data together by organism name
org_freq=lapply(org_freq, function(x) merge(speciesInfo,x, by="Organism"))
# remove phylum that we're not interested in
org_freq=lapply(org_freq, function(x) x[grep("Chordata|Arthropoda|Mollusca|Cnidaria|Streptophyta|Echinodermata", x$phylum, invert=T),])
# We filtered out organisms we're not interested in by phylum, however there are still some NAs in phylum and we may have missed some organisms to filer out
# we can get all the 'NA' entries in phylum and get the corresponding family name
unique(melt(lapply(org_freq, function(x) unique(x[which(is.na(x$phylum)),"family"])))[,1])
#[1] Flaviviridae       Reticulomyxidae    Oxytrichidae       Peronosporaceae   
#[5] Holostichidae      Sphaerozoidae      Pseudokeronopsidae Salpingoecidae    
#[9] Noelaerhabdaceae  

# Noelaerhabdaceae is an algae, so we'll remove it
org_freq=lapply(org_freq, function(x) x[grep("Noelaerhabdaceae", x$family, invert=T),])

# Plot filtered data --------------------------------

# re-order by accession frequency
z=lapply(org_freq, function(x) x[order(-x[,"Accession_Frequency"]),])
# just get the top ten entries ordered by frequency of accession number
z=lapply(z, function(x) x[1:10,])
# melt data
z=melt(z)
genus.df=data.frame(z$genus,z$value,z$L1)
colnames(genus.df)=c("Genus","Accession_Frequency","Sample_Name")

# set up colours
colours=rep(1,length(unique(genus.df$Genus)))
colours[c(2,6,14,17)]=2
colours[c(11,20,23:25)]=4
colours[c(15:16)]=3

colours2=rep(1,length(unique(genus.df$Genus)))
colours2[c(12,13)]=4
colours2[17]=2

pdf("readsMappingToPathogens_within1Percent.pdf", height=10, width=20)
ggplot(genus.df,aes(x=factor(Sample_Name), y=Accession_Frequency)) + geom_col(aes(fill=factor(genus.df$Genus))) + scale_fill_manual(values=colours) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
ggplot(genus.df,aes(x=factor(Sample_Name), y=Accession_Frequency)) + geom_col(aes(fill=factor(genus.df$Genus))) + scale_fill_manual(values=colours2) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()