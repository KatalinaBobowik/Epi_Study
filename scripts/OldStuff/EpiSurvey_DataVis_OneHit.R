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

speciesInfo=read.table(paste0(dir,"/data/speciesInfo_OneHit.txt"),stringsAsFactors=F)

# check dimensions of dataframe
dim(speciesInfo)
# [1] 697   7

# since 'superkingdom' has the least amount of missing information, we'll use this column to check if there are any NA values
speciesInfo[which(is.na(speciesInfo$superkingdom)),]$Organism

 #[1] "synthetic construct"                          
 #[2] "Plasmodium falciparum Tanzania (2000708)"     
 #[3] "Plasmodium falciparum Vietnam Oak-Knoll (FVO)"
 #[4] "Mn"                                           
 #[5] "Cloning vector pInSRT-GFPV1"                  
 #[6] "trans-4-(4-methylpiperazin-1-yl)cyclohexyl"   
 #[7] "Entamoeba histolytica HM-1:IMSS"              
 #[8] "ADP-ribose"                                   
 #[9] "NADP(+)"                                      
#[10] "GB_virus"                                     
#[11] "Knock-in vector FBDki"                        
#[12] "unidentified cloning vector"                  
#[13] "Synthetic plasmid pRK96"                      
#[14] "glutamine-hydrolyzing"  

# let's get rid of the entries where nothing has been found and are not interesting to us
unwantedOrganisms=c(1,4:6,8,9,11:14)
speciesInfo[which(is.na(speciesInfo$superkingdom)),][unwantedOrganisms,]$Organism

# [1] "synthetic construct"                       
# [2] "Mn"                                        
# [3] "Cloning vector pInSRT-GFPV1"               
# [4] "trans-4-(4-methylpiperazin-1-yl)cyclohexyl"
# [5] "ADP-ribose"                                
# [6] "NADP(+)"                                   
# [7] "Knock-in vector FBDki"                     
# [8] "unidentified cloning vector"               
# [9] "Synthetic plasmid pRK96"                   
# [10] "glutamine-hydrolyzing"

# get rid of unwanted organisms from df
speciesInfo=speciesInfo[-which(is.na(speciesInfo$superkingdom))[unwantedOrganisms],]

# The organisms that are coming up as NA do not have proper names recognised by NCBI. Let's rename them
missingOrgReplacement=c("Plasmodium falciparum","Plasmodium falciparum","Entamoeba histolytica","GB virus C")
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
# [1] 13

# for all entries with 'NA' values, we can fill in information from the family and phylum columns
# first fill in information from family
speciesInfo[which(is.na(speciesInfo$genus)),"genus"]=speciesInfo[which(is.na(speciesInfo$genus)),"family"]
# now fill in infomration from phylum
speciesInfo[which(is.na(speciesInfo$genus)),"genus"]=speciesInfo[which(is.na(speciesInfo$genus)),"phylum"]
# finally, we can fill in the last field with he superkingdom
speciesInfo[which(is.na(speciesInfo$genus)),"genus"]=speciesInfo[which(is.na(speciesInfo$genus)),"superkingdom"]

# add in all species information by merging data together by organism name
org_freq=lapply(org_freq, function(x) merge(speciesInfo,x, by="Organism"))
# remove phylum that we're not interested in
org_freq=lapply(org_freq, function(x) x[grep("Chordata|Arthropoda|Mollusca|Cnidaria|Streptophyta|Echinodermata", x$phylum, invert=T),])
# We filtered out organisms we're not interested in by phylum, however there are still some NAs in phylum and we may have missed some organisms to filer out
# we can get all the 'NA' entries in phylum and get the corresponding family name
unique(melt(lapply(org_freq, function(x) unique(x[which(is.na(x$phylum)),"family"])))[,1])
# [1] <NA>            Reticulomyxidae Oxytrichidae    Flaviviridae   
# [5] Peronosporaceae
# Levels: Oxytrichidae Reticulomyxidae Flaviviridae Peronosporaceae 

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
colours[c(4,8,19,20,23)]=4
colours[c(16)]=3
colours[c(17)]=2

pdf("readsMappingToPathogens_oneHit.pdf", height=10, width=20)
ggplot(genus.df,aes(x=factor(Sample_Name), y=Accession_Frequency)) + geom_col(aes(fill=factor(genus.df$Genus))) + scale_fill_manual(values=colours, name="Organism") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + labs(x="Sample Name")
dev.off()