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

# read in rds file. This file was created on the server from the script 'collateFiles_EpiSurvey.R'
unmappedReads=readRDS(file = paste0(dir,"/data/UnmappedReads_Metadata.rds"))

# merge duplicated organism together and get counts
org_freq=lapply(unmappedReads, function(x) aggregate(Accession_Frequency ~ Organism, data=x, FUN=sum))
# sort by accession frequency order (descending)
org_freq=lapply(org_freq, function(x) x[order(-x[,2]),])

# Get information for species --------------

# first, melt the list of all samples and merge by organism name. Doing this will remove redundant organism names so getting the phylum information won't take as long
allOrganisms=melt(org_freq,id.vars=c('Accession_Frequency','Organism'))
# Get species information from NCBI database with 'taxize' package
species=unique(as.character(allOrganisms[,2]))
# This takes a while- around 20 minutes. Time to get a cup of coffee :)
speciesInfo=tax_name(query = c(species), get = c("superkingdom","phylum", "order", "family", "genus"), db = "ncbi")
# rename column "query" to "Organism"
colnames(speciesInfo)[2]="Organism"

# some information is still missing in the speciesInfo dataframe. Most of these empty entries are due to the accession number being in the organism column. We can get all of the 'NA' entries in the speciesInfo dataframe and try and fill it in with organism info from the 'Ape' package
index=which(is.na(speciesInfo$genus) & is.na(speciesInfo$phylum))
na.speciesInfo=speciesInfo$Organism[index]

counter=0
for(entry in na.speciesInfo){
	counter=counter+1
	print(paste(counter,"out of",length(index), sep=" "))
	possibleError <- tryCatch(
      read.GenBank(entry),
      error=function(e) e
  )
	# if there is an error, go to the next entry
	if(inherits(possibleError, "error")){
		next # speciesInfo[index[counter],"phylum"]=NA
	}
	# otherwise, rename the organism name to the species name found by the 'ape' package
	else {
		speciesInfo$Organism[index[counter]]=attr(read.GenBank(entry), "species")
	}
}

# get information for species which haven't been picked up by NCBI by getting names from the organism field which was incorporated in the previous step
newSpeciesNames=speciesInfo$Organism[index]
counter=0
for (i in index){
	counter=counter+1
	print(paste(counter,"out of",length(index), sep=" "))
	possibleError <- tryCatch(
      read.GenBank(entry),
      error=function(e) e
  )
	# if there is an error, go to the next entry
	if(inherits(possibleError, "error")){
		next # speciesInfo[index[counter],"phylum"]=NA
	}
	# otherwise, rename the organism name to the species name found by the 'ape' package
	else {
		speciesInfo[index[counter],]=tax_name(query = newSpeciesNames[counter], get = c("superkingdom","phylum", "order", "family", "genus"), db = "ncbi")
	}
}

# save table
write.table(speciesInfo,file = paste0(dir,"/data/speciesInfo.txt"))
