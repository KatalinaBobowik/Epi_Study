# created by KSB, 20.01.19
# Script created to collate files from reads mapping to other organisms on Spartan server
# this can be saved and then used locally

# set working directory
setwd("/data/cephfs/punim0586/kbobowik/Diamond/matches_output")


# find all files that match accession number + organism + accession count frequency file, then turn into a list
# on Spartan server
temp = list.files(pattern="matches_oneHit")
myfiles = lapply(temp, read.table, sep="\t")

# assign sample names
names(myfiles)=paste(sapply(strsplit(temp,"[_.]"), `[`, 5), sapply(strsplit(temp,"[_.]"), `[`, 3), sep="_")
# assign column names
myfiles=lapply(myfiles, setNames, nm = c("Accession_Number", "Organism", "Accession_Frequency"))
# sort by frequency of hits
myfiles=lapply(myfiles, function(x) x[order(-x$Accession_Frequency),])
# Get the name of the genus
myfiles=lapply(myfiles, function(x) cbind(x, sapply(strsplit(as.character(x$Organism),"[ ]"), `[`, 1)))
# rename last column to 'genus'
myfiles=lapply(myfiles, setNames, nm = c("Accession_Number", "Organism", "Accession_Frequency", "Genus"))

# save as rds file and work with this locally
saveRDS(myfiles, file = "UnmappedReads_Metadata_OneHit.rds")