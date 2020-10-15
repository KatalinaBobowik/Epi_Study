# created by KSB, 20.01.19
# Script created to inspect where reads are mapping to. 

# load packages
require(ggplot2)
require(RColorBrewer)
library(dplyr)
library(plyr)
library(reshape2)
library(ggpubr)
library(vegan)
library(ggfortify)
library(microbiome)
library(microbiomeutilities)
library(viridis)
library(tibble)
library(knitr)

# set ggplot colour theme to white
theme_set(theme_bw())

# set up directories
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/All_Reads/"
outputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/All_Reads/"

# load in allReads Indonesian phyloseq counts object
load(paste0(inputdir,"allReads_Counts_physeq.Rda"))
# load in allreads Indonesian phyloseq RPM object
load(paste0(inputdir,"allReads_RPM_physeq.Rda"))

# Remove Metazoa and Viridiplantae
allReads_Counts_physeq=subset_taxa(allReads_Counts_physeq, (Kingdom!="Viridiplantae"))
allReads_Counts_physeq=subset_taxa(allReads_Counts_physeq, (Kingdom!="Metazoa"))
# Remove Metazoa and Viridiplantae for RPM
allReads_RPM_physeq=subset_taxa(allReads_RPM_physeq, (Kingdom!="Viridiplantae"))
allReads_RPM_physeq=subset_taxa(allReads_RPM_physeq, (Kingdom!="Metazoa"))

# filter reads 1 and below
otu_table(allReads_Counts_physeq)[otu_table(allReads_Counts_physeq)<=1]<-0
otu_table(allReads_RPM_physeq)[otu_table(allReads_RPM_physeq)<=1]<-0

# Data preprocessing ------------------------------------

# Check distribution of how many reads/samples?
SeqDepth = colSums(otu_table(allReads_Counts_physeq))
sample_data(allReads_Counts_physeq)$SeqDepth = SeqDepth

# make a histogran of sequencing depth
pdf(paste0(outputdir,"HistogramSeqDepth_AllReadsIndonesian.pdf"), width=15)
ggplot(meta(allReads_Counts_physeq)) +
    geom_histogram(aes(x = SeqDepth), alpha= 0.6)
dev.off()

# get the minimum and maximum sequencing depth
min(SeqDepth)
# 0
max(SeqDepth)
# 23403

table(colSums(otu_table(allReads_Counts_physeq)))
  #   0     2     3     4     6     7     8     9    10    11    12    13    14 
  #  69    12     2     4     1     2     2     2     3     1     1     1     2 
 
  #  17    18    22    26    27    34    35    51    58   171   209   218   247 
  #   1     3     1     2     1     1     1     1     1     1     1     1     1 
 
  # 726  1632  3790 14307 23403 
  #   1     1     1     1     1 

 # You can see that 69 samples have 0 values, and 54 have values

# which samples actually have data?
names(which(colSums(otu_table(allReads_Counts_physeq))>0))
#  [1] "MPI-025"           "MPI-061"           "MPI-065"          
#  [4] "MPI-296"           "MPI-333"           "MPI-334"          
#  [7] "MPI-335"           "MPI-345"           "MPI-376"          
# [10] "MPI-389"           "MTW-MDB-003"       "MTW-MDB-009"      
# [13] "MTW-MDB-020"       "MTW-MDB-024"       "MTW-MDB-034"      
# [16] "MTW-MDB-035"       "MTW-TLL-002"       "MTW-TLL-006"      
# [19] "MTW-TLL-008"       "MTW-TLL-010"       "MTW-TLL-011"      
# [22] "MTW-TLL-013Batch2" "MTW-TLL-017"       "MTW-TLL-024"      
# [25] "MTW-TLL-027"       "MTW-TLL-029"       "MTW-TLL-030"      
# [28] "MTW-TLL-032"       "MTW-TLL-034"       "MTW-TLL-035"      
# [31] "MTW-TLL-037"       "MTW-TLL013Batch3"  "SMB-ANK-003"      
# [34] "SMB-ANK-005"       "SMB-ANK-011"       "SMB-ANK-013"      
# [37] "SMB-ANK-015"       "SMB-ANK-027Batch1" "SMB-ANK-028"      
# [40] "SMB-ANK027Batch3"  "SMB-HPM006"        "SMB-HPM018"       
# [43] "SMB-PDT007"        "SMB-PTB028"        "SMB-RIN003"       
# [46] "SMB-RIN009"        "SMB-RIN016"        "SMB-WNG-003"      
# [49] "SMB-WNG-004"       "SMB-WNG-007"       "SMB-WNG-018"      
# [52] "SMB-WNG-022"       "SMB-WNG-023"       "SMB-WNG021Batch3" 

# barplot of library sizes
pdf(paste0(outputdir,"BarplotSeqDepth_counts_AllReadsIndonesian.pdf"), width=15)
ggbarplot(meta(allReads_Counts_physeq), "SampleName", "SeqDepth", fill = "SamplePop") + rotate_x_text()
dev.off()

# barplot of library sizes in RPM
SeqDepth_RPM = colSums(otu_table(allReads_RPM_physeq))
sample_data(allReads_RPM_physeq)$SeqDepth = SeqDepth_RPM

pdf(paste0(outputdir,"BarplotSeqDepth_RPM_AllReadsIndonesian.pdf"), width=15)
ggbarplot(meta(allReads_RPM_physeq), "SampleName", "SeqDepth", fill = "SamplePop") + rotate_x_text()
dev.off()

##############################
### Pathogen Investigation ###
##############################

for (rank in c("Phylum","Class","Order","Family")) {
	p=plot_bar(allReads_RPM_physeq, fill = rank)
	fig <- p + geom_bar(stat="identity", position="stack")
	pdf(paste0(outputdir,"AllReads_",rank,".pdf"), width=15)
	fig
	dev.off()
}


