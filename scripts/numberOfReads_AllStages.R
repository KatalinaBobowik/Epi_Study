# Script for creating barplot of number of reads at all different processing stages
# Code developed by Katalina Bobowik, 10.12.2018

# load library and colour palette
library(viridis)
library(Rcmdr)
library(reshape2)
library(ggplot2)

# set working directory
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/QC/"
outputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/QC/"

# set ggplot colour theme to white
theme_set(theme_bw())

# read in summary file and tidy up 
a=read.table(paste0(inputdir,"allReads_AllStages_Ordered.txt"), header=T, sep="\t")
preSubSample=a[,1:4]
preSubSample$Raw=preSubSample$Raw-preSubSample$Trimmed
preSubSample$Trimmed=preSubSample$Trimmed-preSubSample$Unmapped
data=melt(preSubSample)

# Make a barplot of the unmapped reads
pdf(paste0(outputdir,"LibrarySize_UnmappedReads.pdf"), width=10)
ggplot(preSubSample, aes(x=Samples, y=Unmapped)) +
  geom_bar(stat="identity", fill = "#20A387FF") +
  theme(axis.text.x = element_text(angle = 90)) + ylab("Library Size") +
  ggtitle("Unmapped Reads")
dev.off()

pdf(paste0(outputdir,"QCofReads_PreSubSample.pdf"), width=10)
ggplot(data, aes(fill=variable, y=value, x=Samples)) + 
    geom_bar(position="fill", stat="identity") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.title = element_blank()) + 
    scale_fill_viridis(direction = -1, discrete = TRUE, alpha = 0.9) +
    ylab("Number of Reads")
dev.off()

postSubSample=a[,c(1,5,6)]
postSubSample$Repeats_noHuman=NA
postSubSample$Repeats_noHuman = 250000 - (postSubSample$Human + postSubSample$noRepeats_noHuman)
# reorder columns
col_order <- c("Samples", "Human", "Repeats_noHuman", "noRepeats_noHuman")
postSubSample <- postSubSample[, col_order]
data=melt(postSubSample)

pdf(paste0(outputdir,"QCofReads_PostSubSample.pdf"), width=10)
ggplot(data, aes(fill=forcats::fct_rev(variable), y=value, x=Samples)) + 
    geom_bar(position="stack", stat="identity") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.title = element_blank()) + 
    scale_fill_viridis(direction = 1, discrete = TRUE, alpha = 0.9) +
    ylab("Number of Reads")
dev.off()

# make a barplot of the libraries before and after sub-sampling
allreads=read.table(paste0(inputdir,"noRepeats_Sorted.txt"), header=F)
colnames(allreads)=c("Samples","LibrarySize")
pdf(paste0(outputdir,"LibrarySize_AllReads.pdf"), width=10)
ggplot(allreads, aes(x=Samples, y=LibrarySize)) +
  geom_bar(stat="identity", fill = "#20A387FF") +
  theme(axis.text.x = element_text(angle = 90)) + ylab("Library Size") +
  ggtitle("All Reads")
dev.off()

pdf(paste0(outputdir,"LibrarySize_SubSampled.pdf"), width=10)
ggplot(postSubSample, aes(x=Samples, y=noRepeats_noHuman)) +
  geom_bar(stat="identity", fill = "#20A387FF") +
  theme(axis.text.x = element_text(angle = 90)) + ylab("Library Size") +
  ggtitle("Sub-Sampled")
dev.off()




